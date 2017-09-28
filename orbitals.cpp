#include "common.hpp"
#include "gto.hpp"
#include "gto_eval.hpp"
#include "orbitals.hpp"
#include "basis_set.hpp"

struct CGTOs
generate_bfs(const System &system, const std::string &basisset_filename)
{
    CGTOs bfs;
    std::vector<Atom>::const_iterator it = system.atom_list_.begin();
    struct BasisSet dat( parse_basisset_file(basisset_filename) );
    for(; it != system.atom_list_.end(); it++) {
        struct AtomBasis atom_basis = dat.get(it->atomic_number);
        std::cout << atom_basis.to_str() << std::endl;
        for(size_t i_shell = 0; i_shell < atom_basis.orbitals.size(); i_shell++) {
            int l = angular_momentum(atom_basis.orbitals[i_shell].type);
            bfs.add_orbitals(l, it->center, 
                    atom_basis.orbitals[i_shell].exponents,
                    atom_basis.orbitals[i_shell].coeffs);
        }
    }
    return bfs;
}


void
CGTOs::add_orbitals(const int l, const Vector3Real center, 
        const std::vector<REAL> &exponent_list, const std::vector<REAL> &coefficient_list)
{
    if (l < 0) {
        std::cerr << "The value of l (angular momentum) must be ZERO or POSITIVE." << std::endl;
        throw;
    }
    if (MAX_L < l) {
        std::cerr << "The angular momentum " << l << " is not supported" << std::endl;
    }
    if (exponent_list.size() != coefficient_list.size() ) {
        std::cerr << "The length of exponent_list and coefficient_list always must be the same." <<std::endl;
        throw;
    }
    size_t contraction = exponent_list.size();
    std::vector<ContractedGTO> v;
    switch (l) {
        case 0: // s orbital
            v.push_back( ContractedGTO(0, 0, 0, center) );
            break;
        case 1:
            v.push_back( ContractedGTO(1, 0, 0, center) );  // px
            v.push_back( ContractedGTO(0, 1, 0, center) );  // py
            v.push_back( ContractedGTO(0, 0, 1, center) );  // pz
            break;
        case 2:
            v.push_back( ContractedGTO(2, 0, 0, center) );  // dx2
            v.push_back( ContractedGTO(0, 2, 0, center) );  // dy2
            v.push_back( ContractedGTO(0, 0, 2, center) );  // dz2
            v.push_back( ContractedGTO(1, 1, 0, center) );  // dxy
            v.push_back( ContractedGTO(0, 1, 1, center) );  // dyz
            v.push_back( ContractedGTO(1, 0, 1, center) );  // dzx
            break;
        default:
            // Never get here
            throw;
            break;
    }
    for(size_t i = 0; i < v.size() ; i++) {
        for(size_t j = 0; j < contraction; j++) {
            v[i].add_primitiveGTO(coefficient_list[j], exponent_list[j]);
        }
        v[i].normalize();
        this->cgtos_.push_back(v[i]);
    }
    return;
}

MatrixXReal
calculate_S(const CGTOs &cgtos)
{
    size_t len = cgtos.size();
    MatrixXReal S = MatrixXReal::Zero(len, len);
    for(size_t i = 0; i < len; i++) {
        for(size_t j = 0; j < len; j++) {
            S(i,j) = overlap_CGTO(cgtos[i], cgtos[j]);
        }
    }
    return S;
}

MatrixXReal
calculate_T(const CGTOs &cgtos)
{
    size_t len = cgtos.size();
    MatrixXReal T = MatrixXReal::Zero(len, len);
    for(size_t i = 0; i < len; i++) {
        for(size_t j = 0; j < len; j++) {
            T(i,j) = kinetic_CGTO(cgtos[i], cgtos[j]);
        }
    }
    return T;
}

MatrixXReal
calculate_K(const CGTOs &bfs, const System &atoms)
{
    size_t len = bfs.size();
    MatrixXReal T = MatrixXReal::Zero(len, len);
    for(size_t atom_idx = 0; atom_idx < atoms.size(); atom_idx++) {
        MatrixXReal v = MatrixXReal::Zero(len, len);
        for (size_t i = 0; i < len; i++) {
            for(size_t j = 0; j < len; j++) {
                v(i,j) = nuclear_attraction_CGTO(bfs[i], bfs[j], atoms[atom_idx]);
            }
        }
        T = T + v;
    }
    return T;
}

void
calculate_G(const CGTOs &bfs, const MatrixXReal& D, MatrixXReal &G_out)
{
    size_t dim = bfs.size();
    // TODO Optimize and reduce the loop
    for(size_t u = 0; u < dim; u++) {
        for(size_t v = 0; v < dim; v++) {
            if (u < v) { continue;  }
            REAL temp = 0.;
            for(size_t p = 0; p < dim; p++) {
                for(size_t q = 0; q < dim; q++) {
                    REAL doubleJ = electron_repulsion_CGTO(bfs[u], bfs[v], bfs[p], bfs[q]);
                    REAL k = 0.5 * electron_repulsion_CGTO(bfs[u], bfs[q], bfs[p], bfs[v]);
                    temp += D(p,q) * (doubleJ - k);
                }
            }
            G_out(u,v) = temp;
            G_out(v,u) = temp;
        }
    }
}   

void
calculate_G_uhf(const CGTOs &bfs, const MatrixXReal& D_alpha, const MatrixXReal &D_beta, 
        MatrixXReal& G_alpha_out, MatrixXReal& G_beta_out)
{
    size_t dim = bfs.size();
    MatrixXReal D_total = D_alpha + D_beta;
    for(size_t u = 0; u < dim; u++) {
        for(size_t v = u; v < dim; v++) {
            REAL temp_alpha = 0.;
            REAL temp_beta  = 0.;
            for(size_t p = 0; p < dim; p++) {
                for(size_t q = 0; q < dim; q++) {
                    REAL J = electron_repulsion_CGTO(bfs[u], bfs[v], bfs[p], bfs[q]);
                    REAL K = electron_repulsion_CGTO(bfs[u], bfs[q], bfs[p], bfs[v]);
                    temp_alpha += D_total(q,p) * J - D_alpha(q,p) * K;
                    temp_beta  += D_total(q,p) * J - D_beta(q,p) * K;
                }
            }
            G_alpha_out(u,v) = temp_alpha;
            G_alpha_out(v,u) = temp_alpha;
            if (u != v) {
                G_beta_out(u,v)  = temp_beta;
                G_beta_out(v,u)  = temp_beta;
            }
        }
    }
    return;
}
