
#include "gto.hpp"
#include "gto_eval.hpp"
#include "orbitals.hpp"

MatrixXReal
calculate_S(CGTOs &cgtos)
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
calculate_T(CGTOs &cgtos)
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
calculate_K(CGTOs &bfs, System &atoms)
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

MatrixXReal
calculate_G(CGTOs &bfs, MatrixXReal& D)
{
    int dim = bfs.size();
    // TODO Optimize and reduce the loop
    MatrixXReal G = MatrixXReal::Zero(dim, dim);
    for(size_t u = 0; u < dim; u++) {
        for(size_t v = 0; v < dim; v++) {
            REAL temp = 0.;
            for(size_t p = 0; p < dim; p++) {
                for(size_t q = 0; q < dim; q++) {
                    REAL doubleJ = electron_repulsion_CGTO(bfs[u], bfs[v], bfs[p], bfs[q]);
                    REAL k = 0.5 * electron_repulsion_CGTO(bfs[u], bfs[q], bfs[p], bfs[v]);
                    temp += D(p,q) * (doubleJ - k);
                }
            }
            G(u,v) = temp;
        }
    }
    return G;
}   

void
calculate_G_uhf(CGTOs &bfs, MatrixXReal& D_alpha, MatrixXReal &D_beta, 
        MatrixXReal& G_alpha_out, MatrixXReal& G_beta_out)
{
    int dim = bfs.size();
    MatrixXReal D_total = D_alpha + D_beta;
    for(size_t u = 0; u < dim; u++) {
        for(size_t v = 0; v < dim; v++) {
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
            G_beta_out(u,v)  = temp_beta;
        }
    }
    return;
}
