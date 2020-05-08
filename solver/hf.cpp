#include "common.hpp"
#include "math_helper.hpp"
#include "gto.hpp"
#include "gto_eval.hpp"
#include "orbitals.hpp"
#include "system.hpp"
#include "default.hpp"

#include <Eigen/Eigen> 
#include <Eigen/Core>   // for Solver

#include <boost/format.hpp>

#include <iostream>

namespace MOSolver {

MatrixXReal
symmetric_orthogonalization(const MatrixXReal& S)
{
    // Szabo. pp.143 (3.167)
    // Overlap Integral Matrix (S) should be the 
    //  self-adjoint matrix by the definition.
    Eigen::SelfAdjointEigenSolver<MatrixXReal> es(S);
    if (es.info() != Eigen::Success) {  throw;  }
    VectorXReal l = es.eigenvalues();
    MatrixXReal U = es.eigenvectors();
    MatrixXReal l_rt = MatrixXReal::Zero(S.rows(), S.cols() );
    
    size_t len = es.eigenvalues().size();
    for(size_t i = 0; i < len; i++) {
        l_rt(i,i) = 1./std::sqrt(l[i]);
    }
    MatrixXReal X = U * (l_rt * (U.conjugate()));

    return X;
}

MatrixXReal
canonical_orthogonalization(const MatrixXReal& S) 
{
    // Szabo. pp.144 (3.169 - 3.172)
    // Overlap Integral Matrix (S) should be the 
    //  self-adjoint matrix by the definition.
    Eigen::SelfAdjointEigenSolver<MatrixXReal> es(S);
    if (es.info() != Eigen::Success) {  throw;  }

    VectorXReal l = es.eigenvalues();
    MatrixXReal U = es.eigenvectors();

    size_t row = S.rows();
    size_t col = S.cols();
    MatrixXReal X = MatrixXReal::Zero(row, col);
    for (size_t i = 0; i < row; i++) {
        for(size_t j = 0; j < col; j++) {
            X(i,j) = U(i,j) / std::sqrt(l[j]);
        }
    }
    return X;
}


REAL 
check_scf_convergence(const MatrixXReal& D, const MatrixXReal& D_prev, REAL *maxdp)
{
    if (D.rows() != D_prev.rows() || D.cols() != D_prev.cols()) {
        throw;
    }
    MatrixXReal D_diff = D - D_prev;
    REAL maxdp_acc = 0.;
    REAL rmsdp_acc = 0.;

    size_t row = D_diff.rows();
    size_t col = D_diff.cols();
    for(size_t i = 0; i < row; i++) {
        for(size_t j = 0; j < col; j++) {
            rmsdp_acc += std::pow(D_diff(i,j), 2);
            maxdp_acc = std::max(maxdp_acc, std::abs(D_diff(i,j)) );
        }
    }
    if (maxdp != NULL) {
        *maxdp = maxdp_acc;
    }
    REAL rmsdp = std::sqrt( rmsdp_acc/(row*col) );
    return rmsdp;
}

REAL
calculate_rmsdp(const MatrixXReal &D, const MatrixXReal D_prev)
{
    MatrixXReal D_diff(D - D_prev);
    REAL norm2 = D_diff.squaredNorm();
    size_t n_elements = D_diff.rows() *D_diff.cols();
    REAL rmsdp( std::sqrt(norm2/n_elements) );
    return rmsdp;
}

MatrixXReal
initial_guess(const CGTOs& bfs, const System /*&atoms*/)
{
    // XXX For Now, Returns the ZERO Maxrix
    return MatrixXReal::Zero(bfs.size(), bfs.size());
}

MatrixXReal
form_D(const MatrixXReal& C, size_t n_occ_orbitals)
{
    // Szabo. pp. 139 (3.145)
    size_t row = C.rows();
    size_t col = C.cols();
    MatrixXReal D = MatrixXReal::Zero(row, col);
    for(size_t u = 0 ; u < row; u++) {
        for(size_t v = 0; v < col; v++) {
            for(size_t a = 0 ; a < n_occ_orbitals; a++) {
                D(u,v) += 2.0 * C(u,a) * C(v,a);
            }
        }
    }
    return D;
}

MatrixXReal
form_D_uhf(const MatrixXReal &C_spin, size_t n_spin_electron)
{
    // Szabo. pp.213 (3.342 and 3.343)
    size_t row = C_spin.rows();
    size_t col = C_spin.cols();
    MatrixXReal D_spin_new = MatrixXReal::Zero(row, col);
    for(size_t u = 0; u < row; u++) {
        for(size_t v = 0; v < col; v++) {
            for(size_t a = 0; a < n_spin_electron; a++) {
                D_spin_new(u,v) += C_spin(u,a) * C_spin(v,a);
            }
        }
    }
    return D_spin_new;
}

REAL
calculate_E0(const MatrixXReal &D, const MatrixXReal &Hcore, const MatrixXReal &F)
{
    // Szabo. pp.150 (3.184): 
    REAL E0 = 0.;
    MatrixXReal H_F = Hcore + F;
    size_t row = H_F.rows();
    size_t col = H_F.cols();
    for(size_t u = 0; u < row; u++) {
        for(size_t v = 0; v < col; v++) {
            E0 += D(v,u) * H_F(u,v);
        }
    }
    E0 *= 0.5;
    return E0;
}

REAL
calculate_E0_uhf(const MatrixXReal &D_alpha, const MatrixXReal &D_beta, const MatrixXReal &Hcore, 
        const MatrixXReal &F_alpha, const MatrixXReal &F_beta)
{
    // Szabo. pp. 215 (Exercise 3.40)
    REAL E0 = 0.;
    size_t row = Hcore.rows();
    size_t col = Hcore.cols();
    for(size_t u = 0; u < row; u++) {
        for(size_t v = 0; v < col; v++) {
            REAL da = D_alpha(v,u);
            REAL db = D_beta(v,u);
            E0 += ((da+db) * Hcore(u,v) + da * F_alpha(u,v) + db * F_beta(u,v));
        }
    }
    E0 *= 0.5;
    return E0;
}

REAL
rhf(const CGTOs& bfs, const System &system,
        const int  nconvergence, 
        const size_t max_iteration)
{
    std::cout << nconvergence << std::endl;
    int tot_electrons = system.total_electrons();
    int n_occ_orbitals = tot_electrons / 2;
    int dim = bfs.size();
    // Check the System
    {
        if (system.total_electrons() % 2 != 0 || system.nspin() != 0) {
            std::cerr << 
                "ERROR: For RHF calculation, the number of electrons must be even, and nspin must be 1" << std::endl;
            throw;
        }
    }
    //==================================================
    //  Building Integral Matrices
    //==================================================
    // 1. Overlap Matrix
    MatrixXReal S = calculate_S(bfs);

    // 2. Orthodiagonalization of S-matrix (|X>)
    //MatrixXReal X = symmetric_orthogonalization(S);
    MatrixXReal X = canonical_orthogonalization(S);
    // (<X|)   => <X|S|X> equals E(Unit Matrix)
    MatrixXReal X_adj = X.adjoint();

    // 3. Kinetic Energy Integral Matrix
    MatrixXReal T = calculate_T(bfs);

    // 4. Nuclear Attraction Integral Matrix
    MatrixXReal V = calculate_K(bfs, system);

    // 5. Hcore Matrix
    MatrixXReal Hcore = T + V;

    // 6. Obtain Initial Density Matrix from Guess
    MatrixXReal D = initial_guess(bfs, system);

    // 7. Nuclear Repulsions;
    REAL NEI = system.nuclear_repulsion();
    std::cout << "NEI: " << NEI << std::endl;

     std::cout << "************************************************************\n";
     std::cout << "  Entering the SCF Loop\n";
     std::cout << "************************************************************\n";
    //==================================================
    //  Entering the SCF Loop
    //==================================================
    bool convergence_flag = false;
    REAL E_conv = 0.;
    for(size_t i = 0; i < max_iteration; i++) {
        std::cout << " [Iteration " << i << " ] " << std::endl;

        // Calculate the multicenter integrals with D;
        MatrixXReal G = MatrixXReal::Zero(dim, dim);
        calculate_G(bfs, D, G);

        // Build current Fock matrix
        MatrixXReal F = Hcore + G;

        // Rotating F matrix 
        //   (F' = <X|F|X>, where <X|S|X> = E)
        MatrixXReal F_prim = X_adj * F * X;

        // Solve the Fock matrix.
        Eigen::SelfAdjointEigenSolver<MatrixXReal> es(F_prim);
        if (es.info() != Eigen::Success) {  throw;  }
        
        // Energies and New coefficients 
        VectorXReal e = es.eigenvalues();
        MatrixXReal C_new_prime = es.eigenvectors();
        MatrixXReal C_new = X * C_new_prime;

        // Calculate the NEW Density matrix ( <C|C> )
        MatrixXReal D_new =  form_D(C_new, n_occ_orbitals);

        // Calculate the Hartree Fock Energy
        REAL E0 = calculate_E0(D,Hcore,F);
        REAL Etot = E0 + NEI;

        // Check the convergence
        REAL rmsdp = 0.;    
        REAL maxdp = 0.;
        rmsdp = check_scf_convergence(D_new, D, &maxdp);
        REAL convergence = std::pow(10.0, -nconvergence);
        REAL maxdp_convergence = std::pow(10.0, -nconvergence+2);
        if (rmsdp < convergence && maxdp < maxdp_convergence) {
            convergence_flag = true;
        }

        // Output current cycle
        std::cout << "RMSDP: "<< rmsdp << ", MAXDP: " << maxdp << std::endl;
        std::cout << boost::format("E0 : %15.10f   Etot: %15.10f\n") % E0 % Etot;
        std::cout << "============================================================\n";

        if (convergence_flag == true) {
            std::cout << "CONVERGENCE ACHIEVED\n";
            E_conv = Etot;
            break;
        } else {
            //D = D * 0.3 + D_new * 0.7;
            D = D_new;
        }
    }
    return E_conv;
}

REAL
uhf(const CGTOs& bfs, const System &system,
        const int  nconvergence, 
        const size_t max_iteration)
{
    int total_electrons = system.total_electrons();
    int n_spin = system.nspin();
    int occ_alpha = (total_electrons + n_spin) / 2;
    int occ_beta  = (total_electrons - n_spin) / 2;
    std::cout << "alpha: " << occ_alpha << std::endl;
    std::cout << "beta: " << occ_beta << std::endl;
    // Check the System
    {
       if (occ_alpha - occ_beta != n_spin) {
           std::cerr << 
               "ERROR: For RHF calculation, the number of electrons must be even, and nspin must be 1" << std::endl;
           throw;
       }
    }

    //==================================================
    //  Building Integral Matrices
    //==================================================
    // 1. Overlap Matrix
    MatrixXReal S = calculate_S(bfs);

    // 2. Orthodiagonalization of S-matrix (|X>)
    //MatrixXReal X = symmetric_orthogonalization(S);
    MatrixXReal X = canonical_orthogonalization(S);
    // (<X|)   => <X|S|X> equals E(Unit Matrix)
    MatrixXReal X_adj = X.adjoint();

    // 3. Kinetic Energy Integral Matrix
    MatrixXReal T = calculate_T(bfs);

    // 4. Nuclear Attraction Integral Matrix
    MatrixXReal V = calculate_K(bfs, system);

    // 5. Hcore Matrix
    MatrixXReal Hcore = T + V;

    // 6. Obtain Initial Density Matrix from Guess
    MatrixXReal D_alpha = initial_guess(bfs, system);
    MatrixXReal D_beta  = initial_guess(bfs, system);
    MatrixXReal D_total = D_alpha + D_beta;

    // 7. Nuclear Repulsions;
    REAL NEI = system.nuclear_repulsion();
    std::cout << "NEI: " << NEI << std::endl;

    std::cout << "************************************************************\n";
    std::cout << "  Entering the SCF Loop\n";
    std::cout << "************************************************************\n";
    bool convergence_flag = false;
    REAL E_conv = 0.;
    for(size_t i = 0; i < max_iteration; i++) {
        std::cout << " [Iteration " << i << " ] " << std::endl;

        // Calculate the multicenter integrals with D;
        int dim = bfs.size();
        MatrixXReal G_alpha = MatrixXReal::Zero(dim, dim);
        MatrixXReal G_beta  = MatrixXReal::Zero(dim, dim);
        calculate_G_uhf(bfs, D_alpha, D_beta, G_alpha, G_beta);

        // Build current Fock matrix
        MatrixXReal F_alpha = Hcore + G_alpha;
        MatrixXReal F_beta  = Hcore + G_beta;

        // Rotating F matrix 
        //   (F' = <X|F|X>, where <X|S|X> = E)
        MatrixXReal F_alpha_prim = X_adj * F_alpha * X;
        MatrixXReal F_beta_prim  = X_adj * F_beta  * X;

        // Solve the Fock matrix.
        Eigen::SelfAdjointEigenSolver<MatrixXReal> es_alpha(F_alpha_prim);
        Eigen::SelfAdjointEigenSolver<MatrixXReal>  es_beta(F_beta_prim);
        if (es_alpha.info() != Eigen::Success || es_beta.info() != Eigen::Success) {  throw;  }

        // Energies and New coefficients 
        VectorXReal e_alpha = es_alpha.eigenvalues();
        VectorXReal e_beta  = es_beta.eigenvalues();
        MatrixXReal C_alpha_new_prime = es_alpha.eigenvectors();
        MatrixXReal C_beta_new_prime  = es_beta.eigenvectors();

        MatrixXReal C_alpha_new = X * C_alpha_new_prime;
        MatrixXReal C_beta_new  = X * C_beta_new_prime;

        // Calculate the NEW Density matrix ( <C|C> )
        MatrixXReal D_alpha_new =  form_D_uhf(C_alpha_new, occ_alpha);
        MatrixXReal D_beta_new  =  form_D_uhf(C_beta_new,  occ_beta);

        // Calculate the Hartree Fock Energy
        REAL E0 = calculate_E0_uhf(D_alpha, D_beta, Hcore, F_alpha, F_beta);
        REAL Etot = E0 + NEI;

        // Check the convergence
        REAL rmsdp_alpha = 0.;    
        REAL maxdp_alpha = 0.;
        REAL rmsdp_beta  = 0.;    
        REAL maxdp_beta  = 0.;
        rmsdp_alpha = check_scf_convergence(D_alpha_new, D_alpha, &maxdp_alpha);
        rmsdp_beta = check_scf_convergence(D_beta_new, D_beta, &maxdp_beta);

        REAL convergence = std::pow(10.0, -nconvergence);
        REAL maxdp_convergence = std::pow(10.0, -nconvergence+2);
        if (rmsdp_alpha + rmsdp_beta < convergence && maxdp_alpha + maxdp_beta < maxdp_convergence) {
            convergence_flag = true;
        }

        std::cout << "RMSDP: "<< rmsdp_alpha + rmsdp_beta << ", MAXDP: " << maxdp_alpha + maxdp_beta << std::endl;
        std::cout << boost::format("E0 : %15.10f   Etot: %15.10f\n") % E0 % Etot;
        std::cout << "============================================================\n";

        if (convergence_flag == true) {
            std::cout << "CONVERGENCE ACHIEVED\n";
            E_conv = Etot;
            break;
        } else {
            //D = D * 0.3 + D_new * 0.7;
            D_alpha = D_alpha_new;
            D_beta = D_beta_new;
            std::cout << std::endl;
        }
    }
    return E_conv;
}




} // namespace MOSolver
