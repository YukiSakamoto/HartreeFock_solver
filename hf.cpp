#include "common.hpp"

#include "math_helper.hpp"
#include "gto.hpp"
#include "gto_eval.hpp"
#include "orbitals.hpp"
#include "system.hpp"

#include <iostream>
#include <Eigen/Eigen> 
#include <Eigen/Core>

MatrixXReal
symmetric_orthogonalization(const MatrixXReal& S)
{
    // Overlap Integral Matrix (S) should be the 
    //  self-adjoint matrix by the definition.
    Eigen::SelfAdjointEigenSolver<MatrixXReal> es(S);
    if (es.info() != Eigen::Success) {  throw;  }

    size_t row = S.rows();
    size_t col = S.cols();
    MatrixXReal l_rt = MatrixXReal::Zero(row,col);

    VectorXReal l = es.eigenvalues();
    MatrixXReal U = es.eigenvectors();
    
    int len = es.eigenvalues().size();
    for(size_t i = 0; i < len; i++) {
        l_rt(i,i) = 1./std::sqrt(l[i]);
    }
    MatrixXReal X = U * (l_rt * (U.conjugate()));

    return X;
}

MatrixXReal
canonical_orthogonalization(const MatrixXReal& S) 
{
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

    int row = D_diff.rows();
    int col = D_diff.cols();
    for(size_t i = 0; i < D_diff.rows(); i++) {
        for(size_t j = 0; j < D_diff.cols(); j++) {
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

REAL
calculate_rmsdp(const MatrixXReal &D_diff)
{
    REAL norm2 = D_diff.squaredNorm();
    size_t n_elements = D_diff.rows() *D_diff.cols();
    REAL rmsdp( std::sqrt(norm2/n_elements) );
    return rmsdp;
}


REAL
calculate_maxdp(const MatrixXReal &D_diff)
{
    REAL retval = 0.;
    size_t rows = D_diff.rows();
    size_t cols = D_diff.cols();
    for(size_t i = 0; i < rows; i++) {
        for(size_t j = 0; j < cols; j++) {
            retval = std::max(retval, std::abs(D_diff(rows, cols)));
        }
    }
    return retval;
}

MatrixXReal
initial_guess(const CGTOs& bfs, const System &atoms)
{
    // XXX For Now, Returns the ZERO Maxrix
    return MatrixXReal::Zero(bfs.size(), bfs.size());
}

MatrixXReal
form_D(const MatrixXReal& C, int n_occ_orbitals)
{
    size_t row = C.rows();
    size_t col = C.cols();
    MatrixXReal D = MatrixXReal::Zero(row, col);
    for(int u = 0 ; u < row; u++) {
        for(int v = 0; v < col; v++) {
            for(int a = 0 ; a < n_occ_orbitals; a++) {
                D(u,v) += 2.0 * C(u,a) * C(v,a);
            }
        }
    }
    return D;
}

REAL
calculate_E0(const MatrixXReal &D, const MatrixXReal &Hcore, const MatrixXReal &F)
{
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
rhf(CGTOs& bfs, System &system)
{
    int tot_electrons = system.total_electrons();
    int n_occ_orbitals = tot_electrons / 2;
    // Check the System
    {
        if (system.total_electrons() % 2 != 0 || system.nspin() != 1) {
            std::cerr << 
                "ERROR: For RHF calculation, the number of electrons must be even, and nspin must be 1" << std::endl;
            throw;
        }
    }
     std::cout << "************************************************************\n";
     std::cout << "  Building Matrices... (S, Hcore, init_D)\n";
     std::cout << "************************************************************\n";
    //==================================================
    //  Building Integral Matrices
    //==================================================
    // 1. Overlap Matrix
    MatrixXReal S = calculate_S(bfs);
    std::cout << "Overlap Matrix:\n" << S << std::endl;

    // 2. Orthodiagonalization of S-matrix (|X>)
    //MatrixXReal X = symmetric_orthogonalization(S);
    //std::cout << "symmetric_orthogonalization:\n" << X << std::endl;
    MatrixXReal X = canonical_orthogonalization(S);
    std::cout << "canonical_orthogonalization:\n" << X << std::endl;
    // (<X|)   => <X|S|X> equals E(Unit Matrix)
    MatrixXReal X_adj = X.adjoint();

    // 3. Kinetic Energy Integral Matrix
    MatrixXReal T = calculate_T(bfs);
    std::cout << "Kinetic Matrix:\n" << T << std::endl;

    // 4. Nuclear Attraction Integral Matrix
    MatrixXReal V = calculate_K(bfs, system);
    std::cout << "V(NuclearAttractionIntegral):\n" << V << std::endl;

    // 5. Hcore Matrix
    MatrixXReal Hcore = T + V;
    std::cout << "Hcore = T + V:\n" << Hcore << std::endl;

    // 6. Obtain Initial Density Matrix from Guess
    MatrixXReal D = initial_guess(bfs, system);
    std::cout << "D(initial guess):\n" << D << std::endl;

    // 7. Nuclear Repulsions;
    REAL NEI = system.nuclear_repulsion();
    std::cout << "NEI: " << NEI << std::endl;

     std::cout << "************************************************************\n";
     std::cout << "  Entering the SCF Loop\n";
     std::cout << "************************************************************\n";
    //==================================================
    //  Entering the SCF Loop
    //==================================================
    const int max_iteration = 20;
    bool convergence_flag = false;
    REAL E_conv = 0.;
    for(size_t i = 0; i < max_iteration; i++) {
        std::cout << " [Iteration " << i << " ] " << std::endl;
        std::cout << "D:\n" << D << std::endl;

        // Calculate the multicenter integrals with D;
        MatrixXReal G = calculate_G(bfs, D);
        std::cout << "G:\n" <<G << std::endl;

        // Build current Fock matrix
        MatrixXReal F = Hcore + G;
        std::cout << "F:\n" << F << std::endl;

        // Rotating F matrix 
        //   (F' = <X|F|X>, where <X|S|X> = E)
        MatrixXReal F_prim = X_adj * F * X;
        std::cout << "F':\n" << F_prim << std::endl;

        // Solve the Fock matrix.
        Eigen::SelfAdjointEigenSolver<MatrixXReal> es(F_prim);
        if (es.info() != Eigen::Success) {  throw;  }
        
        // Energies and New coefficients 
        VectorXReal e = es.eigenvalues();
        MatrixXReal C_new_prime = es.eigenvectors();
        std::cout << "e: \n" << e << std::endl;
        std::cout << "C_new':\n" << C_new_prime << std::endl;
        MatrixXReal C_new = X * C_new_prime;
        std::cout << "C_new:\n" << C_new << std::endl;

        // Calculate the NEW Density matrix ( <C|C> )
        //MatrixXReal D_new =  C_new_prime * C_new_prime.transpose();
        MatrixXReal D_new =  form_D(C_new, n_occ_orbitals);
        std::cout << "D_new:\n" << D_new << std::endl;

        // Calculate the Hartree Fock Energy
        REAL E0 = calculate_E0(D,Hcore,F);
        REAL Etot = E0 + NEI;

        // Check the convergence
        REAL rmsdp = 0.;    
        REAL maxdp = 0.;
        rmsdp = check_scf_convergence(D_new, D, &maxdp);
        int nconv =  5;
        REAL convergence = std::pow(10.0, -nconv);
        if (rmsdp < convergence) {
            convergence_flag = true;
        }

        // Output current cycle
        std::cout << "RMSDP: "<< rmsdp << ", MAXDP: " << maxdp << std::endl;
        std::cout << "E0: " << E0  << "  Etot: " << Etot << std::endl;
        std::cout << "============================================================\n";

        if (convergence_flag == true) {
            std::cout << "CONVERGENCE ACHIEVED\n";
            E_conv = Etot;
            break;
        } else {
            //D = D * 0.3 + D_new * 0.7;
            D = D_new;
            std::cout << std::endl;
        }
    }
    return E_conv;
}


