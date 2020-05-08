#include "common.hpp"
#include "default.hpp"
#include "math_helper.hpp"
#include "gto.hpp"
#include "gto_eval.hpp"
#include "orbitals.hpp"
#include "system.hpp"
#include "Eigen/Core"

namespace MOSolver {

MatrixXReal
symmetric_orthogonalization(const MatrixXReal& S);

MatrixXReal
canonical_orthogonalization(const MatrixXReal& S);

REAL 
check_scf_convergence(const MatrixXReal& D, const MatrixXReal& D_prev, REAL *maxdp);

REAL
calculate_rmsdp(const MatrixXReal &D, const MatrixXReal D_prev);

REAL
calculate_rmsdp(const MatrixXReal &D_diff);

REAL
calculate_maxdp(const MatrixXReal &D_diff);

MatrixXReal
form_D(const MatrixXReal& C, int n_occ_orbitals);

REAL
calculate_E0(const MatrixXReal &D, const MatrixXReal &Hcore, const MatrixXReal &F);

MatrixXReal
initial_guess(const CGTOs& bfs, const System &atoms);

REAL
rhf(const CGTOs& bfs, const System &system,
        const int  nconvergence  = DefaultNConvergence,
        const size_t max_iteration = DefaultMaxIteration );

REAL
uhf(const CGTOs& bfs, const System &system,
        const int  nconvergence  = DefaultNConvergence,
        const size_t max_iteration = DefaultMaxIteration );
}   // namespace MOSolver
