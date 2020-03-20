#pragma once

#include "common.hpp"
#include "gto.hpp"
#include "system.hpp"
#include "Eigen/Core"
#include <vector>

namespace MOSolver {

struct CGTOs {
    typedef std::vector<ContractedGTO> CGTO_container;

    const ContractedGTO& operator[](size_t idx) const {
        return cgtos_[idx];
    }
    ContractedGTO& operator[](size_t idx) {
        return cgtos_[idx];
    }

    int push_back(ContractedGTO& cgto)
    {
        this->cgtos_.push_back(cgto);
        return 0;
    }
    size_t size(void) const
    {
        return cgtos_.size();
    }
    CGTO_container::iterator begin() {
        return cgtos_.begin();
    }
    CGTO_container::iterator end() {
        return cgtos_.end();
    }

    void add_orbitals( const int l, const Vector3Real center, 
        const std::vector<REAL> &exponent_list, const std::vector<REAL> &coefficient_list);
    
    CGTO_container cgtos_;
};

struct CGTOs
generate_bfs(const System &system, const std::string &basisset_filename);

MatrixXReal
calculate_S(const CGTOs &cgtos);

MatrixXReal
calculate_T(const CGTOs &cgtos);

// Nuclear Attraction Integral
MatrixXReal
calculate_K(const CGTOs &cgtos, const System &atoms);

// two-electron integral
void
calculate_G(const CGTOs &bfs, const MatrixXReal& D, MatrixXReal &G_out);

// two-electron integral for UFH
void
calculate_G_uhf(const CGTOs &bfs, const MatrixXReal& D_alpha, const MatrixXReal &D_beta, 
        MatrixXReal& G_alpha_out, MatrixXReal& G_beta_out);

}   // namespace MOSolver