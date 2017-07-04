#pragma once

#include "common.hpp"
#include "gto.hpp"
#include "system.hpp"
#include "Eigen/Core"
#include <vector>

struct CGTOs {
    typedef std::vector<ContractedGTO> CGTO_container;

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
    
    CGTO_container cgtos_;
};


MatrixXReal
calculate_S(CGTOs &cgtos);

MatrixXReal
calculate_T(CGTOs &cgtos);

// Nuclear Attraction Integral
MatrixXReal
calculate_K(CGTOs &cgtos, System &atoms);

// two-electron integral
MatrixXReal
calculate_G(CGTOs &bfs, MatrixXReal& D);

// two-electron integral for UFH
void
calculate_G_uhf(CGTOs &bfs, MatrixXReal& D_alpha, MatrixXReal &D_beta, 
        MatrixXReal& G_alpha_out, MatrixXReal& G_beta_out);
