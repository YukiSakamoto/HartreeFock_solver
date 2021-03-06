#include "system.hpp"

namespace MOSolver {

REAL 
System::nuclear_repulsion() const
{
    REAL t = 0.;
    for(size_t i = 0; i < this->size(); i++) {
        for(size_t j = i + 1; j < this->size(); j++) {
            if (i != j) {
                REAL r = (atom_list_[i].center - atom_list_[j].center).norm();
                t += (atom_list_[i].atomic_number * atom_list_[j].atomic_number) / r;
            }
        }
    }
    return t;
}

} // namespace MOSolver 