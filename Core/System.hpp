#include "Definition.hpp"
#include <vector>
#include <numeric>

namespace MOSolver2 {


struct System {
    System(int charge = 0, int nspin = 0) : charge(charge), nspin(nspin)
    {;}

    void add_atom(int Z, REAL pos_x, REAL pos_y, REAL pos_z)
    {
        this->atomic_number.push_back(Z);
        this->x.push_back(pos_x);
        this->y.push_back(pos_y);
        this->z.push_back(pos_z);
    }
    int size() const
    {   return this->atomic_number.size(); }

    int n_electron() const
    {   
        return std::accumulate(atomic_number.begin(), atomic_number.end(), 0) - charge;   
    }

    void get_nth_atom(const int n, int &atomic_number_out, REAL &x_out, REAL &y_out, REAL &z_out)
    {
        atomic_number_out = this->atomic_number[n];
        x_out = this->x[n];
        y_out = this->y[n];
        z_out = this->z[n];
    }

    int charge;
    int nspin;
    std::vector<REAL> x;
    std::vector<REAL> y;
    std::vector<REAL> z;
    std::vector<int>  atomic_number;
};
    
} //namespace MOSolver 
