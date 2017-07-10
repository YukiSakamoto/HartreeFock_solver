#pragma once

#include "common.hpp"
#include "Eigen/Core"
#include <vector>


struct AtomType {
    char atom[3];
    int  Z;
};

const AtomType atom_table[] = {
    { "H",  1 },
    { "He", 2 },
    { "Li", 3 },
    { "Be", 4 },
    { "B",  5 },
    { "C",  6 },
    { "N",  7 },
    { "O",  8 },
    { "F",  9 },
    { "Ne", 10},
};


struct Atom {
    Atom(int atomic_number, REAL x, REAL y, REAL z):
        atomic_number(atomic_number), center(x,y,z)
    {
        if (atomic_number <= 0) {   throw;  }
    }
    Atom(int atomic_number, Vector3Real pos):
        atomic_number(atomic_number), center(pos)
    {
        if (atomic_number <= 0) {   throw;  }
    }
    int atomic_number;
    Vector3Real center;
};



struct System {
    typedef std::vector<Atom> AtomContainer;
    System(int charge = 0, int nspin = 0) :
        charge_(charge), nspin_(nspin)
    {;}
    void set_charge(int charge)
    {
        charge_ = charge;
    }
    void set_nspin(int nspin)
    {
        nspin_ = nspin;
    }
    int nspin(void) const 
    {
        return nspin_;
    }
    int total_electrons(void) const 
    {
        int tot_e = 0;
        for(size_t i = 0; i < atom_list_.size(); i++) {
            tot_e += atom_list_[i].atomic_number;
        }
        tot_e += (-1)*charge_;  // charge_ < 0 means that the number of total electrons is bigger than the that of neutral

        return tot_e;
    }
    
    Atom& operator[](size_t idx) {
        return atom_list_[idx];
    }
    const Atom& operator[](size_t idx) const {
        return atom_list_[idx];
    }
    size_t size(void) const
    {
        return atom_list_.size();
    }
    void add_atom(int atomic_number, REAL x, REAL y, REAL z) {
        atom_list_.push_back(Atom(atomic_number, x, y, z));
    }
    void add_atom(int atomic_number, Vector3Real pos) {
        atom_list_.push_back(Atom(atomic_number, pos));
    }
    void add_atom(const Atom &atom) {
        atom_list_.push_back(atom);
    }
    void push_back(const Atom &atom) {
        atom_list_.push_back(atom);
    }
    void generate_basis();
    REAL nuclear_repulsion() const;
    std::vector<Atom> atom_list_;
    int charge_;
    int nspin_;
};




