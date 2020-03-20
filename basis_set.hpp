#pragma once

#include "common.hpp"
#include "boost/format.hpp"

namespace MOSolver {

inline
std::string vec_real_to_string(const std::vector<REAL> &v) {
    std::stringstream ss;
    ss << "[ ";
    for(std::vector<REAL>::const_iterator it = v.begin(); it != v.end(); it++) {
        ss << boost::format("%12.8f ") % *it;
    }
    ss << " ]";
    return ss.str();
}


struct Shell {
    char type;  // S,P,D,F
    std::vector<REAL> coeffs;
    std::vector<REAL> exponents;
    //constructor
    Shell(char type) :type(type) {;}

    std::string to_str(int indentlevel = 1) const 
    {
        std::stringstream ss;
        std::string indent;
        for(int i = 0; i < indentlevel; i++) {  indent.push_back('\t'); }
        ss << indent << boost::format("%1s  Exp:   %s") % type % vec_real_to_string(exponents) << std::endl;
        ss << indent << boost::format("%1s  Coeff: %s") % type % vec_real_to_string(coeffs) << std::endl;
        return ss.str();
    }
};

struct AtomBasis {
    int atomic_number;  // Atom name is upto 2 characters currently.
    std::vector<Shell> orbitals;

    std::string to_str(int indentlevel = 0) const 
    {
        std::stringstream ss;
        ss << boost::format("%2s (%3d)") % Elements[atomic_number] % atomic_number << std::endl;
        for(std::vector<Shell>::const_iterator it = orbitals.begin(); it != orbitals.end(); it++) {
            ss << it->to_str(indentlevel + 1);
        }
        return ss.str();
    }
};

struct BasisSet {
    typedef std::vector<AtomBasis> basis_set_container_type;
    
    std::string name;
    std::vector<AtomBasis> container;

    //constructor
    BasisSet(std::string name) : name(name){;}
    void append(const AtomBasis &b) {
        this->container.push_back(b);
    }
    const struct AtomBasis &get(int atomic_number) {
        for(size_t i = 0; i < container.size(); i++) {
            if (container[i].atomic_number == atomic_number) {
                return container[i];
            }
        }
        // Never get here
        throw;
    }
};

struct BasisSet 
parse_basisset_file(const std::string &filename);

} //namespace MOSolver 