#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

//for boost::tokenizer
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

#include <cstdlib>
#include <algorithm>

#include "common.hpp"
#include "basis_set.hpp"

namespace MOSolver {

//============================================================
//  This file contains the routine for parse basis set data
//   formatted for Gaussian program.
//   The data can be obtained from EMSL basis set exchange.)
//============================================================

std::vector<std::string> split_by_space(const std::string &str)
{
    std::vector<std::string> retval;
    boost::char_separator<char> sep(" \t\n");
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    tokenizer tokens(str,sep);
    for(tokenizer::iterator it = tokens.begin(); it != tokens.end(); it++) {
        retval.push_back(*it);
    }
    return retval;
}


struct BasisSet 
parse_basisset_file(const std::string &filename)
{
    std::ifstream ifs(filename.c_str() ,std::ios::in);
    std::string line_buf;
    struct BasisSet retval(filename);

    if(!ifs.is_open())
    {
        std::cerr << "Error! maybe the file doesn't exist." << std::endl;
        throw;
    }
    
    // Ignore the comments and find first line of basis ste data;
    while(! ifs.eof() ) {
        std::getline(ifs, line_buf);
        if (line_buf == "****") {   break;  }
    }

    while(true) {
        struct AtomBasis new_atom;
        // Parse Atom_Name
        std::getline(ifs, line_buf);
        std::vector<std::string> a = split_by_space(line_buf);
        if (a.size() == 0) {    break;  }
        int atomic_num = get_atomic_number(a[0]);
        new_atom.atomic_number = atomic_num;
        int parse_done = 0;

        while(true) {
            // Parse L, Contraction, Scaling
            std::getline(ifs, line_buf);
            if (line_buf == "****") {   
                retval.append(new_atom);
                // Next Atom
                break;
            }
            a = split_by_space(line_buf);
            if (a.size() < 2) {
                throw;
            }

            std::string &L = a[0];
            int n_contraction = boost::lexical_cast<int>(a[1]);
            //REAL scaling_parameter = boost::lexical_cast<REAL>(a[2]);

            transform (L.begin (), L.end (), L.begin (), toupper);
            size_t n_shell_types = 0;
            for(std::string::iterator it = L.begin(); it != L.end(); it++) {
                char l = *it;
                if (l == 'S' || l == 'P' || l == 'D' || l == 'F') {
                    struct Shell new_orbital(l);
                    new_atom.orbitals.push_back(new_orbital);
                    n_shell_types++;
                } else {
                    throw;
                }
            }
            // Parse Exponent and Coefficients
            for(int i = 0; i < n_contraction; i++) {
                std::getline(ifs, line_buf);
                std::vector<std::string> b = split_by_space(line_buf);
                REAL exp = boost::lexical_cast<REAL>(b[0]);
                if (b.size() != n_shell_types + 1) {throw;}

                for(size_t idx = 1; idx < b.size(); idx++) {
                    REAL coeff = boost::lexical_cast<REAL>(b[idx]);
                    new_atom.orbitals[parse_done + idx - 1].exponents.push_back(exp);
                    new_atom.orbitals[parse_done + idx - 1].coeffs.push_back(coeff);
                }
            }
            parse_done = new_atom.orbitals.size();
        }
    }
    return retval;
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        std::cerr << "Nothing to Read. Stopping...\n" ;
        exit(1);
    }
    struct BasisSet result(parse_basisset_file(argv[1]) );
    for(BasisSet::basis_set_container_type::iterator it = result.container.begin();
            it != result.container.end(); it++) {
        std::cout << it->to_str(); 
    }
    return 0;
}

} //namespace MOSolver
