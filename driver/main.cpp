#include "../solver/hf_solver.hpp"

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

std::vector<std::string> split_line_by_space(const std::string &str)
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

std::vector<std::string> split_line_by_eql(const std::string &str)
{
    std::vector<std::string> retval;
    boost::char_separator<char> sep("= \t\n");
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    tokenizer tokens(str,sep);
    for(tokenizer::iterator it = tokens.begin(); it != tokens.end(); it++) {
        retval.push_back(*it);
    }
    return retval;
}


void
read_basisset_file(const std::string &filename, MOSolver::BasisSet &retval)
{
    std::ifstream ifs(filename.c_str() ,std::ios::in);
    std::string line_buf;

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
        struct MOSolver::AtomBasis new_atom;
        // Parse Atom_Name
        std::getline(ifs, line_buf);
        std::vector<std::string> a = split_line_by_space(line_buf);
        if (a.size() == 0) {    break;  }
        int atomic_num = MOSolver::get_atomic_number(a[0]);
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
            a = split_line_by_space(line_buf);
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
                    struct MOSolver::Shell new_orbital(l);
                    new_atom.orbitals.push_back(new_orbital);
                    n_shell_types++;
                } else {
                    throw;
                }
            }
            // Parse Exponent and Coefficients
            for(int i = 0; i < n_contraction; i++) {
                std::getline(ifs, line_buf);
                std::vector<std::string> b = split_line_by_space(line_buf);
                MOSolver::REAL exp = boost::lexical_cast<MOSolver::REAL>(b[0]);
                if (b.size() != n_shell_types + 1) {throw;}

                for(size_t idx = 1; idx < b.size(); idx++) {
                    MOSolver::REAL coeff = boost::lexical_cast<MOSolver::REAL>(b[idx]);
                    new_atom.orbitals[parse_done + idx - 1].exponents.push_back(exp);
                    new_atom.orbitals[parse_done + idx - 1].coeffs.push_back(coeff);
                }
            }
            parse_done = new_atom.orbitals.size();
        }
    }
}

void
read_xyz_file(const std::string xyz_filename, MOSolver::System &system, std::string &title)
{
    std::ifstream ifs(xyz_filename.c_str(), std::ios::in);
    if (!ifs.is_open()) {
        std::cerr << "Error! read_xyz_file(): maybe ths file doesn't exist: " << xyz_filename << std::endl;
        throw;
    }

    size_t num_atoms = 0;
    size_t line_count = 0;
    while ( !ifs.eof() ) {
        std::string line_buf;
        std::vector<std::string> tokens;
        std::getline(ifs, line_buf); 
        line_count++;

        if (line_count == 1) {
            tokens = split_line_by_space(line_buf);
            // the first line is the number of atoms
            if (tokens.size()  == 1) { 
                num_atoms = boost::lexical_cast<int>(tokens[0]);
            } else {
                throw; 
            }
        } else if (line_count == 2) {
            title = line_buf;
        } else {
            tokens = split_line_by_space(line_buf);
            if (tokens.size() == 4) {
                int atomic_number = MOSolver::get_atomic_number(tokens[0]);
                MOSolver::REAL x = MOSolver::Angstrom2Bohr( boost::lexical_cast<MOSolver::REAL>(tokens[1]) );
                MOSolver::REAL y = MOSolver::Angstrom2Bohr( boost::lexical_cast<MOSolver::REAL>(tokens[2]) );
                MOSolver::REAL z = MOSolver::Angstrom2Bohr( boost::lexical_cast<MOSolver::REAL>(tokens[3]) );
                system.add_atom(atomic_number, x, y, z);
            } else {
                continue;
            }
        }
    }
    if (system.size() != num_atoms) {
        std::cerr << "Number of atoms are not matched" << std::endl;
        throw;
    }

    return;
}

std::string remove_comment(const std::string &line)
{
    std::string ret;
    for(std::string::const_iterator it = line.begin(); it != line.end(); it++) {
        if (*it == '#') { break; }
        ret.push_back(*it);
    }
    return ret;
}

struct execute_context {
    std::string method;
    std::string input_file;
    std::string geometry_file;
    std::string basisset_file;
    int charge;
    int multiplicity;
    std::string title_in_geometry_file;
};

void read_input(const std::string input_filename, struct execute_context &params)
{
    std::ifstream ifs(input_filename.c_str(), std::ios::in);
    if (!ifs.is_open()) {
        std::cerr << "Error! read_input(): maybe the file doesn't exist: " << input_filename << std::endl;
        throw;
    }
    size_t line_count = 0;
    while ( !ifs.eof() ) {
        std::string line_buf;
        std::getline(ifs, line_buf); 
        line_count++;

        line_buf =  remove_comment(line_buf);
        std::vector<std::string> tokens = split_line_by_eql(line_buf);
        if (tokens.size() == 2) {
            std::string name, val;
            std::transform(tokens[0].begin(), tokens[0].end(), std::back_inserter(name), ::tolower);
            std::transform(tokens[1].begin(), tokens[1].end(), std::back_inserter(val) , ::tolower);

            if (name == "method") {
                params.method = "HF";
            } else if(name == "system") {
                params.geometry_file = val;
            } else if(name == "basis") {
                params.basisset_file = val;
            } else if(name == "charge") {
                params.charge = boost::lexical_cast<int>(val);
            } else if(name == "nspin") {
                params.multiplicity = boost::lexical_cast<int>(val);
            } else {
                std::cerr << boost::format("line %3d in %s: Unknown parameter: %s.") 
                    % line_count % input_filename % tokens[0]  << std::endl;
                throw;
            }
        }
    }
    return;
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        std::cerr << "Input file is not specifies" << std::endl;
        throw;
    }
    std::string input_filename(argv[1]);
    struct execute_context params;
    read_input(input_filename, params);

    MOSolver::System system(params.charge, params.multiplicity);
    read_xyz_file(params.geometry_file, system, params.title_in_geometry_file);

    MOSolver::BasisSet basis_set(params.basisset_file);
    read_basisset_file(params.basisset_file, basis_set);
    MOSolver::CGTOs bfs = generate_bfs(system, basis_set);

    for(size_t i = 0; i < system.atom_list_.size(); i++) {
        std::cout << boost::format("%-2d\t%-3.6f\t%-3.6f\t%-3.6f") 
            % system.atom_list_[i].atomic_number 
            % system.atom_list_[i].center[0]
            % system.atom_list_[i].center[1]
            % system.atom_list_[i].center[2] << std::endl;
    }
    if (system.nspin() == 0) {
        MOSolver::rhf(bfs,system);
    } else {
        MOSolver::uhf(bfs,system);
    }
    return 0;
}
