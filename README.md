Simple Implementation of Hartree Fock method
=======================

A simple implementation of Hartree Fock Method (just for my learning!).

The implementation is based on the "Modern Quantum Chemistry" by Szabo and Ostland"

Since the molecular integrals are evaluated by the methods of Taketa, Huzinaga and Oohata, this program could calculate not only the S orbitals, but also P and D orbitals (but TOO SLOW!). 

Build and Run
------------

### Requirements

- cmake (build system)
- boost (boost unit_test_framework is used for the test)
- gsl   (for special function)
- Eigen (temlate library for the linear algebra)

### Linux

Before running the CMake, expand the Eigen template library into the same directory.


 ``` shell
 cmake .
 make 
 make test
 ```

### Other Environment

To Be Confirmed


Running
-------
 ```shell
 cd example
 ./run_rhf
 ```
 
 
References
------------

* Hartree Fock Methods: 

	Szabo, Attila, and Neil S. Ostlund. "Modern quantum chemistry: introduction to advanced electronic theory." NY: McGraw-Hill (1989).

* Molecular Integral

	Taketa, Hiroshi, Sigeru Huzinaga, and Kiyosi O-ohata. "Gaussian-Expansion Methods for Molecular Integrals." Journal of the Physical Society of Japan 21.11 (1966): 2313-2324.
