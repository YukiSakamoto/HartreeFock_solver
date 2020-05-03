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
- Eigen (temlate library for the linear algebra)

### Linux and MacOS

 ``` shell
 mkdir build	# out-of-source build is recommended
 cd build
 cmake ..
 make 
 make test
 ```

### Options

* Loop optimization (default is OFF)
```
cmake .. -DLOOP_OPT=ON
```

### Other Environment

To Be Confirmed


Running
-------
 ```shell
 cd build/examples
 ../MOSolver methane.inp
 ```
 
Input format
-----------
To run this program, we have to prepare 3 files: method, geometry, basis-set.

### Method file

Specify the calculation method, parameters, and the filenames of gemetry and basis-set.

```
# Methane
Method =  hf			# Currently, only hf is available
system =  methane.xyz	# methane
basis  =  sto3g.dat		# Specify the basis_set file
nspin  =  0				# Spin multiplicity 
charge =  0		
```

### Geometry file

Molecular geometry is specified by xyz format
The unit is in angstrom.

```
5
Methane
C    0. 0. 0.
H    0.000000    0.000000    1.083010
H    0.000000    1.021071   -0.361003
H    0.884274   -0.510536   -0.361003
H   -0.884274   -0.510536   -0.361003
```

### Basis set file

Basis set is specified by the following format.
This format is the same as that used for GEN keyword in Gaussian.
The data can be obtained from [Basis Set Exchange](https://www.basissetexchange.org/).

```
****
H     0 
S   3   1.00
      3.42525091             0.15432897       
      0.62391373             0.53532814       
      0.16885540             0.44463454       
****
He     0 
S   3   1.00
      6.36242139             0.15432897       
      1.15892300             0.53532814       
      0.31364979             0.44463454       
****
```

 
References
------------

* Hartree Fock Methods: 

	Szabo, Attila, and Neil S. Ostlund. "Modern quantum chemistry: introduction to advanced electronic theory." NY: McGraw-Hill (1989).

* Molecular Integral

	Taketa, Hiroshi, Sigeru Huzinaga, and Kiyosi O-ohata. "Gaussian-Expansion Methods for Molecular Integrals." Journal of the Physical Society of Japan 21.11 (1966): 2313-2324.
