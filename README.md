## **Linear system conjugate gradient solver**

<p align="center">
  <img src="img/sparse-matrix.png" alt="<icon here>" width="256"/>
</p>

Imlementation of linear system solver using Jacobi preconditioner. Using OpenMP for parrallelizing.

#### Compilation 

First you must ensure you have OpenMP library installed in your system. For example, the g++ compiler option `-fopenmp` is dedicated for using the library. Project is built with the help of GNU Make. Current version allows usage of IBM and GNU compilers for C++: xlc++\_r (thread-safe instance) and g++. Code is intended for running on IBM BlueGene/P, IBM Polus and generic computers supporting at least C++98 and OpenMP 3.0+. Compilation:

For IBM BlueGene/P: `make bgp`;
for IBM Polus: `make polus`;
for computer with g++ compiler: `make gnu`.

The `Makefile` also supports C++ standard switching for GNU compilation, which is possible via changing `GNU_STANDARD` variable while starting `make` command. For instance, the C++98 standard option selecting will have the following form (note the lowercase 'c' in the name):

`make gnu GNU_STANDARD=c++98`

#### Running

Use the result of installation make command, i.e. `./solver`.

#### Cleaning out

`make clean`
