## **Linear system conjugate gradient solver**

<p align="center">
  <img src="img/sparse-matrix.png" alt="<icon here>" width="384"/>
</p>

Imlementation of linear system solver using Jacobi preconditioner. Using OpenMP for thread parallelizing and MPI for distributed version of solver.

#### Compilation 

First you must ensure you have OpenMP library installed in your system. Also, make sure you are possible to compile and run MPI programs. The `mpicxx` compiler option `-fopenmp` is dedicated for using the OpenMP library. Project is built with the help of GNU Make. Current version allows usage of IBM and GNU compilers for MPI C++: `mpixlcxx_r` (thread-safe instance), `mpixlC` and `mpicxx`. Code is intended for running on IBM BlueGene/P, IBM Polus and generic computers supporting at least C++98, OpenMP 3.0+, MPI 2.1+. Compilation:

For IBM BlueGene/P: `make bgp`;
for IBM Polus: `make polus`;
for computer with g++ compiler: `make gnu`.

The `Makefile` also supports C++ standard switching for GNU compilation, which is possible via changing `GNU_STANDARD` variable while starting `make` command. For instance, the C++98 standard option selecting will have the following form (note the lowercase 'c' in the name):

`make gnu GNU_STANDARD=c++98`

#### Running

In order to run MPI programs, one needs to clarify the launch manager preferred for target system. It may be, for instance, `mpirun` in most cases but on IBM Polus and IBM BlueGene/P specialized commands should be used instead. For the case of `mpirun`: 

`mpirun -np <process number> ./solver`

There is also a possibility that the environment demands kernel module loading for MPI utilizing. Such work can be done via `module load <module name>` command. Look for target system MPI requirements before proceeding.

#### Cleaning out

`make clean`
