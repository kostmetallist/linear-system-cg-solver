## **Linear system conjugate gradient solver**

Imlementation of linear system solver using Jacobi preconditioner. Using OpenMP for parrallelizing.

#### Compilation 

First you must ensure you have OpenMP library installed in your system. The g++ compiler option `-fopenmp` is dedicated for using the library. Source files compiling command line string (assuming you are calling it from the project root directory):

```
g++ -fopenmp source/*
```

Later such compilation workaround will be changed with GNU make assembly.
