#include "specialops.h"
#include <omp.h>
#include <getopt.h>
#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <stdio.h>

int param_nx = 100;
int param_ny = 100;
int param_nz = 100;
// residual parameter
double param_tol = 9.8E-9;
// maximum iterations number
int param_maxit = 1000;
// thread number 
int param_nt = 1;
// testing mode flag
bool param_qa = false;
// input matrix filename
std::string param_input_filename = "";


int main(int argc, char *argv[]) {

    int option_char = -1, 
        option_idx = 0;

    while (true) {

        static struct option long_options[] = {
            { "nx",    required_argument, 0,  1  }, 
            { "ny",    required_argument, 0,  2  },
            { "nz",    required_argument, 0,  3  },
            { "qa",    no_argument,       0,  4  },
            { "tol",   required_argument, 0, 't' },
            { "maxit", required_argument, 0, 'i' },
            { "nt",    required_argument, 0, 'n' }, 
            { "help",  no_argument,       0, 'h' },
            { 0,       0,                 0,  0  }
        };

        option_char = getopt_long(argc, argv, "t:i:n:h", 
            long_options, &option_idx);
        if (option_char == -1) {
            break;
        }

        switch (option_char) {
            case 1:
                param_nx = atoi(optarg);
                break;

            case 2:
                param_ny = atoi(optarg);
                break;

            case 3:
                param_nz = atoi(optarg);
                break;

            case 4:
                param_qa = true;
                break;

            case 't':
                param_tol = atof(optarg);
                break;

            case 'i':
                param_maxit = atoi(optarg);
                break;

            case 'n':
                param_nt = atoi(optarg);
                break;

            case 'h':
                printf(
                    "Usage: <binary> [--nx=<generated matrix x size>] "
                    "[--ny=<generated matrix y size>] "
                    "[--nz=<generated matrix z size>] "
                    "[--tol=<acceptable residual> | -t <...>] " 
                    "[--maxit=<maximum solver iterations number> | -i <...>] "
                    "[--nt=<thread number> | -n <...>] "
                    "[--qa] [--help | -h] [matrix ELLPACK input file]"
                    "\n"
                );
                break;

            case '?':
                break;
        }
    }

    // if anything was specified after named parameters, 
    // first one will be treated as an input filename, while 
    // others will be ignored
    if (optind < argc) {
        param_input_filename = argv[optind];
    }

    printf(
        "Updated params: "
        "  int param_nx = %d\n"
        "  int param_ny = %d\n"
        "  int param_nz = %d\n"
        "  double param_tol = %lf\n"
        "  int param_maxit = %d\n"
        "  int param_nt = %d\n"
        "  bool param_qa = %d\n",
        param_nx, 
        param_ny, 
        param_nz, 
        param_tol, 
        param_maxit, 
        param_nt, 
        (int) param_qa
    );

    std::cout << "  string param_input_filename = " <<
        param_input_filename << std::endl;

    std::vector<double> v1(3, 2.0), v2(3, 3.0);
    std::vector<double> vec = so::axpby(v1, 2.5, v2, 2);

    return 0;
}