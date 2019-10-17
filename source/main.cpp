#include "specialops.h"
#include <omp.h>
#include <getopt.h>
#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <stdio.h>

int param_nx = 10;
int param_ny = 10;
int param_nz = 10;
// acceptable residual parameter
double param_tol = 9.8E-9;
// maximum iterations number
int param_maxit = 1000;
// thread number 
int param_nt = 1;
// testing mode flag
bool param_qa = false;
// input matrix filename
std::string param_input_filename = "";
const bool DEBUG_MODE = true;


bool validate_parameters() {

    bool is_valid = true;
    if (param_nx <= 0) {
        std::cerr << "input parameters error:" 
            " --nx must be positive integer, but "<< 
            param_nx << " was given" << std::endl;
        is_valid = false;
    }

    if (param_ny <= 0) {
        std::cerr << "input parameters error:" 
            " --ny must be positive integer, but "<< 
            param_ny << " was given" << std::endl;
        is_valid = false;
    }

    if (param_nz <= 0) {
        std::cerr << "input parameters error:" 
            " --nz must be positive integer, but "<< 
            param_nz << " was given" << std::endl;
        is_valid = false;
    }

    if (param_tol < 0) {
        std::cerr << "input parameters error:" 
            " --tol (-t) must be positive real number, but "<< 
            param_tol << " was given" << std::endl;
        is_valid = false;
    }

    if (param_maxit <= 0) {
        std::cerr << "input parameters error:" 
            " --maxit (-i) must be positive integer, but "<< 
            param_maxit << " was given" << std::endl;
        is_valid = false;
    }

    if (param_nt <= 0) {
        std::cerr << "input parameters error:" 
            " --nt (-n) must be positive integer, but "<< 
            param_nt << " was given" << std::endl;
        is_valid = false;
    }

    return is_valid;
}

int main(int argc, char *argv[]) {

    int option_char = -1, 
        option_idx = 0;
    std::string help_message = 
        "Usage: <binary> [OPTIONS], where OPTIONS include\n"
        "[--nx=<generated matrix x size>] \n"
        "[--ny=<generated matrix y size>] \n"
        "[--nz=<generated matrix z size>] \n"
        "[--tol=<acceptable residual> | -t <...>] \n" 
        "[--maxit=<maximum solver iterations number> | -i <...>] \n"
        "[--nt=<thread number> | -n <...>] \n"
        "[--qa] \n"
        "[--help | -h] \n"
        "[<path to ELLPACK input file>]";

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
                std::cout << help_message << std::endl;
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

    if (argc == 1) {
        std::cout << help_message << std::endl;
    }

    if (DEBUG_MODE) {

        printf(
            "Updated params: \n"
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
    }

    if (!validate_parameters()) {
        exit(1);
    }

    // ellpack_matrix em = so::read_ellpack_matrix(param_input_filename);
    // plain_matrix pm = so::ellpack2plain(em, 4);
    // std::cout << "input matrix: " << std::endl;
    // for (int i = 0; i < pm.rows.size(); ++i) {
    //     for (int j = 0; j < pm.rows[i].size(); ++j) {
    //         std::cout << pm.rows[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // std::vector<double> x(4);
    // x[0] = 1.0;
    // x[1] = 2.0;
    // x[2] = 1.0;
    // x[3] = 1.0;
    // std::cout << "x vector: " << std::endl;
    // for (int i = 0; i < x.size(); ++i) {
    //     std::cout << x[i] << " ";
    // }
    // std::cout << std::endl;

    // const std::vector<double> &result = so::spmv(em, x);
    // std::cout << "result vector: " << std::endl;
    // for (int i = 0; i < result.size(); ++i) {
    //     std::cout << result[i] << " ";
    // }
    // std::cout << std::endl;

    // ellpack_matrix em = so::generate_diag_dominant_matrix(param_nx, param_ny, param_nz);
    // plain_matrix pm = so::ellpack2plain(em, 7);
    // std::cout << "generated matrix: " << std::endl;
    // for (int i = 0; i < pm.rows.size(); ++i) {
    //     for (int j = 0; j < pm.rows[i].size(); ++j) {
    //         printf("% 05.3lf ", pm.rows[i][j]);
    //     }
    //     std::cout << std::endl;
    // }

    exit(0);
}