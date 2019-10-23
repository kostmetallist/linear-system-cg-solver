#include "specialops.h"
#include <omp.h>
#include <getopt.h>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>

int param_nx = 10;
int param_ny = 10;
int param_nz = 10;
// acceptable residual parameter
double param_tol = 1E-9;
// maximum iterations number
int param_maxit = 100;
// thread number 
int param_nt = 1;
// basic operations testing flag
bool param_qa = false;
// input matrix filename
std::string param_input_filename = "";
const bool DEBUG_INFO = true;


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

    if (DEBUG_INFO) {

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

        std::cout << "  string param_input_filename = \"" <<
            param_input_filename << "\"" << std::endl;
    }

    if (!validate_parameters()) {
        exit(1);
    }

    omp_set_num_threads(param_nt);
    ellpack_matrix em;
    std::vector<double> b;
    std::size_t N;

    // perform data generating if input file haven't been specified
    if (param_input_filename.empty()) {

        N = param_nx * param_ny * param_nz;
        em = so::generate_diag_dominant_matrix(param_nx, param_ny, param_nz);
        b.reserve(N);
        for (std::size_t i = 0; i < N; ++i) {
            b[i] = std::cos(i);
        }
    }

    // time measurements storage value
    double t;
    // basic operations testing
    // executing three times each for minimizing randomness
    // TODO add "and GENERATION MODE"
    if (param_qa) {
 
        std::vector<double> x(N), y(N);
        for (std::size_t i = 0; i < N; ++i) {
            x[i] = std::cos(i*i);
            y[i] = std::sin(i*i);
        }

        t = omp_get_wtime();
        double dot_result = so::dot(x, y);
               dot_result = so::dot(x, y);
               dot_result = so::dot(x, y);
        t = omp_get_wtime() - t;
        std::cout << "Dot operation has been done in " << (t/3.0)*1000 << 
            " ms" << std::endl;
        double dot_test = 0;

        for (std::size_t i = 0; i < N; ++i) {
            dot_test += x[i] * y[i];
        }

        if (std::abs(dot_result - dot_test) > EPS) {
            std::cerr << "--qa: dot assertion error: result (" << 
                dot_result << ") is much different from expected (" << 
                dot_test << ") value" << std::endl;
        }

        t = omp_get_wtime();
        std::vector<double> axpby_result = so::axpby(x, 1.0, y, -1.0);
                            axpby_result = so::axpby(x, 1.0, y, -1.0);
                            axpby_result = so::axpby(x, 1.0, y, -1.0);
        t = omp_get_wtime() - t;
        std::cout << "Axpby operation has been done in " << (t/3.0)*1000 << 
            " ms" << std::endl;
        std::vector<double> axpby_test(N);
        for (std::size_t i = 0; i < N; ++i) {
            axpby_test[i] = x[i] - y[i];
        }

        for (std::size_t i = 0; i < N; ++i) {
            if (std::abs(axpby_result[i] - axpby_test[i]) > EPS) {
                std::cerr << "--qa: axpby assertion error: result in " << i <<
                    "th position (" << axpby_result[i] << ") is much different"
                    " from expected (" << axpby_test[i] << ")" << std::endl;
                break;
            }
        }

        t = omp_get_wtime();
        std::vector<double> spmv_result = so::spmv(em, x);
                            spmv_result = so::spmv(em, x);
                            spmv_result = so::spmv(em, x);
        t = omp_get_wtime() - t;
        std::cout << "Spmv operation has been done in " << (t/3.0)*1000 << 
            " ms" << std::endl;
        const std::vector<double> spmv_test = so::spmv_consecutive(em, x);
        const double norm_result = std::sqrt(so::dot(spmv_result, spmv_result));
        const double norm_test = std::sqrt(so::dot(spmv_test, spmv_test));
        if (std::abs(norm_result - norm_test) > EPS) {
            std::cerr << "--qa: spmv assertion error: result norm (" << 
                norm_result << ") is much different from expected norm (" << 
                norm_test << ") value" << std::endl;
        }

        t = omp_get_wtime();
        std::vector<double> copy_test1(N);
        so::copy_vector(copy_test1, x);
        t = omp_get_wtime() - t;
        std::cout << "Parallel copy vector operation has been done in " << 
            t*1000 << " ms" << std::endl; 

        t = omp_get_wtime();
        std::vector<double> copy_test2(N);
        copy_test2 = x;
        t = omp_get_wtime() - t;
        std::cout << "Ordinary copy vector operation has been done in " << 
            t*1000 << " ms" << std::endl; 

        for (std::size_t i = 0; i < N; ++i) {
            if (copy_test1[i] != copy_test2[i]) {
                std::cerr << "--qa: copy assertion error: result in " << i <<
                    "th position of first vector (" << copy_test1[i] << 
                    ") is not equal to appropriate element of second vector (" 
                    << axpby_test[i] << ")" << std::endl;
                break;
            }
        }
    }

    t = omp_get_wtime();
    std::vector<double> solution = so::cg_solve(em, b, param_tol, param_maxit);
    t = omp_get_wtime() - t;
    std::cout << "Solver is finished in " << t*1000 << " ms" << std::endl;
    exit(0);
}