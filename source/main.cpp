#include "specialops.h"
#include <mpi.h>
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
int param_px = 1;
int param_py = 1;
int param_pz = 1;
// acceptable residual parameter
double param_tol = 1E-8;
// maximum iterations number
int param_maxit = 100;
// thread number 
int param_nt = 1;
// basic operations testing flag
bool param_qa = false;
// path to file containing input ELLPACK matrix 
std::string param_matrix_filename = "";
// path to file containing input vector 
std::string param_vector_filename = "";

const bool DEBUG_INFO = true;
const int  MASTER_PROCESS = 0;
// total number of processes involved
int nproc, 
// process id
    rank, 
// 3 'spacial' coordinates for each process according to rank
    proc_x, proc_y, proc_z;

typedef enum {
    INVALID_PARAMETERS = 1,
    AMBIGUOUS_PROCESS_NUMBER
} errcode;

typedef struct {
    int _1;
    int _2;
} pair;


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

    if (param_px <= 0) {
        std::cerr << "input parameters error:" 
            " --px must be positive integer, but "<< 
            param_px << " was given" << std::endl;
        is_valid = false;
    }

    if (param_py <= 0) {
        std::cerr << "input parameters error:" 
            " --py must be positive integer, but "<< 
            param_py << " was given" << std::endl;
        is_valid = false;
    }

    if (param_pz <= 0) {
        std::cerr << "input parameters error:" 
            " --pz must be positive integer, but "<< 
            param_pz << " was given" << std::endl;
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

// returns pair like [i_begin, i_end) -- the second index is exclusive
pair get_index_range(const int total_elem_number, const int rank, 
    const int nproc) {

    pair result;
    result._1 = rank*(total_elem_number/nproc) + std::min(rank, 
        total_elem_number%nproc);
    result._2 = result._1 + total_elem_number/nproc + 
        ((rank<total_elem_number%nproc)? 1: 0);

    return result;
}


int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    // used for param_p(x|y|z), param_n(x|y|z) and param_nt send/receive
    int shared_params[7];
    // used for mapping `cell number` -> `process-holder`
    int *part;
    // used for mapping `local cell number` -> `global cell number`
    int *l2g;

    if (rank == MASTER_PROCESS) {

        int option_char = -1, 
            option_idx = 0;
        std::string help_message = 
            "Usage: <binary> [OPTIONS], where OPTIONS include\n"
            "[--nx=<generated matrix x size>] \n"
            "[--ny=<generated matrix y size>] \n"
            "[--nz=<generated matrix z size>] \n"
            "[--px=<processes number for x axis>] \n"
            "[--py=<processes number for y axis>] \n"
            "[--pz=<processes number for z axis>] \n"
            "[--tol=<acceptable residual> | -t <...>] \n" 
            "[--maxit=<maximum solver iterations number> | -i <...>] \n"
            "[--nt=<thread number> | -n <...>] \n"
            "[--qa] \n"
            "[--help | -h] \n"
            "[<path to ELLPACK matrix file> <path to vector file>]";

        while (true) {

            static struct option long_options[] = {
                { "nx",    required_argument, 0,  1  }, 
                { "ny",    required_argument, 0,  2  },
                { "nz",    required_argument, 0,  3  },
                { "px",    required_argument, 0,  4  },
                { "py",    required_argument, 0,  5  },
                { "pz",    required_argument, 0,  6  },
                { "qa",    no_argument,       0,  7  },
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
                    param_px = atoi(optarg);
                    break;

                case 5:
                    param_py = atoi(optarg);
                    break;

                case 6:
                    param_pz = atoi(optarg);
                    break;

                case 7:
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

        // If some options were specified after named parameters, 
        // first one will be treated as a path for matrix input, 
        // second will be used as a path for vector input, 
        // anything else will be ignored.
        // If user has not provided two filenames (i.e. has set only one), 
        // input filenames won't be processed.
        if (optind < argc and (argc-optind) >= 2) {
            param_matrix_filename = argv[optind];
            param_vector_filename = argv[optind+1];
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
                "  int param_px = %d\n"
                "  int param_py = %d\n"
                "  int param_pz = %d\n"
                "  double param_tol = %lf\n"
                "  int param_maxit = %d\n"
                "  int param_nt = %d\n"
                "  bool param_qa = %d\n",
                param_nx, 
                param_ny, 
                param_nz, 
                param_px, 
                param_py, 
                param_pz, 
                param_tol, 
                param_maxit, 
                param_nt, 
                (int) param_qa
            );

            std::cout << "  string param_matrix_filename = \"" <<
                param_matrix_filename << "\"" << std::endl;
            std::cout << "  string param_vector_filename = \"" <<
                param_vector_filename << "\"" << std::endl;
        }

        if (!validate_parameters()) {
            MPI_Abort(MPI_COMM_WORLD, INVALID_PARAMETERS);
            exit(INVALID_PARAMETERS);
        }

        if (nproc != param_px*param_py*param_pz) {

            std::cerr << "error: process number dedicated" <<
                " for program is not equal to px*py*pz" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, AMBIGUOUS_PROCESS_NUMBER);
            exit(AMBIGUOUS_PROCESS_NUMBER);
        }

        shared_params[0] = param_px;
        shared_params[1] = param_py;
        shared_params[2] = param_pz;
        shared_params[3] = param_nx;
        shared_params[4] = param_ny;
        shared_params[5] = param_nz;
        shared_params[6] = param_nt;

        ellpack_matrix em;
        std::vector<double> b;
        std::size_t N;

        // generate data if input files haven't been specified correctly
        if (param_matrix_filename.empty() or param_vector_filename.empty()) {

            N = param_nx * param_ny * param_nz;
            em = so::generate_diag_dominant_matrix(param_nx, 
                param_ny, param_nz);
            b.reserve(N);
            for (std::size_t i = 0; i < N; ++i) {
                b.push_back(std::cos(i));
            }

        } else {

            em = so::read_ellpack_matrix(param_matrix_filename);
            b = so::read_vector(param_vector_filename);
            N = b.size();
            if (DEBUG_INFO) {

                plain_matrix pm = so::ellpack2plain(em, 4);
                std::cout << "input matrix: " << std::endl;
                for (std::size_t i = 0; i < pm.rows.size(); ++i) {
                    for (std::size_t j = 0; j < pm.rows[i].size(); ++j) {
                        std::cout << pm.rows[i][j] << " ";
                    }
                    std::cout << std::endl;
                }

                std::cout << "input vector: " << std::endl;
                for (std::size_t j = 0; j < b.size(); ++j) {
                    std::cout << b[j] << " " << std::endl;
                }
            }
        }

        // time measurements storage value
        double t;
        // basic operations testing
        // executing three times each for minimizing randomness
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
                    std::cerr << "--qa: axpby assertion error: result in " << 
                        i << "th position (" << axpby_result[i] << 
                        ") is much different from expected (" << 
                        axpby_test[i] << ")" << std::endl;
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
            const double norm_result = 
                std::sqrt(so::dot(spmv_result, spmv_result));
            const double norm_test = std::sqrt(so::dot(spmv_test, spmv_test));
            if (std::abs(norm_result - norm_test) > EPS) {
                std::cerr << "--qa: spmv assertion error: result norm (" << 
                    norm_result << ") is much different from expected norm (" << 
                    norm_test << ") value" << std::endl;
            }

            // t = omp_get_wtime();
            // std::vector<double> copy_test1(N);
            // so::copy_vector(copy_test1, x);
            // t = omp_get_wtime() - t;
            // std::cout << "Parallel copy vector operation has been done in " << 
            //     t*1000 << " ms" << std::endl; 

            // t = omp_get_wtime();
            // std::vector<double> copy_test2(N);
            // copy_test2 = x;
            // t = omp_get_wtime() - t;
            // std::cout << "Ordinary copy vector operation has been done in " << 
            //     t*1000 << " ms" << std::endl; 

            // for (std::size_t i = 0; i < N; ++i) {
            //     if (copy_test1[i] != copy_test2[i]) {
            //         std::cerr << "--qa: copy assertion error: result in " << 
            //             i << "th position of first vector (" << copy_test1[i] << 
            //             ") is not equal to appropriate element of second "
            //             "vector (" << axpby_test[i] << ")" << std::endl;
            //         break;
            //     }
            // }
        }

        // t = omp_get_wtime();
        // std::vector<double> solution = 
        //     so::cg_solve(em, b, param_tol, param_maxit);
        // t = omp_get_wtime() - t;
        // std::cout << "Solver is finished in " << t*1000 << " ms" << std::endl;
    }

    // MPI_Bcast(static_cast<void *>(&param_px), 1, MPI_INT, 
    //     MASTER_PROCESS, MPI_COMM_WORLD);
    // MPI_Bcast(static_cast<void *>(&param_py), 1, MPI_INT, 
    //     MASTER_PROCESS, MPI_COMM_WORLD);
    // MPI_Bcast(static_cast<void *>(&param_pz), 1, MPI_INT, 
    //     MASTER_PROCESS, MPI_COMM_WORLD);

    MPI_Bcast(shared_params, 7, MPI_INT, MASTER_PROCESS, MPI_COMM_WORLD);

    param_px = shared_params[0];
    param_py = shared_params[1];
    param_pz = shared_params[2];
    param_nx = shared_params[3];
    param_ny = shared_params[4];
    param_nz = shared_params[5];
    param_nt = shared_params[6];

    omp_set_num_threads(param_nt);
    proc_x = rank % param_px;
    proc_y = (rank % (param_px * param_py)) / param_px;
    proc_z = rank / (param_px * param_py);

    // *_range contains two indices by each axis representing start (inc) and 
    // finish (excl) index of elements that current process possess
    pair i_range = get_index_range(param_nx, proc_x, param_px);
    pair j_range = get_index_range(param_ny, proc_y, param_py);
    pair k_range = get_index_range(param_nz, proc_z, param_pz);

    const int cells_by_x = i_range._2-i_range._1;
    const int cells_by_y = j_range._2-j_range._1;
    const int cells_by_z = k_range._2-k_range._1;

    const int internal_num = cells_by_x * cells_by_y * cells_by_z;
    int halo_num = 0;
    // flags for indicating neighbour presence
    bool down = false, left = false, back = false, 
        front = false, right = false, up = false;

    if (proc_z) {
        halo_num += cells_by_x * cells_by_y;
        down = true; 
    }

    if (proc_y) {
        halo_num += cells_by_x * cells_by_z;
        left = true; 
    }

    if (proc_x) {
        halo_num += cells_by_y * cells_by_z;
        back = true;
    }

    if (proc_x < param_px-1) {
        halo_num += cells_by_y * cells_by_z;
        front = true;
    }

    if (proc_y < param_py-1) {
        halo_num += cells_by_x * cells_by_z;
        right = true;
    }

    if (proc_z < param_pz-1) {
        halo_num += cells_by_x * cells_by_y;
        up = true;
    }

    const int extended_num = internal_num + halo_num;
    printf("extended for process #%d is %d\n", rank, extended_num);
    part = new int[extended_num];
    l2g  = new int[extended_num];

    int idx = 0;
    for (idx; idx < internal_num; ++idx) {
        part[idx] = rank;
    }

    if (down) {
        int neigh_rank = rank - param_px*param_py;
        for (int i = 0; i < cells_by_x * cells_by_y; ++i) {
           part[idx+i] = neigh_rank;
        }
        idx += cells_by_x * cells_by_y;
    }

    if (left) {
        int neigh_rank = rank - param_px;
        for (int i = 0; i < cells_by_x * cells_by_z; ++i) {
           part[idx+i] = neigh_rank;
        }
        idx += cells_by_x * cells_by_z;
    }

    if (back) {
        int neigh_rank = rank - 1;
        for (int i = 0; i < cells_by_y * cells_by_z; ++i) {
           part[idx+i] = neigh_rank;
        }
        idx += cells_by_y * cells_by_z;
    }

    if (front) {
        int neigh_rank = rank + 1;
        for (int i = 0; i < cells_by_y * cells_by_z; ++i) {
           part[idx+i] = neigh_rank;
        }
        idx += cells_by_y * cells_by_z;
    }

    if (right) {
        int neigh_rank = rank + param_px;
        for (int i = 0; i < cells_by_x * cells_by_z; ++i) {
           part[idx+i] = neigh_rank;
        }
        idx += cells_by_x * cells_by_z;
    }

    if (up) {
        int neigh_rank = rank + param_px*param_py;
        for (int i = 0; i < cells_by_x * cells_by_y; ++i) {
           part[idx+i] = neigh_rank;
        }
        idx += cells_by_x * cells_by_y;
    }

    printf("process #%d (%d,%d,%d)\n", rank, proc_x, proc_y, proc_z);
    // std::string part_representation = "";
    // for (int i = 0; i < extended_num; ++i) {
    //     part_representation += std::string(part[i]);
    // }
    // printf("process #%d part: %s\n", rank, part_representation.c_str());

    delete[] part;
    delete[] l2g;
    MPI_Finalize();
    exit(0);
}
