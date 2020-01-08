#include "specialops.h"
#include <mpi.h>
#include <omp.h>
#include <getopt.h>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
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
const int  IDXS_TAG = 1;
const int  DATA_TAG = 2;
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
std::pair<int, int> get_index_range(const int total_elem_number, 
    const int rank, const int nproc) {

    std::pair<int, int> result;
    result.first = rank*(total_elem_number/nproc) + std::min(rank, 
        total_elem_number%nproc);
    result.second = result.first + total_elem_number/nproc + 
        ((rank<total_elem_number%nproc)? 1: 0);

    return result;
}

// "gids" is for "global ids of cells"
void fill_internal_gids(int *l2g, const std::pair<int, int> i_range, 
    const std::pair<int, int> j_range, const std::pair<int, int> k_range, 
    const int param_nx, const int param_ny) {

    const int by_layer = param_nx * param_ny;
    int index = 0;
    for (int k = k_range.first; k < k_range.second; ++k) {
        for (int j = j_range.first; j < j_range.second; ++j) {
            for (int i = i_range.first; i < i_range.second; ++i) {
                l2g[index++] = k*by_layer + j*param_nx + i;
            }
        }
    }
}

// TODO remove this in production
std::string intarray2string(const int *arr, const int size) {
    std::stringstream out;
    for (int i = 0; i < size; ++i) {
        out << arr[i] << " ";
    }

    return out.str();
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

    // input values handling
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
                "  double param_tol = %.8lf\n"
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

        // time measurements storage value
        // double t;

        // basic operations testing
        // executing three times each for minimizing randomness
        // if (param_qa) {
        //     std::vector<double> x(N), y(N);
        //     for (std::size_t i = 0; i < N; ++i) {
        //         x[i] = std::cos(i*i);
        //         y[i] = std::sin(i*i);
        //     }

        //     t = omp_get_wtime();
        //     double dot_result = so::dot(x, y);
        //            dot_result = so::dot(x, y);
        //            dot_result = so::dot(x, y);
        //     t = omp_get_wtime() - t;
        //     std::cout << "Dot operation has been done in " << (t/3.0)*1000 << 
        //         " ms" << std::endl;
        //     double dot_test = 0;

        //     for (std::size_t i = 0; i < N; ++i) {
        //         dot_test += x[i] * y[i];
        //     }

        //     if (std::abs(dot_result - dot_test) > EPS) {
        //         std::cerr << "--qa: dot assertion error: result (" << 
        //             dot_result << ") is much different from expected (" << 
        //             dot_test << ") value" << std::endl;
        //     }

        //     t = omp_get_wtime();
        //     std::vector<double> axpby_result = so::axpby(x, 1.0, y, -1.0);
        //                         axpby_result = so::axpby(x, 1.0, y, -1.0);
        //                         axpby_result = so::axpby(x, 1.0, y, -1.0);
        //     t = omp_get_wtime() - t;
        //     std::cout << "Axpby operation has been done in " << (t/3.0)*1000 << 
        //         " ms" << std::endl;
        //     std::vector<double> axpby_test(N);
        //     for (std::size_t i = 0; i < N; ++i) {
        //         axpby_test[i] = x[i] - y[i];
        //     }

        //     for (std::size_t i = 0; i < N; ++i) {
        //         if (std::abs(axpby_result[i] - axpby_test[i]) > EPS) {
        //             std::cerr << "--qa: axpby assertion error: result in " << 
        //                 i << "th position (" << axpby_result[i] << 
        //                 ") is much different from expected (" << 
        //                 axpby_test[i] << ")" << std::endl;
        //             break;
        //         }
        //     }

        //     t = omp_get_wtime();
        //     std::vector<double> spmv_result = so::spmv(em, x);
        //                         spmv_result = so::spmv(em, x);
        //                         spmv_result = so::spmv(em, x);
        //     t = omp_get_wtime() - t;
        //     std::cout << "Spmv operation has been done in " << (t/3.0)*1000 << 
        //         " ms" << std::endl;
        //     const std::vector<double> spmv_test = so::spmv_consecutive(em, x);
        //     const double norm_result = 
        //         std::sqrt(so::dot(spmv_result, spmv_result));
        //     const double norm_test = std::sqrt(so::dot(spmv_test, spmv_test));
        //     if (std::abs(norm_result - norm_test) > EPS) {
        //         std::cerr << "--qa: spmv assertion error: result norm (" << 
        //             norm_result << ") is much different from expected norm (" << 
        //             norm_test << ") value" << std::endl;
        //     }
        // }

        // t = omp_get_wtime();
        // std::vector<double> solution = 
        //     so::cg_solve(em, b, param_tol, param_maxit);
        // t = omp_get_wtime() - t;
        // std::cout << "Solver is finished in " << t*1000 << " ms" << std::endl;
    }

    // receiving entered values from MASTER
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
    std::pair<int, int> i_range = get_index_range(param_nx, proc_x, param_px);
    std::pair<int, int> j_range = get_index_range(param_ny, proc_y, param_py);
    std::pair<int, int> k_range = get_index_range(param_nz, proc_z, param_pz);

    const int cells_by_x = i_range.second-i_range.first;
    const int cells_by_y = j_range.second-j_range.first;
    const int cells_by_z = k_range.second-k_range.first;

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
    // printf("extended for process #%d is %d\n", rank, extended_num);
    part = new int[extended_num];
    l2g  = new int[extended_num];

    // `part` and `l2g` filling activity
    fill_internal_gids(l2g, i_range, j_range, k_range, param_nx, param_ny);
    int idx = 0;
    for (idx; idx < internal_num; ++idx) {
        part[idx] = rank;
    }

    if (down) {
        int neigh_rank = rank - param_px*param_py;
        for (int i = 0; i < cells_by_x * cells_by_y; ++i) {
            part[idx+i] = neigh_rank;
            l2g[idx+i]  = l2g[i] - param_nx*param_ny;
        }
        idx += cells_by_x * cells_by_y;
    }

    if (left) {
        int neigh_rank = rank - param_px;
        for (int i = 0; i < cells_by_x * cells_by_z; ++i) {
            part[idx+i] = neigh_rank;
            l2g[idx+i]  = l2g[cells_by_x*cells_by_y*(i/cells_by_x) + 
                i%cells_by_x] - param_nx;
        }
        idx += cells_by_x * cells_by_z;
    }

    if (back) {
        int neigh_rank = rank - 1;
        for (int i = 0; i < cells_by_y * cells_by_z; ++i) {
            part[idx+i] = neigh_rank;
            l2g[idx+i]  = l2g[i*cells_by_x] - 1;
        }
        idx += cells_by_y * cells_by_z;
    }

    if (front) {
        int neigh_rank = rank + 1;
        for (int i = 0; i < cells_by_y * cells_by_z; ++i) {
            part[idx+i] = neigh_rank;
            l2g[idx+i]  = l2g[i*cells_by_x] + cells_by_x;
        }
        idx += cells_by_y * cells_by_z;
    }

    if (right) {
        int neigh_rank = rank + param_px;
        for (int i = 0; i < cells_by_x * cells_by_z; ++i) {
            part[idx+i] = neigh_rank;
            l2g[idx+i]  = l2g[cells_by_x*cells_by_y*(i/cells_by_x) + 
                i%cells_by_x] + cells_by_y*param_nx;
        }
        idx += cells_by_x * cells_by_z;
    }

    if (up) {
        int neigh_rank = rank + param_px*param_py;
        for (int i = 0; i < cells_by_x * cells_by_y; ++i) {
            part[idx+i] = neigh_rank;
            l2g[idx+i]  = l2g[i] + cells_by_z*(param_nx*param_ny);
        }
        idx += cells_by_x * cells_by_y;
    }

    // `matrix` will be filled only in the MASTER, for others it's just a stub
    ellpack_matrix matrix;
    // initial matrix retrieving (via generating or reading from file)
    if (rank == MASTER_PROCESS) {
        std::vector<double> b;
        std::size_t N;

        // generate data if input files haven't been specified correctly
        if (param_matrix_filename.empty() or param_vector_filename.empty()) {

            N = param_nx * param_ny * param_nz;
            matrix = so::generate_diag_dominant_matrix(param_nx, 
                param_ny, param_nz);
            b.reserve(N);
            for (std::size_t i = 0; i < N; ++i) {
                b.push_back(std::cos(i));
            }

        } else {

            matrix = so::read_ellpack_matrix(param_matrix_filename);
            b = so::read_vector(param_vector_filename);
            N = b.size();
        }
    }

    printf("process #%d part: %s\n", rank, 
        intarray2string(part, extended_num).c_str());
    printf("process #%d l2g: %s\n", rank, 
        intarray2string(l2g, extended_num).c_str());

    ellpack_matrix local_data;
    // seven elements corresponds to max number of neighbors for cell 
    // (including itself)in 3-d space 
    const char SEVEN = 7;
    local_data.idxs = std::vector< std::vector<int> >(internal_num, 
        std::vector<int>(SEVEN));
    local_data.data = std::vector< std::vector<double> >(internal_num, 
        std::vector<double>(SEVEN));

    int *idxs_storage = new int[internal_num*SEVEN];
    double *data_storage = new double[internal_num*SEVEN];
    MPI_Request matrix_req, idxs_req, data_req;
    MPI_Isend(l2g, internal_num, MPI_INT, 0, rank, MPI_COMM_WORLD, &matrix_req);
    MPI_Irecv(idxs_storage, internal_num*SEVEN, MPI_INT, MASTER_PROCESS, 
        IDXS_TAG, MPI_COMM_WORLD, &idxs_req);
    MPI_Irecv(data_storage, internal_num*SEVEN, MPI_DOUBLE, MASTER_PROCESS, 
        DATA_TAG, MPI_COMM_WORLD, &data_req);

    if (rank == MASTER_PROCESS) {
        for (int i = 0; i < nproc; ++i) {

            MPI_Status stat;
            MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);

            int process = stat.MPI_SOURCE;
            int buf_size;
            MPI_Get_count(&stat, MPI_INT, &buf_size);
            int *claimed_rows = new int[buf_size];
            MPI_Recv(claimed_rows, buf_size, MPI_INT, process, 0, 
                MPI_COMM_WORLD, &stat);
            std::vector<int> idxs_to_send = 
                so::unroll_matrix_rows(matrix.idxs, claimed_rows, buf_size);
            std::vector<double> data_to_send = 
                so::unroll_matrix_rows(matrix.data, claimed_rows, buf_size);

            MPI_Send(&idxs_to_send[0], buf_size, MPI_INT, process, 
                IDXS_TAG, MPI_COMM_WORLD);
            MPI_Send(&data_to_send[0], buf_size, MPI_DOUBLE, process, 
                DATA_TAG, MPI_COMM_WORLD);
            delete[] claimed_rows;
        }

        // cleaning out all global matrix data
        const int row_number = matrix.idxs.size();
        for (int i = 0; i < row_number; ++i) {
            matrix.idxs[i].clear();
            matrix.data[i].clear();
        }

        matrix.idxs.clear();
        matrix.data.clear();
    }

    MPI_Status stat;
    MPI_Wait(&idxs_req, &stat);
    MPI_Wait(&data_req, &stat);
    for (int i = 0; i < internal_num; ++i) {
        for (int j = 0; j < SEVEN; ++j) {
            local_data.idxs[i][j] = idxs_storage[i*SEVEN + j];
            local_data.data[i][j] = data_storage[i*SEVEN + j];
        }
    }

    delete[] idxs_storage;
    delete[] data_storage;

    delete[] part;
    delete[] l2g;
    MPI_Finalize();
    exit(0);
}
