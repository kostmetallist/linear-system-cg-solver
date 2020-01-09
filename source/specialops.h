#ifndef SPECIALOPS_H
#define SPECIALOPS_H
#include <iostream>
#include <vector>
#include <string>

#define EPS 1E-7
#define PI  3.141592

typedef struct {
    std::vector< std::vector<double> > rows;
} plain_matrix;

typedef struct {
    std::vector< std::vector<int> > idxs;
    std::vector< std::vector<double> > data;
} ellpack_matrix;


namespace so {

    std::vector<double> axpby(const std::vector<double> &x, const double a, 
        const std::vector<double> &y, const double b);
    std::vector<double> axpby_consecutive(const std::vector<double> &x, 
        const double a, const std::vector<double> &y, const double b);
    void axpby_mpi(std::vector<double> &res, const std::vector<double> &x, 
        const double a, const std::vector<double> &y, const double b);
    double dot(const std::vector<double> &x, const std::vector<double> &y);
    double dot_consecutive(const std::vector<double> &x, 
        const std::vector<double> &y);
    std::vector<double> spmv(const ellpack_matrix &matrix, 
        const std::vector<double> &vec);
    std::vector<double> spmv_consecutive(const ellpack_matrix &matrix, 
        const std::vector<double> &vec);
    ellpack_matrix plain2ellpack(const plain_matrix &matrix, 
        const std::size_t max_nonzero);
    plain_matrix ellpack2plain(const ellpack_matrix &matrix, 
        const std::size_t resulting_column_num);
    ellpack_matrix derive_diagonal(const ellpack_matrix &matrix, 
        const bool inverse_elements);
    std::vector<double> cg_solve(const ellpack_matrix &matrix, 
        const std::vector<double> &right_side, const double tolerance, 
        const int max_iterations);

    ellpack_matrix read_ellpack_matrix(const std::string path);
    std::vector<double> read_vector(const std::string path);
    ellpack_matrix generate_diag_dominant_matrix(const int nx, const int ny, 
        const int nz);

    template <class T> 
    std::vector<T> unroll_matrix_rows( const std::vector< std::vector<T> > &mat, 
        const int *claimed_rows, const int claimed_num) {

        if (mat.empty()) {
            std::cerr << "unroll_matrix: empty matrix given" << std::endl;
            return std::vector<T>();
        }

        // assuming all nested vectors have the same length
        const int W = mat[0].size();
        const int H = mat.size();
        std::vector<T> linear = std::vector<T>(W*claimed_num);
        int claimed_idx = 0;
        for (int i = 0; i < H and claimed_idx < claimed_num; ++i) {
            if (i == claimed_rows[claimed_idx]) {
                for (int j = 0; j < W; ++j) {
                    linear[claimed_idx*W+j] = mat[i][j];
                }
                claimed_idx++;
            }
        }

        return linear;
    };
}

#endif
