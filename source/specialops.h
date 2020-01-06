#ifndef SPECIALOPS_H
#define SPECIALOPS_H
#include <vector>
#include <string>

#define EPS 1E-7
#define PI  3.141592

typedef struct {
    std::vector< std::vector<double> > rows;
} plain_matrix;

typedef struct 
{
    std::vector< std::vector<int> > idxs;
    std::vector< std::vector<double> > data;
} ellpack_matrix;


namespace so {

    void copy_vector(std::vector<double> &to, const std::vector<double> &from);

    std::vector<double> axpby(const std::vector<double> &x, const double a, 
        const std::vector<double> &y, const double b);
    std::vector<double> axpby_consecutive(const std::vector<double> &x, 
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

    template <class T> T *unrollVector(std::vector<T> vec);
}

#endif
