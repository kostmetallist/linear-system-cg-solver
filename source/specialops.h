#ifndef SPECIALOPS_H
#define SPECIALOPS_H
#include <vector>
#include <string>

#define EPS 1E-6

typedef struct {
    std::vector< std::vector<double> > rows;
} plain_matrix;

typedef struct 
{
    std::vector< std::vector<int> > idxs;
    std::vector< std::vector<double> > data;
} ellpack_matrix;

namespace so {

    std::vector<double> axpby(std::vector<double> &x, const double a, 
        std::vector<double> &y, const double b);
    double dot(std::vector<double> &x, std::vector<double> &y);
    std::vector<double> spmv(const ellpack_matrix &matrix, 
        const std::vector<double> &vec);
    ellpack_matrix plain2ellpack(const plain_matrix &matrix, 
        const std::size_t max_nonzero);
    plain_matrix ellpack2plain(const ellpack_matrix &matrix, 
        const std::size_t resulting_column_num);

    ellpack_matrix read_ellpack_matrix(const std::string path);
    ellpack_matrix generate_diag_dominant_matrix(const int nx, const int ny, 
        const int nz);
}

#endif