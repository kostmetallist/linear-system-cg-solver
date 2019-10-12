#ifndef SPECIALOPS_H
#define SPECIALOPS_H
#include <vector>

#define EPS 1E-6

typedef struct {
    std::vector<std::vector<double>> rows;
} plain_matrix;

typedef struct 
{
    std::vector<std::vector<int>> idxs;
    std::vector<std::vector<double>> data;
} ellpack_matrix;

namespace so {

    std::vector<double> axpby(std::vector<double> &x, double a, 
        std::vector<double> &y, double b);
    double dot(std::vector<double> &x, std::vector<double> &y);
    ellpack_matrix plain2ellpack(plain_matrix &matrix, std::size_t max_nonzero);
}

#endif