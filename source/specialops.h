#ifndef SPECIALOPS_H
#define SPECIALOPS_H
#include <vector>

typedef struct {
    /* data */
} sparse_matrix;

namespace so {

    std::vector<double> axpby(std::vector<double> &x, double a, 
        std::vector<double> &y, double b);
    double dot(std::vector<double> &x, std::vector<double> &y);
}

#endif