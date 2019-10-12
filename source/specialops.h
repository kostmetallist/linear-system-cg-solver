#ifndef SPECIALOPS_H
#define SPECIALOPS_H
#include <vector>

namespace so {

    std::vector<double> axpby(std::vector<double> &x, double a, 
        std::vector<double> &y, double b);
}

#endif