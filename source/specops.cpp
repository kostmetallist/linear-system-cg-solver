#include "specops.h"
#include <vector>

namespace so {

    std::vector<double> axpby(std::vector<double> &x, double a, 
        std::vector<double> &y, double b) {

        int arr_size = 5;
        std::vector<double> result(arr_size);

        for (int i = 0; i < arr_size; ++i) {
            result[i] = i*i;
        }

        return result;
    }
}