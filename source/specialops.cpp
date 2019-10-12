#include "specialops.h"
#include <vector>
#include <iostream>
#include <cmath>

namespace so {

    std::vector<double> axpby(std::vector<double> &x, const double a, 
        std::vector<double> &y, const double b) {

        if (x.size() != y.size()) {
            std::cerr << "so::axpby: different vector sizes" << std::endl;
            return std::vector<double>();
        }

        std::size_t size = x.size();
        std::vector<double> result(size);
        for (std::size_t i = 0; i < size; ++i) {
            result[i] = a*x[i] + b*y[i];
        }

        return result;
    }

    double dot(std::vector<double> &x, std::vector<double> &y) {

        if (x.size() != y.size()) {
            std::cerr << "so::dot: different vector sizes" << std::endl;
            return 0;
        }

        double dot_result = 0;
        for (std::size_t i = 0; i < x.size(); ++i) {
            dot_result += x[i] * y[i];
        }

        return dot_result;
    }

    ellpack_matrix plain2ellpack(const plain_matrix &matrix, 
        const std::size_t max_nonzero) {

        std::size_t rows_number = matrix.rows.size();
        std::vector< std::vector<int> > idxs(rows_number);
        std::vector< std::vector<double> > data(rows_number);

        for (std::size_t i = 0; i < rows_number; ++i) {

            idxs[i] = std::vector<int>(max_nonzero);
            data[i] = std::vector<double>(max_nonzero);

            std::vector<double> origin_row = matrix.rows[i];
            std::vector<int> idxs_row = idxs[i];
            std::vector<double> data_row = data[i];
            // index for iterating over idxs and data
            std::size_t k = 0;
            for (std::size_t j = 0; 
                 k < max_nonzero && j < origin_row.size(); 
                 ++j) {

                double value = origin_row[j];
                // checking whether the value is nonzero
                if (abs(value) > EPS) {
                    idxs_row[k] = j;
                    data_row[k] = value;
                    k++;
                }
            }

            // if all the elements were filled, ok;
            // otherwise, create the padding containing 
            // last consistent value in the vector
            if (k != max_nonzero) {

                int tail_idxs_elem = idxs_row[k-1];
                int tail_data_elem = data_row[k-1];
                while (k < max_nonzero) {

                    // padding 
                    idxs_row[k] = tail_idxs_elem;
                    data_row[k] = tail_data_elem;
                    k++;
                }
            }
        }

        ellpack_matrix result;
        result.idxs = idxs;
        result.data = data;
        return result;
    }

    plain_matrix ellpack2plain(const ellpack_matrix &matrix, 
        const std::size_t resulting_column_num) {

        std::size_t row_number = matrix.idxs.size();
        std::size_t col_number = resulting_column_num;
        std::vector< std::vector<double> > rows(row_number);
        for (std::size_t i = 0; i < row_number; ++i) {

            std::vector<double> row(col_number);
            const std::vector<int> ell_idxs = matrix.idxs[i];
            const std::vector<double> ell_data = matrix.data[i];
            // index for ellpack structure elements referencing
            std::size_t k = 0;
            for (std::size_t j = 0; j < col_number; ++j) {

                if (k < ell_idxs.size() && j == ell_idxs[k]) {
                    row[j] = ell_data[k++];
                } else {
                    row[j] = 0;
                }
            }

            rows[i] = row;
        }

        plain_matrix result;
        result.rows = rows;
        return result;
    }
}