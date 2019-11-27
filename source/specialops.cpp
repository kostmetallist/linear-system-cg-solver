#include "specialops.h"
#include <mpi.h>
#include <omp.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>

#define  CUSTOM_COPY 0

namespace so {

    std::vector<double> axpby(const std::vector<double> &x, const double a, 
        const std::vector<double> &y, const double b) {

        if (x.size() != y.size()) {
            std::cerr << "so::axpby: different vector sizes" << std::endl;
            return std::vector<double>();
        }

        std::size_t size = x.size();
        std::vector<double> result(size);
        #pragma omp parallel for
        for (int i = 0; i < size; ++i) {
            result[i] = a*x[i] + b*y[i];
        }

        return result;
    }

    std::vector<double> axpby_consecutive(const std::vector<double> &x, const double a, 
        const std::vector<double> &y, const double b) {

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

    void axpby_mpi(std::vector<double> &res, const std::vector<double> &x, 
        const double a, const std::vector<double> &y, const double b) {

        if (x.size() != y.size()) {
            std::cerr << "so::axpby_mpi: different x and y vector sizes" 
                << std::endl;
            return;
        }

        if (res.size() != x.size()) {
            std::cerr << "so::axpby_mpi: result and input vector sizes"
                " do not match" << std::endl;
            return;
        }

        const std::size_t size = x.size();
        #pragma omp parallel for
        for (int i = 0; i < size; ++i) {
            res[i] = a*x[i] + b*y[i];
        }
    }

    double dot(const std::vector<double> &x, const std::vector<double> &y) {

        if (x.size() != y.size()) {
            std::cerr << "so::dot: different vector sizes" << std::endl;
            return 0;
        }

        double dot_result = 0;
        #pragma omp parallel for reduction(+:dot_result)
        for (int i = 0; i < x.size(); ++i) {
            dot_result += x[i] * y[i];
        }

        return dot_result;
    }

    double dot_consecutive(const std::vector<double> &x, const std::vector<double> &y) {

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

    // defines A*x operation, where A represents matrix with N columns and
    // x is a N-row [transposed] vector
    std::vector<double> spmv(const ellpack_matrix &matrix, 
        const std::vector<double> &vec) {

        if (matrix.idxs.size() < 1) {
            std::cerr << "so::spmv: cannot process empty matrix" << std::endl;
            return std::vector<double>();
        } else if (matrix.idxs[0].size() > vec.size()) {
            std::cerr << "so::spmv: matrix column number a priori"
                " exceeds vector size" << std::endl;
            return std::vector<double>();
        }

        std::vector<double> result(matrix.idxs.size());
        #pragma omp parallel for 
        for (int i = 0; i < matrix.idxs.size(); ++i) {

            double result_i = 0;
            int last_extracted = -1;
            for (std::size_t j = 0; j < matrix.idxs[i].size(); ++j) {

                int idx_extracted = matrix.idxs[i][j];
                if (idx_extracted == last_extracted) {
                    break;
                } else {
                    result_i += matrix.data[i][j] * vec[idx_extracted];
                    last_extracted = idx_extracted;
                }
            }

            result[i] = result_i;
        }

        return result;
    }

    std::vector<double> spmv_consecutive(const ellpack_matrix &matrix, 
        const std::vector<double> &vec) {

        if (matrix.idxs.size() < 1) {
            std::cerr << "so::spmv_consecutive: cannot process"
                " empty matrix" << std::endl;
            return std::vector<double>();
        } else if (matrix.idxs[0].size() > vec.size()) {
            std::cerr << "so::spmv_consecutive: matrix column number a priori"
                " exceeds vector size" << std::endl;
            return std::vector<double>();
        }

        std::vector<double> result(matrix.idxs.size());
        for (std::size_t i = 0; i < matrix.idxs.size(); ++i) {

            double result_i = 0;
            int last_extracted = -1;
            for (std::size_t j = 0; j < matrix.idxs[i].size(); ++j) {

                int idx_extracted = matrix.idxs[i][j];
                if (idx_extracted == last_extracted) {
                    break;
                } else {
                    result_i += matrix.data[i][j] * vec[idx_extracted];
                    last_extracted = idx_extracted;
                }
            }

            result[i] = result_i;
        }

        return result;
    }

    ellpack_matrix plain2ellpack(const plain_matrix &matrix, 
        const std::size_t max_nonzero) {

        std::size_t rows_number = matrix.rows.size();
        std::vector< std::vector<int> > idxs(rows_number);
        std::vector< std::vector<double> > data(rows_number);

        for (std::size_t i = 0; i < rows_number; ++i) {

            std::vector<int> idxs_row(max_nonzero);
            std::vector<double> data_row(max_nonzero);
            const std::vector<double> &origin_row = matrix.rows[i];
            // index for iterating over idxs and data
            std::size_t k = 0;
            for (std::size_t j = 0; 
                 k < max_nonzero && j < origin_row.size(); 
                 ++j) {

                double value = origin_row[j];
                // checking whether the value is nonzero
                if (std::abs(value) > EPS) {
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

            idxs[i] = idxs_row;
            data[i] = data_row;
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
            const std::vector<int> &ell_idxs = matrix.idxs[i];
            const std::vector<double> &ell_data = matrix.data[i];
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

    // gets the diagonal elements from the matrix and form a 
    // new instance of ellpack_matrix with only one non-zero element per row;
    // inverse_elements specifies whether the resulting matrix will be inverted
    ellpack_matrix derive_diagonal(const ellpack_matrix &matrix, 
        const bool inverse_elements) {

        ellpack_matrix empty_matrix;
        empty_matrix.idxs = std::vector< std::vector<int> >(); 
        empty_matrix.data = std::vector< std::vector<double> >();
        if (matrix.idxs.size() < 1) {
            std::cerr << "so::derive_diagonal:"
                " cannot process empty matrix" << std::endl;
            return empty_matrix;
        }

        const std::size_t n_rows = matrix.idxs.size();
        std::vector< std::vector<int> > idxs(n_rows, std::vector<int>(1));
        std::vector< std::vector<double> > data(n_rows, std::vector<double>(1));
        for (std::size_t i = 0; i < n_rows; ++i) {

            const std::vector<int> &i_idxs = matrix.idxs[i];
            for (std::size_t j = 0; j < i_idxs.size(); ++j) {

                std::size_t col_idx = i_idxs[j];
                if (col_idx >= i) {

                    if (col_idx == i)
                        data[i][0] = inverse_elements? 
                            (1.0 / matrix.data[i][j]): matrix.data[i][j];
                    if (col_idx > i)
                        data[i][0] = 0;
                    break;
                }
            }

            idxs[i][0] = i;
        }

        ellpack_matrix result;
        result.idxs = idxs; 
        result.data = data;
        return result;
    }

    std::vector<double> cg_solve(const ellpack_matrix &matrix, 
        const std::vector<double> &right_side, const double tolerance, 
        const int max_iterations) {

        if (tolerance < 0) {
            std::cerr << "so::cg_solve: tolerance must be"
                " positive real number" << std::endl;
            return std::vector<double>();
        }

        if (max_iterations < 1) {
            std::cerr << "so::cg_solve: max_iterations must be an integer"
                " greater or equal to 1" << std::endl;
            return std::vector<double>();
        }

        if (matrix.idxs.size() != right_side.size()) {
            std::cerr << "so::cg_solve: given matrix row number must match the"
                " right side size" << std::endl;
            return std::vector<double>();
        }

        const std::size_t n_rows = right_side.size();
        const ellpack_matrix &inv_diag = derive_diagonal(matrix, true);

        std::vector<double> x(n_rows, 0);
        // in general, if initial x is not 0-vector, r will be equal 
        // `axpby(right_side, 1, spmv(matrix, x), -1)` but for 
        // optimisation reasons it is omitted in this case 
        std::vector<double> r = right_side;
        // std::vector<double> r = axpby(right_side, 1, spmv(matrix, x), -1);
        std::vector<double> p(n_rows, 0);

        double ro_prev = 0, ro_curr = 0;
        bool convergence = false;
        int k = 1;
        do {

            const std::vector<double> &z = spmv(inv_diag, r);
            ro_curr = dot(r, z);

            if (k == 1) {
                p = z;
            } else {
                double beta = ro_curr / ro_prev;
                p = axpby(z, 1, p, beta);
            }

            const std::vector<double> &q = spmv(matrix, p);
            const double alpha = ro_curr / dot(p, q);
            x = axpby(x, 1, p, alpha);
            r = axpby(r, 1, q, -alpha);

            if (ro_curr < tolerance or k >= max_iterations) {

                std::cout << "ro_curr < tolerance: " << 
                    ((ro_curr < tolerance)? "true": "false") << std::endl;
                std::cout << "k >= max_iterations: " << 
                    ((k >= max_iterations)? "true": "false") << ", k = " << 
                    k << std::endl;
                convergence = true;
            } else {
                k++;
                ro_prev = ro_curr;
            }

        } while (not convergence);

        return x;
    }

    ellpack_matrix read_ellpack_matrix(const std::string path) {

        ellpack_matrix empty_matrix = { 
            std::vector< std::vector<int> >(), 
            std::vector< std::vector<double> >()
        };

        std::ifstream in(path.c_str());
        if (!in) {
            std::cerr << "so::read_ellpack_matrix: an issue"
                " occurred while creating ifstream" << std::endl;
            return empty_matrix;
        }

        std::string input;
        if (!std::getline(in, input)) {
            std::cerr << "so::read_ellpack_matrix: first string should be"
                " valid integer representing row number" << std::endl;
            return empty_matrix;
        }

        int rows_input;
        std::istringstream(input) >> rows_input;
        if (!std::getline(in, input)) {
            std::cerr << "so::read_ellpack_matrix: second string should be"
                " valid integer representing column number" << std::endl;
            return empty_matrix;
        }

        int cols_input;
        std::istringstream(input) >> cols_input;
        if (rows_input <= 0 || cols_input <= 0) {
            std::cerr << "so::read_ellpack_matrix: row and column numbers"
                " must be positive integers" << std::endl;
            return empty_matrix;
        }

        const std::size_t rows_number = static_cast<std::size_t>(rows_input);
        const std::size_t cols_number = static_cast<std::size_t>(cols_input);
        std::vector< std::vector<int> > idxs(rows_number, 
            std::vector<int>(cols_number));
        std::vector< std::vector<double> > data(rows_number, 
            std::vector<double>(cols_number));
        int line_counter = 0;
        const int lines_by_chunk = rows_number * cols_number;
        bool is_data_chunk = false;
        while (std::getline(in, input)) {

            if (is_data_chunk) {

                double elem;
                std::istringstream(input) >> elem;
                data[line_counter/cols_number][line_counter%cols_number] = elem;

            } else {

                int idx;
                std::istringstream(input) >> idx;
                if (idx < 0) {
                    std::cerr << "so::read_ellpack_matrix: found negative"
                        " ellpack index which is illegal" << std::endl;
                    return empty_matrix;
                }

                idxs[line_counter/cols_number][line_counter%cols_number] = idx;
                if (line_counter == lines_by_chunk - 1) {
                    is_data_chunk = true;
                    // line_counter will be increased by 1 after exiting to 
                    // external conditional branch. Hence, when execution reach
                    // `if (is_data_chunk)` block, line_counter will be zeroed.
                    line_counter = -1;
                }
            }

            line_counter++;
        }

        // checking that both ellpack parts were processed AND 
        // data section was read completely
        if (!is_data_chunk && line_counter < lines_by_chunk) {
            std::cerr << "so::read_ellpack_matrix: input file contains" 
                " insufficient amount of data (expected " << 
                2*lines_by_chunk+2 << " lines in file)" << std::endl;
            return empty_matrix;
        }

        ellpack_matrix result;
        result.idxs = idxs;
        result.data = data;
        return result;
    }

    // opens a file by path specified as an argument and parses every 
    // string to appropriate `double` element of resulting vector
    std::vector<double> read_vector(const std::string path) {

        std::vector<double> empty_vector;
        std::ifstream in(path.c_str());
        if (!in) {
            std::cerr << "so::read_vector: an issue"
                " occurred while creating ifstream" << std::endl;
            return empty_vector;
        }

        std::string input;
        if (!std::getline(in, input)) {
            std::cerr << "so::read_vector: first string must be valid positive"
                " integer representing vector elements number" << std::endl;
            return empty_vector;
        }

        int n_rows, line_counter = 0;
        std::istringstream(input) >> n_rows;
        std::vector<double> result(n_rows);
        while (std::getline(in, input)) {

            double elem;
            std::istringstream(input) >> elem;
            result[line_counter++] = elem;
        }

        if (line_counter != n_rows) {
            std::cerr << "so::read_vector: elements number specified in header"
                " and elements number been read are not matching" << std::endl;
            return empty_vector;
        }

        return result;

    }

    ellpack_matrix generate_diag_dominant_matrix(const int nx, 
        const int ny, const int nz) {

        if (nx <= 0 || ny <= 0 || nz <= 0) {
 
            std::cerr << "so::generate_diag_dominant_matrix: input parameters" 
                " must be positive integers" << std::endl;
            ellpack_matrix empty_matrix;
            empty_matrix.idxs = std::vector< std::vector<int> >();
            empty_matrix.data = std::vector< std::vector<double> >();
            return empty_matrix;
        }

        // Three-dimensional regular mesh topology implies that 
        // every node can have at most 6 neighbours. Links with another
        // nodes are mapped onto matrix elements. Thus, every row may 
        // have max 7 non-zero elements, including diagonal one.
        const int max_nonzero = 7;
        const int total_mesh_nodes = nx * ny * nz;
        std::vector< std::vector<int> > idxs(total_mesh_nodes, 
            std::vector<int>(max_nonzero));
        std::vector< std::vector<double> > data(total_mesh_nodes, 
            std::vector<double>(max_nonzero));
        for (int i = 0; i < total_mesh_nodes; ++i) {

            int filled_indices = 0, 
                last_inserted = 0, 
                diagonal_idx = 0;
            if (i/(nx*ny) > 0) {
                last_inserted = i-nx*ny;
                idxs[i][filled_indices++] = last_inserted;
            }

            if (i%(nx*ny) > nx-1) {
                last_inserted = i-nx;
                idxs[i][filled_indices++] = last_inserted;
            }

            if (i%nx > 0) {
                last_inserted = i-1;
                idxs[i][filled_indices++] = last_inserted;
            }

            last_inserted = i;
            diagonal_idx = filled_indices;
            idxs[i][filled_indices++] = last_inserted;

            if (i%nx < nx-1) {
                last_inserted = i+1;
                idxs[i][filled_indices++] = last_inserted;
            }

            if (i%(nx*ny) < nx*(ny-1)) {
                last_inserted = i+nx;
                idxs[i][filled_indices++] = last_inserted;
            }

            if (i/(nx*ny) < nz-1) {
                last_inserted = i+nx*ny;
                idxs[i][filled_indices++] = last_inserted;
            }

            // in case neighbours number is less than max allowed, 
            // make a padding for indices vector
            if (filled_indices < max_nonzero) {
                for (int j = 0; j < max_nonzero-filled_indices; ++j) {
                    idxs[i][filled_indices+j] = last_inserted;
                }
            }

            double non_diag_abs_sum = 0;
            for (int j = 0; j < filled_indices; ++j) {

                if (j == diagonal_idx) { continue; }
                double generated_value = cos(PI + i*idxs[i][j]);
                data[i][j] = generated_value;
                non_diag_abs_sum += std::abs(generated_value);
            }

            // making generated matrix diagonal dominant
            data[i][diagonal_idx] = 1.5 * non_diag_abs_sum;
        }

        ellpack_matrix result;
        result.idxs = idxs;
        result.data = data;
        return result;
    }
}