#ifndef _GOAT_MATRIX_HPP
#define _GOAT_MATRIX_HPP

#include "ztk.hpp"

#include <vector>
#include <iostream>
#include <array>

namespace goat {

using ztk::GR;
using std::vector;

template<size_t K>
std::array<GR<K, 4>, 4> gr_deg4_vec_mul(const GR<K, 4>& x, const std::array<GR<K, 4>, 4>& y) {

    std::array<GR<K, 4>, 4> r;

    auto x0 = x[0]; auto x1 = x[1]; auto x2 = x[2]; auto x3 = x[3];
    auto y0 = y[0]; auto y1 = y[1]; auto y2 = y[2]; auto y3 = y[3];

    r[0] = x0*y0 - x3*y1 - x2*y2 - x1*y3;
    r[1] = x1*y0 + x0*y1 - x3*(y1 + y2) - x2*y2 - x2*y3 - x1*y3;
    r[2] = x2*y0 - x2*y3 + x1*y1 + y2*(x0 - x3) - x3*y3;
    r[3] = x3*y0 + x2*y1 + x1*y2 + x0*y3 - x3*y3;

    return r;
}

template<size_t K, size_t D>
void init_vand_matrix(const vector<GR<K,D>> &excep_seq, const size_t r, const size_t c, vector<vector<GR<K,D>>>& vand_matrix){
    // r : number of rows
    // c : number of columns
    // vandermonde matrix is of the form
    //    1   x1  x1^2  ...   x1^{c-1}
    //    1   x2  x2^2  ...   x2^{c-1}
    //    ... ... ...   ...     ...
    //    1   xr  xr^2  ...   xr^{c-1}

    assert(r <= excep_seq.size());
    // assert(c <= excep_seq.size()/2);

#ifdef TESTING
    std::cout << "Initializing vandermonde matrix\n";
#endif

    vand_matrix.resize(r);

    for(size_t i = 0; i < r; i++){
        vand_matrix[i].resize(c);
        vand_matrix[i][0] = GR<K,D>::one(); // We'll initialize this as 1 for now

        for(size_t j = 1; j < c; j++){
            vand_matrix[i][j] = excep_seq[i] * vand_matrix[i][j-1];
        }
    }
}

template<size_t K, size_t D>
void init_inv_vand_matrix(const vector<vector<GR<K,D>>> &mat, const size_t d, vector<vector<GR<K,D>>>& inv_vand_matrix){
    // Only ever going to deal with a (n-t)-(n-t) matrix
    // This is specifically for public reconstruction, so going to hardcode the dimensions
    // This currently does not handle if the matrix is not invertible... It should be

    // Need to make sure the matrix we're inverting is the right size
    assert(mat.size() >= d);
    assert(mat[0].size() >= d);

#ifdef TESTING
    std::cout << "Initializing vandermonde matrix inverse\n";
#endif

    // augmented matrix for storing computation
    vector<vector<GR<K,D>>> aug_matrix(d);
    inv_vand_matrix.resize(d);

    for(size_t i = 0; i < d; i++){
        aug_matrix[i].resize(d);
        inv_vand_matrix[i].resize(d); // answer will be stored here
        for(size_t j = 0; j < d; j++){
            aug_matrix[i][j] = mat[i][j];
            if(i==j){
                inv_vand_matrix[i][j] = GR<K,D>::one();
            } else {
                inv_vand_matrix[i][j] = GR<K,D>::zero();
            }
        }
    }

    // Do gauss elim, and put the answer directly into inv_vand_matrix
    // matrix we're going to use is [aug_matrix | inv_vand_matrix]

    // step 1: get in row echelon form
    for(size_t p = 0; p < d-1; p++){
        // p for pivot (row)
        GR<K,D> a_pp = aug_matrix[p][p];

        // If the pivot point is zero, then we need to switch it with a nonzero row
        if(a_pp == GR<K,D>::zero()){
            for(size_t new_p = p+1; new_p < d; new_p++){
                // replace the current row with a row that doesn't have a leading zero for the column we're looking at
                if(aug_matrix[new_p][p] != GR<K,D>::zero()){
                    std::cout << "Replacing a row\n";
                    vector<GR<K,D>> row_holder = aug_matrix[new_p];
                    aug_matrix[new_p] = aug_matrix[p];
                    aug_matrix[p] = row_holder;

                    row_holder = inv_vand_matrix[new_p];
                    inv_vand_matrix[new_p] = inv_vand_matrix[p];
                    inv_vand_matrix[p] = row_holder;

                    a_pp = aug_matrix[new_p][p];

                    // go back to the main loop once we've switched
                    break;
                }
            }

            // If we can't find a nonzero row, then this is an issue
            if(a_pp == GR<K,D>::zero()){
                throw std::runtime_error("issue trying to do gaussian elimination\n");
            }
        }

        // try to get it in row echelon form
        // i : row index (not pivot row)
        for(size_t i = p+1; i < d; i++){
            GR<K,D> b_ip = aug_matrix[i][p] / a_pp;
            aug_matrix[i][p] = GR<K,D>::zero();
            // j : column index
            for(size_t j = p+1; j < d; j++){
                // M[i,j] = M[i,j] - M[p,j]*M[p,p]
                aug_matrix[i][j] = aug_matrix[i][j] - (b_ip * aug_matrix[p][j]);
            }
            // inv_vand_matrix[i][p] = inv_vand_matrix[i][p] - (b_ip * inv_vand_matrix[p][p]);
            for(size_t j = 0; j < d; j++){
                inv_vand_matrix[i][j] = inv_vand_matrix[i][j] - (b_ip * inv_vand_matrix[p][j]);
            }
        }

    }

    // step 2: Back substitution
    // Going to start with the bottom row
    //      Subtract starting with the top row
    //      Then keep subtracting going down
    // Gonna start from the bottom, now we here
    for(ssize_t i = d-1; i >= 0; i--){
        GR<K,D> aii = aug_matrix[i][i];
        // for the i-th row, make the left column a 1
        if(aii != GR<K,D>::one()){
            for(size_t k = i; k < d; k++){
                aug_matrix[i][k] = aug_matrix[i][k] / aii;
            }
            for(size_t k = 0; k < d; k++){
                inv_vand_matrix[i][k] = inv_vand_matrix[i][k] / aii;
            }
        }
        // current row
        for(ssize_t j = 0; j < i; j++){
            GR<K,D> bij = aug_matrix[j][i];
            //current column
            for(size_t k = i; k < d; k++){
                aug_matrix[j][k] = aug_matrix[j][k] - (bij * aug_matrix[i][k]);
            }
            for(size_t k = 0; k < d; k++){
                inv_vand_matrix[j][k] = inv_vand_matrix[j][k] - (bij * inv_vand_matrix[i][k]);
            }
        }
    }

// #ifdef TESTING
//     std::cout << "Testing if the matrix is correctly computed\n";
//     vector<vector<GR<K,D>>> test_matrix(d);
//     bool is_correct = true;

//     // Check the matrix is correctly computed
//     for(size_t i = 0; i < d; i++){
//         test_matrix[i].resize(d);
//         for(size_t j = 0; j < d; j++){
//             test_matrix[i][j] = mat[i][0] * inv_vand_matrix[0][j];
//             for(size_t k = 1; k < d; k++){
//                 test_matrix[i][j] = test_matrix[i][j] + mat[i][k] * inv_vand_matrix[k][j];
//             }

//             is_correct = is_correct && ((i==j && test_matrix[i][j] == GR<K,D>::one()) || (i!=j && test_matrix[i][j] == GR<K,D>::zero()));
//         }
//     }

//     if(!is_correct){
//         std::cout << "There was a problem with inverting the matrix!\n";
//     }else{
//         std::cout << "vandermonde matrix was successfully inverted\n";
//     }
// #endif
}

// gaussian elimination, in case we need it
//      (to be tested)
// [M][a] = [x] solve for [a] by doing gaussian elim on [M | x]
template<size_t K, size_t D>
vector<GR<K,D>> solve_lin_eq(const vector<vector<GR<K,D>>> &mat, const vector<GR<K,D>> xs, const size_t d){

    // Need to make sure the matrix we're inverting is the right size
    assert(mat.size() >= d);
    assert(mat[0].size() >= d);
    assert(xs.size() >= d);

    // augmented matrix for storing computation
    vector<vector<GR<K,D>>> aug_matrix(d);
    vector<GR<K,D>> as(d);

    // initialize the aug_matrix as
    // [M | x]
    for(size_t i = 0; i < d; i++){
        aug_matrix[i].resize(d);
        for(size_t j = 0; j < d; j++){
            aug_matrix[i][j] = mat[i][j];
        }
        as[i] = xs[i];
    }

    // gauss_elim(aug_matrix, d, d+1);
    // step 1: get in row echelon form
    for(size_t p = 0; p < d-1; p++){
        // p for pivot (row)
        GR<K,D> a_pp = aug_matrix[p][p];

        // If the pivot point is zero, then we need to switch it with a nonzero row
        if(a_pp == GR<K,D>::zero()){
            for(size_t new_p = p+1; new_p < d; new_p++){
                // replace the current row with a row that doesn't have a leading zero for the column we're looking at
                if(aug_matrix[new_p][p] != GR<K,D>::zero()){
                    vector<GR<K,D>> row_holder = aug_matrix[new_p];
                    aug_matrix[new_p] = aug_matrix[p];
                    aug_matrix[p] = row_holder;

                    GR<K,D> elem_holder = as[new_p];
                    as[new_p] = as[p];
                    as[p] = elem_holder;

                    a_pp = aug_matrix[new_p][p];

                    // go back to the main loop once we've switched
                    break;
                }
            }

            // If we can't find a nonzero row, then this is an issue
            if(a_pp == GR<K,D>::zero()){
                throw std::runtime_error("issue trying to do gaussian elimination\n");
            }
        }

        // try to get it in row echelon form
        // i : row index (not pivot row)
        for(size_t i = p+1; i < d; i++){
            GR<K,D> b_ip = aug_matrix[i][p] / a_pp;
            aug_matrix[i][p] = GR<K,D>::zero();
            // j : column index
            for(size_t j = p+1; j < d; j++){
                // M[i,j] = M[i,j] - M[p,j]*M[p,p]
                aug_matrix[i][j] = aug_matrix[i][j] - (b_ip * aug_matrix[p][j]);
            }
            as[i] = as[i] - (b_ip * as[p]);
        }

    }

    // step 2: Back substitution
    // Going to start with the bottom row
    //      Subtract starting with the top row
    //      Then keep subtracting going down
    // Gonna start from the bottom, now we here
    for(ssize_t i = d-1; i >= 0; i--){
        GR<K,D> aii = aug_matrix[i][i];
        // for the i-th row, make the left column a 1
        if(aii != GR<K,D>::one()){
            for(size_t k = i; k < d; k++){
                aug_matrix[i][k] = aug_matrix[i][k] / aii;
            }
            as[i] = as[i] / aii;
        }
        // current row
        for(ssize_t j = 0; j < i; j++){
            GR<K,D> bij = aug_matrix[j][i];
            //current column
            for(size_t k = i; k < d; k++){
                aug_matrix[j][k] = aug_matrix[j][k] - (bij * aug_matrix[i][k]);
            }
            as[j] = as[j] - (bij * as[i]);
        }
    }

    return as;
}

// Initializes HIM
template<size_t K, size_t D>
void init_hyper_matrix(const vector<GR<K,D>> &excep_seq, const size_t r, const size_t c, const size_t n, vector<vector<GR<K,D>>>& hyper_matrix){
    // r : number of rows
    // c : number of columns

    // p^d >= 2n
    // assert((r+c) <= excep_seq.size());

#ifdef TESTING
    std::cout << "Initializing hyper invertible matrix\n";
#endif

    // We'll let alpha1, ..., alphan be the first n elem of excep_seq
    //      and beta1, ..., betan be the second n elem of excep_seq

    hyper_matrix.resize(r);

    for(size_t i = 0; i < r; i++){
        hyper_matrix[i].resize(c);
        GR<K,D> beta_i = excep_seq[i+c];

        for(size_t j = 0; j < c; j++){
            GR<K,D> mij = GR<K,D>::one(); // We'll initialize this as 1 for now
            GR<K,D> alpha_j = excep_seq[j];
            for(size_t k = 0; k < n; k++){
                if(k==j){
                    continue;
                }
                // m_{i, j} = \prod_{k \neq j} (beta_i - alpha_k)/(alpha_j - alpha_k)
                GR<K,D> alpha_k = excep_seq[k];
                mij = mij * (beta_i - alpha_k)/(alpha_j - alpha_k);
            }
            hyper_matrix[i][j] = mij;
        }
    }
}

} // goat

#endif // _GOAT_MATRIX_HPP
