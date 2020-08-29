#ifndef _GOAT_ENCODING_HPP
#define _GOAT_ENCODING_HPP

#include "ztk.hpp"
#include <vector>
#include <array>

using ztk::GR;
using ztk::Z2k;
using std::vector;
using std::array;

namespace goat {

template<size_t K, size_t D>
GR<K, D> SIMDEin(const vector<Z2k<K>>& a){
    const array<size_t, 3> set_I = {0, 1}; // index sets for SIMD Ein
    const size_t simd_delta = (D+1)/4+1; //dividing truncates, adding 1 before dividing rounds correctly

    ztk::gr_coeff<K, D> coeffs;
    for(size_t i = 0; i < simd_delta; i++){
        coeffs[set_I[i]] = a[i];
    }
    return GR<K,D>{coeffs};
}

template<size_t K, size_t D>
GR<K, D> SIMDEout(const vector<Z2k<K>>& a){
    const array<size_t, 3> set_J = {0, 2}; // index sets for SIMD Eout
    const size_t simd_delta = (D+1)/4+1; //dividing truncates, adding 1 before dividing rounds correctly

    ztk::gr_coeff<K, D> coeffs;

    size_t j = 0;

    // fill all the elements not in J will random elements
    for(size_t i = 0; i < D; i++){
        if(j >= simd_delta || i != set_J[j]){
            coeffs[i] = randomize<Z2k<K>>();
        }
        else{
            coeffs[i] = a[j];
            j++;
        }
    }
    return GR<K,D>{coeffs};
}

template<size_t K, size_t D>
GR<K, D> SIMDEout_const(const vector<Z2k<K>>& a){
    const array<size_t, 3> set_J = {0, 2}; // index sets for SIMD Eout
    const size_t simd_delta = (D+1)/4+1; //dividing truncates, adding 1 before dividing rounds correctly

    ztk::gr_coeff<K, D> coeffs;

    // fill all the elements not in J will random elements
    for(size_t i = 0; i < simd_delta; i++){
        coeffs[set_J[i]] = a[i];
    }
    return GR<K,D>{coeffs};
}

template<size_t K, size_t D>
vector<GR<K, D>> SIMDEout2Ein(const vector<GR<K, D>>& x_out){
    const array<size_t, 3> set_I = {0, 1}; // index sets for SIMD Ein
    const array<size_t, 3> set_J = {0, 2}; // index sets for SIMD Eout
    const size_t simd_delta = (D+1)/4+1; //dividing truncates, adding 1 before dividing rounds correctly

    vector<GR<K,D>> x_in(x_out.size());

    for(size_t j = 0; j < x_out.size(); j++){
        ztk::gr_coeff<K, D> coeffs;
        for(size_t i = 0; i < simd_delta; i++){
            coeffs[set_I[i]] = x_out[j][set_J[i]];
        }
        x_in[j] = GR<K,D>{coeffs};
    }

    return x_in;
}

template<size_t K, size_t D>
bool is_I_encoding(const GR<K,D>&x){
    const array<size_t, 3> set_I = {0, 1, 3}; // index sets for SIMD Ein
    const size_t simd_delta = (D+1)/4+1; //dividing truncates, adding 1 before dividing rounds correctly

    size_t p = 0;
    for(size_t k=0; k < D; k++){
        if(p < simd_delta && set_I[p]==k){
            p++;
        } else if(x[k]!=Z2k<K>::zero){
            return false;
        }
    }
    return true;
}

template<size_t K, size_t D>
bool SIMDencodings_match(const GR<K,D>& x_in, const GR<K,D>& x_out){
    const array<size_t, 3> set_I = {0, 1, 3}; // index sets for SIMD Ein
    const array<size_t, 3> set_J = {0, 2, 6}; // index sets for SIMD Eout
    const size_t simd_delta = (D+1)/4+1; //dividing truncates, adding 1 before dividing rounds correctly

    if(!is_I_encoding<K,D>(x_in)){
#ifdef TESTING
        std::cout << "xs_ein is not an I-encoding\n";
#endif
        return false;
    }

    for(size_t k = 0; k < simd_delta; k++){
        if(x_in[set_I[k]]!= x_out[set_J[k]]){
#ifdef TESTING
            std::cout << "preimages of xs_ein and xs_eout do not match\n";
#endif
            return false;
        }
    }

    return true;
}

template<size_t K, size_t D>
bool are_zero_encodings(const vector<GR<K,D>>& x_eout){
    const array<size_t, 3> set_J = {0, 2, 6};
    const size_t simd_delta = (D+1)/4+1;

    const size_t num_elem = x_eout.size();

    for(size_t i = 0; i < num_elem; i++){
        for(size_t k = 0; k < simd_delta; k++){
            if(x_eout[i][set_J[k]]!=Z2k<K>::zero){
                return false;
            }
        }
    }

    return true;
}
}

#endif // _GOAT_ENCODING_HPP
