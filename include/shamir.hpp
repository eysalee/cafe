#ifndef _SHAMIR_HPP
#define _SHAMIR_HPP

#include "ztk.hpp"
#include "ncomm.hpp"
#include "crypto.hpp"
#include "matrix.hpp"
#include "timing.hpp"
#include "f2d.hpp"
#include "encoding.hpp"

#include <vector>
#include <array>

#include <iostream>

namespace goat {

using ncomm::Network;
using ztk::GR;
using ztk::Z2k;
using std::vector;
using std::array;
using goat::StatefulPRG;

#define SAME_SIZE(a, b) assert((a).size() == (b).size())

template<size_t K, size_t D>
class ShamirProtocol {
public:

    // This is a bit confusing. For share_t and share_2t the 't' denotes the
    // degree. But for double_share_t the 't' is just to say that it's a type.
    typedef GR<K, D> share_t;
    typedef GR<K, D> share_2t;
    typedef array<GR<K, D>, 2> double_share_t;
    typedef array<GR<K,D>, D+1> flex_double_share_t;
    typedef array<GR<K,D>, 2> innerprod_double_share_t;
    typedef array<GR<K,D>, 2> simd_double_share_t;

    // residue filed types
    typedef F2d<D> fshare_t;
    typedef F2d<D> fshare_2t;

    size_t number_of_parties() const {
	return n;
    };

    size_t threshold() const {
	return t;
    };

    ShamirProtocol(Network& network) :
	_network{network}, n{network.size()}, t{(n-1)/3}, stat_sec{40}
	{
	    // note: outside of benchmarking setting this should be replaced with coin tossing
        vector <unsigned char> seed(PRG_KEY_SIZE, 42);
#ifdef TESTING
	std::cout << "seed: ";
	for (const auto &b: seed)
	    std::cout << b << " ";
	std::cout << "\n";
#endif
        _prg = new StatefulPRG(seed.data());
	    init_excep_seq();
        init_matrices();

// Note: This seems to mess with timings of local computations
#ifdef WITH_SYNC
        sync();
#endif
	};

    size_t get_n() const { return n; };
    size_t get_t() const { return t; };

    size_t get_number_of_dshares() const { return _dshares.size(); };
    size_t get_number_of_flex_dshares() const { return _flexdshares.size(); };
    size_t get_number_of_innerprod_dshares() const { return _innerproddshares.size(); };
    size_t get_number_of_simd_dshares() const { return _simddshares.size(); };

    vector<share_t> input(const vector<int>& inputs);

    // helper function. calls recon_publ
    vector<GR<K, D>> open(const vector<GR<K, D>>& shares, const size_t d);

    inline share_t add(const share_t& x, const share_t& y) {
	return x + y;
    };

    inline fshare_t add(const fshare_t& x, const fshare_t& y) {
	return x + y;
    }

    inline share_t sub(const share_t& x, const share_t& y) {
	return x - y;
    };

    inline fshare_t sub(const fshare_t& x, const fshare_t& y) {
	return x - y;
    }

    inline share_2t mul(const share_t& x, const share_t& y) {
	// NB: this doubles the degree of the share
	return x * y;
    };

    inline fshare_2t mul(const fshare_t& x, const fshare_t& y) {
	return x * y;
    };

    inline share_2t dotp(const vector<share_t>& xs, const vector<share_t>& ys) {
	SAME_SIZE(xs, ys);
	share_2t r = mul(xs[0], ys[0]);
	for (size_t i = 1; i < xs.size(); i++)
	    r = r + mul(xs[i], ys[i]);
	return r;
    };

    inline vector<vector<share_2t>> matmul(const vector<vector<share_t>>& x, const vector<vector<share_t>>& y) {
	const auto n = x.size();
	const auto m = x[0].size();
	const auto m_ = y.size();
	const auto p = y[0].size();

	assert (m == m_);

	vector<vector<share_2t>> z (n);
	for (size_t i = 0; i < n; i++) {
	    z[i].resize(p);
	    for (size_t j = 0; j < p; j++) {
		z[i][j] = x[i][0]*y[0][j];
		for (size_t k = 1; k < m; k++) {
		    z[i][j] += x[i][k]*y[k][j];
		}
	    }
	}

	return z;
    }

    vector<share_t> add(const vector<share_t>& xs, const vector<share_t>& ys) {
	SAME_SIZE(xs, ys);

	vector<share_t> r;
	r.reserve(xs.size());
	for (size_t i = 0; i < xs.size(); i++)
	    r.emplace_back(add(xs[i], ys[i]));
	return r;
    };

    vector<share_t> sub(const vector<share_t>& xs, const vector<share_t>& ys) {
	SAME_SIZE(xs, ys);
	vector<share_t> r;
	r.reserve(xs.size());
	for (size_t i = 0; i < xs.size; i++)
	    r.emplace_back(sub(xs[i], ys[i]));
	return r;
    };

    // r : number of rows in output
    //      (allows us to reuse matrix for different size inputs)
    vector<share_t> matrix_mult(const vector<vector<GR<K,D>>>& mat, const vector<share_t>& xs, const size_t r){
        assert(mat.size() > 0);
        assert(mat.size() >= r);
        assert(xs.size() > 0);
        assert(mat[0].size() >= xs.size());

        vector<share_t> ys(r);

        for(size_t i = 0; i < r; i++){
            ys[i] = mat[i][0] * xs[0];
            for(size_t j = 1; j < xs.size(); j++){
                ys[i] += mat[i][j] * xs[j];
            }
        }

        return ys;
    };

    // This is when we interpret M as Z2k and do special mults
    vector<share_t> matrix_mult_special(const vector<vector<GR<K,D>>>& mat, const vector<std::array<share_t, D>>& xs, const size_t r){

	static_assert(D == 4);
        // if(D!=4){
        //     throw std::runtime_error("Only support this particular matrix mult for D=4");
        // }

        const size_t c = xs.size();

        assert(mat.size() > 0);
        assert(mat.size() >= r);
        assert(c > 0);
        assert(mat[0].size() >= c);

        vector<share_t> ys;

        for(size_t i = 0; i < r; i++){
            std::array<GR<K,4>, 4> y4 = gr_deg4_vec_mul<K>(mat[i][0], xs[0]);
            for(size_t k = 0; k < 4; k++){
                ys.emplace_back(y4[k]);
            }

            for(size_t j = 1; j < c; j++){
                y4 = gr_deg4_vec_mul<K>(mat[i][j], xs[j]);

                for(size_t k = 0; k < 4; k++){
                    ys[i*4 + k] += y4[k];
                }
            }
        }

        return ys;
    }

    vector<share_t> reduce_degree_regular(const vector<share_2t>& xs);
    vector<share_t> reduce_degree_simd(const vector<share_2t>& xs);
    vector<share_t> reduce_degree_flex(const vector<share_2t>& xs);
    vector<share_t> reduce_degree_innerprod(const vector<share_2t>& xs);

    // naive double-shares from Abspoel et al TCC 2019 paper
    void preproc_dshares_regular(const size_t m);
    bool check_dshares_regular(const vector<double_share_t>& x);

    void preproc_dshares_flex(const size_t m);
    bool check_dshares_flex(vector<flex_double_share_t>& x);

    // void preproc_dshares_innerprod(const size_t m);
    // bool check_dshares_innerprod(vector<innerprod_double_share_t>& x);

    void preproc_dshares_simd(const size_t m);
    bool check_dshares_simd(vector<simd_double_share_t>& x);

    vector<share_t> rand_bits_simd(const size_t m);
    vector<GR<K,D>> gen_simd_const_shares(const size_t m, const vector<Z2k<K>> a);
    vector<share_t> gen_simd_zero_shares(const size_t m);
    bool check_simd_zero_shares(vector<share_t>& x);

    // non-robust broadcast n-t elements
    void broadcast(const unsigned int send_id, vector<GR<K,D>>& xs);

    void recon_priv(GR<K,D> share, const size_t d, const unsigned int recv_id, GR<K,D>& secret);

    vector<GR<K,D>> recon_publ(vector<GR<K,D>> &shares, const size_t d);

    // wait for parties that are still connecting
    void sync(){
    TIMER_BEGIN(sync);

    const size_t last_party = n-1;
    const size_t bool_size = sizeof(bool);

    if(_network.id()==last_party){
        bool sync_bit = true;
        vector<unsigned char> bit_sbuf(bool_size);
        memcpy(bit_sbuf.data(), &sync_bit, bool_size);

        _network.broadcast_send(bit_sbuf);
    } else {
        vector<unsigned char> bit_rbuf(bool_size);

        _network.broadcast_recv(last_party, bit_rbuf);
    }

    TIMER_END(sync);
    };

#ifdef TESTING
    void test();

    void test_bcast();

    void test_flex_lincombo();

    void test_flex_shares();

    vector<share_t> test_make_share();

    void test_prg();

    void test_simd();

    void test_deg_reduce();

#endif

    inline vector<share_t> share(const GR<K,D> v, const size_t d) {
	// generate a sharing of v

    assert(d < n);

	vector<share_t> points (d+1);
	vector<share_t> shares (n);
	points[0] = v;

	// This is polynomial evaluation via. Horners method.
	for (size_t i = 1; i < d+1; i++)
	    points[i] = randomize<share_t>();

	for (size_t i = 0; i < n; i++) {
	    const auto x = excep_seq[i];
	    auto px = points[d] * x + points[d-1];

	    // note: can't use size_t because size_t doesn't have negative values
	    for (ssize_t j = d-2; j >= 0; j--)
		px = px * x + points[j];

	    shares[i] = px;
	}

	return shares;
    };

private:

    // reconstructs polynomial and only returns the "secret" (polynomial evaluated at 0)
    // calling reconstruct with only 2 parameters checks all n shares are consistent
    GR<K, D> reconstruct(const vector<GR<K, D>> &shares, const size_t d){ return reconstruct(shares, d, n); };

    // reconstruct with a variable number of shares to check are consistent
    GR<K, D> reconstruct(const vector<GR<K, D>> &shares, const size_t d, const size_t num_check);

    // reconstructs polynomial and returns all of the coefficients
    vector<GR<K,D>> recon_coeffs(const vector<GR<K,D>> &shares, const size_t d, const size_t num_check);

    // helper functions for sharing
    inline vector<share_t> share(const int v, const size_t d){ return share(GR<K,D>{Z2k<K>{v}}, d); };
    inline vector<share_t> share(const Z2k<K> v, const size_t d){ return share(GR<K, D>{v}, d); };


    void init_excep_seq_d4(vector<GR<K, 4>> &seq);

    void init_excep_seq() {
	if (D == 4){
	    init_excep_seq_d4(excep_seq);
    }
	else
	    throw std::runtime_error("only support D=4 exceptional sequences");
    };

    void init_matrices(){
        init_hyper_matrix<K,D>(excep_seq, n, n, n, hyper_matrix);
        init_vand_matrix<K,D>(excep_seq, n, n, vand_matrix);
        init_inv_vand_matrix<K, D>(vand_matrix, n-t, inv_vand_matrix);
    };

    Network _network;

    size_t n;  // number of parties
    size_t t;  // threshold
    size_t stat_sec; // statistical security parameter

    // degree t items in position 0 and degree 2t in position 1
    vector<double_share_t> _dshares; // preproc double shares

    // t-share will go in positions 0...D-1 and 2t-share in D position
    vector<flex_double_share_t> _flexdshares;

    vector<simd_double_share_t> _simddshares;

    vector<innerprod_double_share_t> _innerproddshares;

    // keeps track of which preprocessed shares we've used
    // (not a good way to pop things off front of vector as far as I know)
    size_t _dcntr;
    size_t _fdcntr;
    size_t _sdcntr;
    size_t _ipdcntr;

    vector<GR<K, D>> excep_seq;

    // hyper-invertible matrix
    vector<vector<GR<K, D>>> hyper_matrix;
    // vandermonde matrix
    vector<vector<GR<K,D>>> vand_matrix;
    // inverse of a vandermonde matrix; for computing coefficients
    vector<vector<GR<K,D>>> inv_vand_matrix;

    // Stateful PRG
    StatefulPRG *_prg;
};

template<size_t K, size_t D>
void ShamirProtocol<K,D>::init_excep_seq_d4(vector<GR<K, 4>> &seq) {
    assert(n <= 15);
    // Hard coded to 16 since we'll possibly need more all 2^4 points
    seq.resize(16);

#define ADD(idx, i0, i1, i2, i3) do {                   \
    seq[idx] = GR<K, 4>{{(int)i3, (int)i2, (int)i1, (int)i0}};  \
    if(idx >= 15){ \
        goto done; \
    } \
    } while(0)

    ADD(0,  0, 0, 0, 1);
    ADD(1,  0, 0, 1, 0);
    ADD(2,  0, 0, 1, 1);
    ADD(3,  0, 1, 0, 0);
    ADD(4,  0, 1, 0, 1);
    ADD(5,  0, 1, 1, 0);
    ADD(6,  0, 1, 1, 1);
    ADD(7,  1, 0, 0, 0);
    ADD(8,  1, 0, 0, 1);
    ADD(9,  1, 0, 1, 0);
    ADD(10, 1, 0, 1, 1);
    ADD(11, 1, 1, 0, 0);
    ADD(12, 1, 1, 0, 1);
    ADD(13, 1, 1, 1, 0);
    ADD(14, 1, 1, 1, 1);

#undef ADD

done:
    return;

}

// create preprocessed double shares
// naive double-shares given in the Abspoel et al TCC 2019 paper
template<size_t K, size_t D>
void ShamirProtocol<K,D>::preproc_dshares_regular(const size_t m) {
    // m : the number of double shares each party should create

    // 1. Sample m random values
    // 2. Create shares of degree t and 2t for each val
    // 3. Exchange with everyone
    // 4. save into _dshares


    TIMER_BEGIN(ds);

    TIMER_BEGIN(gen);

    const size_t num_batches = m / (D*t);
    const size_t share_size = GR<K, D>::byte_size();

    std::cout << "num_batches: " << num_batches << "\n";

    vector<vector<unsigned char>> sbufs (n);
    vector<vector<unsigned char>> rbufs (n);

    for (size_t i = 0; i < n; i++) {
	sbufs[i].resize( 2 * D * num_batches*share_size );
	rbufs[i].resize( 2 * D * num_batches*share_size );
    }

    TIMER_BEGIN(local1);

// #pragma omp parallel for
    for (size_t i = 0; i < num_batches; i++) {
	const auto x = randomize<GR<K,D>>();

	vector<vector<share_t>> s_t (D);
	vector<vector<share_2t>> s_2t (D);

	for (size_t k = 0; k < D; k++) {
	    s_t[k] = share(x[k], t);
	    s_2t[k] = share(x[k], 2*t);
	}

	for (size_t j = 0; j < n; j++) {
	    for (size_t k = 0; k < D; k++) {
		s_t[k][j].pack(sbufs[j].data() + ((2*D*i + k)*share_size));
		s_2t[k][j].pack(sbufs[j].data() + ((2*D*i + D + k)*share_size));
	    }
	}
    }

    TIMER_END(local1);

    _network.exchange_all(sbufs, rbufs);

    TIMER_BEGIN(local2);

    vector<double_share_t> shares_to_check ((n-t)*D*num_batches);
    _dshares.resize(t*D*num_batches);

    for (size_t i = 0; i < num_batches; i++) {
	vector<std::array<share_t, D>> x_t (n);
	vector<std::array<share_2t, D>> x_2t (n);

	for (size_t j = 0; j < n; j++) {
	    for (size_t k = 0; k < D; k++) {
		x_t[j][k] = GR<K,D>{rbufs[j].data() + ((2*D*i + k)*share_size)};
		x_2t[j][k] = GR<K,D>{rbufs[j].data() + ((2*D*i + k + D)*share_size)};
	    }
	}

	vector<share_t> r_t = matrix_mult_special(hyper_matrix, x_t, n);
	vector<share_2t> r_2t = matrix_mult_special(hyper_matrix, x_2t, n);

	for (size_t j = 0; j < t; j++) {
	    for (size_t k = 0; k < D; k++) {
		_dshares[i*t*D + k] = {r_t[j*D + k], r_2t[j*D + k]};
	    }
	}
	for (size_t j = t; j < n; j++) {
	    const auto off = j - t;
	    for (size_t k = 0; k < D; k++) {
		shares_to_check[i*off*D + k] = {r_t[j*D + k], r_2t[j*D + k]};
	    }
	}
    }


    TIMER_END(local2);
    TIMER_END(gen);

    TIMER_BEGIN(check);


    if (!check_dshares_regular(shares_to_check)) {
#ifdef TESTING
	std::cout << "share check failed\n";
#else
    	throw std::runtime_error("checking shares failed");
#endif
    }

    TIMER_END(check);

    _dcntr = 0;

    TIMER_END(ds);
}

template<size_t K, size_t D>
bool ShamirProtocol<K, D>::check_dshares_regular(const vector<double_share_t>& x) {

    const size_t batch_size = n-t;
    const auto m = x.size()/batch_size;
    
    // For each n shares, construct i = 1,...,n-t to Pi. Pi verifies the secret is correct
    // y and y2 will be empty for parties n-t,...,n
    vector<GR<K,D>> y(m);
    vector<GR<K,D>> y2(m);

    for(size_t j = 0; j < m; j++){
        for(size_t i = 0; i < batch_size; i++){
            // Party i has double share privately reconstructed to them
            //   and stores the values in y[j] and y2[j]
            recon_priv(x[j*batch_size + i][0], t, i, y[j]);
            recon_priv(x[j*batch_size + i][1], 2*t, i, y2[j]);
        }
    }

    bool is_good = true;
    // only parties 0,...,n-t-1 will have had values reconstructed to them, so they are the only ones that need to check
    if(_network.id() < batch_size){
        for(size_t j = 0; j < m; j++){
            // Check each double shared value is (1) Z2k element and (2) matches
            for(size_t k = 1; k < D; k++){
                if(y[j][k]!=Z2k<K>::zero || y2[j][k]!=Z2k<K>::zero){
                    is_good = false;
                    goto endloop_dshare;
                }
            }

            if(y[j][0]!=y2[j][0]){
                is_good = false;
                break;
            }
        }
    endloop_dshare:

        // send the is_ok bit to everyone
        vector<unsigned char> bit_sbuf(sizeof(bool));
        memcpy(bit_sbuf.data(), &is_good, sizeof(bool));

    #ifdef TESTING
        std::cout << "sending bit " << is_good << "\n";
    #endif

        _network.broadcast_send(bit_sbuf);
    }

    vector<vector<unsigned char>> bit_rbufs(batch_size);

    for(size_t i = 0; i < batch_size; i++){
        bit_rbufs[i].resize(sizeof(bool));
        _network.broadcast_recv(i, bit_rbufs[i]);

        bool r_bit;
        memcpy(&r_bit, bit_rbufs[i].data(), sizeof(bool));

#ifdef TESTING
        std::cout << "bit received from party "<< i <<": " << r_bit << "\n";
#endif

        is_good = is_good && r_bit;
    }

    // indicates if everyone was happy or not
    return is_good;
}

// FLEX double shares
template<size_t K, size_t D>
void ShamirProtocol<K,D>::preproc_dshares_flex(const size_t m){
    // note: m will be the total number of shares
    //      we make batches of size n-t, so m needs to be divisble by this

    const size_t share_size = GR<K, D>::byte_size();
    const size_t batch_size = n-t;

    assert(m % batch_size == 0);

    const size_t num_batches = m / batch_size;

    vector<vector<unsigned char>> sbufs(n);
    vector<vector<unsigned char>> rbufs(n);

    for(size_t i = 0; i < n; i++){
        sbufs[i].resize(num_batches * share_size * (D + 1));
        rbufs[i].resize(num_batches * share_size * (D + 1));
    }

    _flexdshares.resize(m);

    // vector<vector<share_t>> s_t (D);
    // vector<share_2t> s_2t (n);

    TIMER_BEGIN(ds);

    TIMER_BEGIN(gen);
    TIMER_BEGIN(local1);
    // i : batch

    vector<GR<K,D>> rands (num_batches);
#pragma omp parallel for
    for (size_t i = 0; i < num_batches; i++)
	rands[i] = randomize<GR<K,D>>();

#pragma omp parallel for
    for(size_t i = 0; i < num_batches; i++){

	vector<vector<share_t>> s_t (D);
	vector<share_2t> s_2t (n);

        // std::array<Z2k<K>, D> s_coeffs;

	// const auto x = randomize<GR<K,D>>();

        for(size_t j = 0; j < D; j++){
            // s_coeffs[j] = randomize<Z2k<K>>();
            s_t[j] = share(rands[i][j], t);
        }

        s_2t = share(rands[i], 2*t);

        for(size_t j = 0; j < n; j++){
            for(size_t k = 0; k < D; k++){
                s_t[k][j].pack(sbufs[j].data() + ((D + 1)*i + k)*share_size);
            }
            s_2t[j].pack(sbufs[j].data() + ((D + 1)*i + D)*share_size);
        }
    }

    TIMER_END(local1);

    _network.exchange_all(sbufs, rbufs);

    TIMER_BEGIN(local2);
    // vector<std::array<share_t, D>> x_t(n);
    // vector<share_2t> x_2t(n);

// #pragma omp parallel for
    for(size_t i = 0; i < num_batches; i++){

	vector<std::array<share_t, D>> x_t(n);
	vector<share_2t> x_2t(n);

        for(size_t j = 0; j < n; j++){
            for(size_t k = 0; k < D; k++){
                x_t[j][k] = GR<K,D>{ rbufs[j].data() + ((D + 1)*i + k)*share_size };
            }
            x_2t[j] = GR<K,D>{ rbufs[j].data() + ((D + 1)*i + D)*share_size};
        }

        vector<share_2t> r_2t = matrix_mult(hyper_matrix, x_2t, batch_size);
        vector<share_t> r_t = matrix_mult_special(hyper_matrix, x_t, batch_size);

        for(size_t j = 0; j < batch_size; j++){
            for(size_t k = 0; k < D; k++){
                _flexdshares[i * batch_size + j][k] = r_t[j*D + k];
            }
            _flexdshares[i * batch_size + j][D] = r_2t[j];
        }
    }

    TIMER_END(local2);
    TIMER_END(gen);

    TIMER_BEGIN(check);

    if(check_dshares_flex(_flexdshares)!=true){
        throw std::runtime_error("inconsistent flex dshare check!");
    }

    TIMER_END(check);

    TIMER_END(ds);

    _fdcntr = 0;
}

template<size_t K, size_t D>
bool ShamirProtocol<K, D>::check_dshares_flex(vector<flex_double_share_t>& x) {
    // 1. check each <r^i>t is a sharing of a Z2k element
    // 2. check that each double sharing is well formed
    // 3. throw out double shares that acted as masks

    const size_t batch_size = n - t;
    const size_t num_shares = x.size();
    assert(num_shares % batch_size == 0);


    assert(D==4); // currently only works for D = 4
    assert(num_shares > stat_sec);
    // assert(stat_sec % batch_size == 0); // can currently only handle this case (because of how we're using reconspubl)
    assert(stat_sec % D == 0);

    const size_t num_mask_shares = stat_sec/D; // going to throw out last stat_sec/D shares
    const size_t num_checked_shares = num_shares - num_mask_shares;
    const size_t num_rand_bits = num_checked_shares*stat_sec;

    // NOTE: IF YOU TRY TO CHECK TOO MANY SHARES AT ONCE, THIS WILL SEGFAULT
    // It's because arrays are allocated on the stack, and at compile time.
    // unsigned char rand_bits[num_rand_bits];
    unsigned char *rand_bits = (unsigned char *)malloc(num_rand_bits);

    if (!rand_bits)
    throw std::runtime_error("malloc error");

    _prg->next(rand_bits, num_rand_bits);

    // Doing vector<vector<share_t>> makes it easier to call recon_publ later
    vector<vector<share_t>> xt(stat_sec/batch_size);
    for(size_t ii = 0; ii < stat_sec/batch_size; ii++){
        xt[ii].resize(batch_size);
        for(size_t jj = 0; jj < batch_size; jj++){
            xt[ii][jj] = GR<K,D>::zero();

            for(size_t i = 0; i < num_checked_shares; i++){
                for(size_t k = 0; k < D; k++){
                    unsigned char bit = (1 & (rand_bits[(ii*batch_size + jj)*num_checked_shares + i] >> k));
                    if(bit==1){
                        xt[ii][jj] += x[i][k];
                    }
                }
            }

            // hardcode these ones to always be added (their rand_bits are always 1)
            for(size_t i = num_checked_shares; i < num_shares; i++){
                for(size_t k = 0; k < D; k++){
                    xt[ii][jj] += x[i][k];
                }
            }
        }
    }

    free(rand_bits);

    for(size_t i = 0; i < xt.size(); i++){
        vector<GR<K,D>> xs = recon_publ(xt[i], t);
        // check each xs is in Z2k
        for(size_t j = 0; j < batch_size; j++){
            for(size_t k = 1; k < D; k++){
                // checking it's a Z2k element by checking all other coefficients are zero
                if(xs[j][k] != 0){
                    // throw std::runtime_error("xt in flex dshare check is not element in Z2k");
#ifdef TESTING
                    std::cout << "xt in flex dshare check is not element in Z2k\n";
#endif
                    return false;
                }
            }
        }
    }

    //  can throw the last lambda/d elements if we hardcode those to always be 1
    x.erase(x.end() - num_mask_shares, x.end());

    // check correct relationship holds for the remaining double shares
    // this we need tau*d > stat_sec
    const size_t tau = batch_size + (stat_sec/D);

    const size_t num_rand_seq = tau*num_checked_shares;
    unsigned char rand_buf[num_rand_seq];
    vector<GR<K,D>> rand_seq(num_rand_seq);

    _prg->next(rand_buf, num_rand_seq);
    // convert random values to random excep_seq elements
    for(size_t i = 0; i < num_rand_seq; i++){
        rand_seq[i] = excep_seq[(0b111 & rand_buf[i])]; // todo: lol double check
    }

    vector<vector<share_2t>> yt(tau/batch_size);

    array<int, 4> powi = {0, 1, 3, 7};

    for(size_t i = 0; i < tau/batch_size; i++){
        yt[i].resize(batch_size);
        for(size_t ii = 0; ii < batch_size; ii++){
            yt[i][ii] = GR<K,D>::zero();

            for(size_t j = 0; j < num_checked_shares; j++){
                share_2t r = x[j][D];

                for(size_t k = 0; k < D; k++){
                    r -= x[j][k] * excep_seq[powi[k]];
                }

                yt[i][ii] += rand_seq[(i * batch_size + ii) * num_checked_shares + j] * r;
            }
        }
    }

    for(size_t i = 0; i < tau/batch_size; i++){
        vector<GR<K,D>> ys = recon_publ(yt[i], 2*t);
        for(size_t j = 0; j < batch_size; j++){
            if(ys[j]!= GR<K,D>::zero()){
                // std::runtime_error("yt in flex dshare check is non-zero (double share is not of the correct form)");
#ifdef TESTING
                std::cout << "yt in flex dshare check is non-zero (double share is not of the correct form)\n";
#endif
                return false;
            }
        }
    }

#ifdef TESTING
    std::cout << "flex dshares passed check!\n";
#endif
    return true;
}

// template<size_t K, size_t D, char Direction>
// static inline GR<K, D> Ein(const vector<Z2k<K>>& a) {

//     // delta < (D + 1)/ 2 ===> D = 4 means we can have delta = 2
//     // delta < D - Dtilde + 1 ===> H(X) = X^4 + X + 1, so D - Dtile + 1 = 4 - 1 + 1 = 4

//     constexpr size_t delta = D/2;

//     static_assert (Direction == 'L' || Direction == 'R');

//     if(Direction == 'L') {
// 	ztk::gr_coeff<K, D> r;
// 	for (size_t i = 0; i < delta; i++)
// 	    r[i] = a[i];
// 	return GR<K, D>{r};
//     } else if (Direction == 'R') {
// 	ztk::gr_coeff<K, D> r;
// 	for (size_t i = 0; i < delta; i++)
// 	    r[D+i-1] = a[i];
// 	return GR<K, D>{r};
//     }
// }

// template<size_t K, size_t D>
// static inline GR<K, D> Eout(const Z2k<K>& c) {
//     ztk::gr_coeff<K,D> x;
//     x[D-1] = c;
//     return GR<K,D>{x};
// }

// inner product -> FLEX double shares.
// template<size_t K, size_t D>
// void ShamirProtocol<K,D>::preproc_dshares_innerprod(const size_t m){

//     // a lot of this stuff follows the same general layout as with FLEX
//     const size_t share_size = GR<K, D>::byte_size();
//     const size_t batch_size = (n-t);
//     const size_t num_batches = m / batch_size;

//     std::cout << "num_batches=" << num_batches << "\n";

//     typedef unsigned char u8;

//     vector<vector<u8>> sbufs(n);
//     vector<vector<u8>> rbufs(n);

//     for (size_t i = 0; i < n; i++) {
//     	sbufs[i].resize(num_batches * share_size * 2*D);
//     	rbufs[i].resize(num_batches * share_size * 2*D);
//     }

//     // For each batch, sample D random elements and share them with degree t and
//     // 2*t.

//     vector<vector<share_t>> s_t (D);
//     vector<vector<share_2t>> s_2t (D);

//     TIMER_BEGIN(gen);
//     TIMER_BEGIN(local1);

//     for (size_t i = 0; i < num_batches; i++) {
//     	for (size_t j = 0; j < D; j++) {
//     	    const auto sij = randomize<Z2k<K>>();
//     	    s_t[j] = share(sij, t);
//     	    s_2t[j] = share(Eout<K,D>(sij), 2*t);
//     	}

// 	for (size_t j = 0; j < n; j++) {
// 	    for (size_t k = 0; k < D; k++) {
// 		s_t[k][j].pack(sbufs[j].data() + (2*D*i + k) * share_size);
// 		s_2t[k][j].pack(sbufs[j].data() + (2*D*i + D + k) * share_size);
// 	    }
// 	}
//     }
//     TIMER_END(local1);

//     TIMER_BEGIN(exchange1);
//     _network.exchange_all(sbufs, rbufs);
//     TIMER_END(exchange1);

//     _innerproddshares.resize(D*batch_size*num_batches);

//     vector<std::array<share_t, D>> x_t (n);
//     vector<vector<share_2t>> x_2t (D);
//     vector<vector<share_2t>> r_2t (D);

//     for (size_t i = 0; i < D; i++) {
// 	x_2t[i].resize(n);
// 	r_2t[i].resize(n);

//     }

//     TIMER_BEGIN(local2);

//     for (size_t i = 0; i < num_batches; i++) {
//     	for (size_t k = 0; k < D; k++) {
//     	    for (size_t j = 0; j < n; j++) {
//     		x_2t[k][j] = GR<K,D>{rbufs[j].data() + (2*D*i + D + k) * share_size};
//     		x_t[j][k] = GR<K,D>{rbufs[j].data() + (2*D*i + k)*share_size};
//     	    }
//     	    r_2t[k] = matrix_mult(hyper_matrix, x_2t[k], batch_size);
//     	}
//     	vector<share_t> r_t = matrix_mult_special(hyper_matrix, x_t, batch_size);

// 	const size_t offset = i * batch_size * D;

//     	for (size_t ii = 0; ii < batch_size; ii++) {
//     	    for (size_t kk = 0; kk < D; kk++) {
//     		_innerproddshares[offset + ii*D + kk][0] = r_t[ii*D + kk];
// 		_innerproddshares[offset + ii*D + kk][1] = r_2t[kk][ii];
//     	    }
//     	}
//     }
//     TIMER_END(local2);

//     TIMER_END(gen);

//     // TIMER_BEGIN(check);
//     // if (check_dshares_innerprod(_innerproddshares)) {
//     // 	throw std::runtime_error("innerprod check failed");
//     // }
//     // TIMER_END(check);

//     _ipdcntr = 0;
// }

// template<size_t K, size_t D>
// bool ShamirProtocol<K,D>::check_dshares_innerprod(vector<innerprod_double_share_t>& x) {
//     (void)x;
//     return false;
// }

// SIMD double shares
template<size_t K, size_t D>
void ShamirProtocol<K,D>::preproc_dshares_simd(const size_t m){

    TIMER_BEGIN(ds);

    TIMER_BEGIN(simd_gen);
    // 1. sample d arrays of delta Z2k elements and share the SIMDin and SIMDout encodings of the value
    // 2. apply M to the degree t and degree 2t sharings

    const size_t simd_delta = (D+1)/4+1;

    const size_t batch_size = D*(n-t);
    assert(m % batch_size == 0);
    const size_t num_batches = m / batch_size;

    const size_t share_size = GR<K, D>::byte_size();
    vector<vector<unsigned char>> sbufs(n);
    vector<vector<unsigned char>> rbufs(n);

    for(size_t i = 0; i < n; i++){
        sbufs[i].resize(num_batches * share_size * D * 2);
        rbufs[i].resize(num_batches * share_size * D * 2);
    }

    TIMER_BEGIN(simd_local1);
#pragma omp parallel for
    for(size_t k = 0; k < num_batches; k++){

	    vector<share_t> s_t(n);
	    vector<share_2t> s_2t(n);
	    vector<Z2k<K>> s(simd_delta);

        for(size_t i = 0; i < D; i++){
            for(size_t j = 0; j < simd_delta; j++){
                s[j] = randomize<Z2k<K>>();
            }

            s_t = share(SIMDEin<K,D>(s), t);
            s_2t = share(SIMDEout<K,D>(s), 2*t);

            for(size_t j = 0; j < n; j++){
                s_t[j].pack(sbufs[j].data() + (2*(k*D + i))*share_size);
                s_2t[j].pack(sbufs[j].data() + (2*(k*D + i)+1)*share_size);
            }
        }
    }
    TIMER_END(simd_local1);

    TIMER_BEGIN(simd_exchange);
    _network.exchange_all(sbufs, rbufs);
    TIMER_END(simd_exchange);

    _simddshares.resize(m);

    vector<array<share_t, D>> s_ein(n);
    vector<array<share_2t, D>> s_eout(n);

    TIMER_BEGIN(simd_local2);
    for(size_t k = 0; k < num_batches; k++){
        for(size_t i = 0; i < n; i++){
            for(size_t j = 0; j < D; j++){
                s_ein[i][j] = GR<K,D>{ rbufs[i].data() + 2*(k*D + j)*share_size };
                s_eout[i][j] = GR<K,D>{ rbufs[i].data() + (2*(k*D + j) + 1)*share_size };
            }
        }

        vector<share_t> r_ein = matrix_mult_special(hyper_matrix, s_ein, n-t);
        vector<share_2t> r_eout = matrix_mult_special(hyper_matrix, s_eout, n-t);

        for(size_t i = 0; i < batch_size; i++){
            _simddshares[k*batch_size + i][0] = r_ein[i];
            _simddshares[k*batch_size + i][1] = r_eout[i];
        }
    }
    TIMER_END(simd_local2);
    TIMER_END(simd_gen);

    TIMER_BEGIN(simd_check);
    if(!check_dshares_simd(_simddshares)){
        #ifdef TESTING
        std::cout << "not good!\n";
        return;
        #endif
        throw std::runtime_error("inconsistent simd dshare check!");
    }
    TIMER_END(simd_check);

    TIMER_END(ds);

    _sdcntr = 0;
}

template<size_t K, size_t D>
bool ShamirProtocol<K,D>::check_dshares_simd(vector<simd_double_share_t>& x){
    const size_t num_shares = x.size();
    const size_t num_masks = stat_sec;
    const size_t num_checked = num_shares - num_masks;
    const size_t num_rand_bits = stat_sec * num_checked;

    unsigned char *rand_bits = (unsigned char *)malloc(num_rand_bits);

    if(!rand_bits)
        throw std::runtime_error("malloc error");

    _prg->next(rand_bits, num_rand_bits);

    // recon_publ expects this many elements
    const size_t batch_size = n-t;
    const size_t num_batches = stat_sec/batch_size;
    vector<vector<share_t>> y_t(num_batches);
    vector<vector<share_2t>> y_2t(num_batches);

    for(size_t ii = 0; ii < num_batches; ii++){
        y_t[ii].resize(batch_size);
        y_2t[ii].resize(batch_size);
        for(size_t jj = 0; jj < batch_size; jj++){
            y_t[ii][jj] = GR<K,D>::zero();
            y_2t[ii][jj] = GR<K,D>::zero();

            for(size_t k = 0; k < num_checked; k++){
                unsigned char bit = 1 & rand_bits[(ii * batch_size + jj)*num_checked + k];
                if(bit==1){
                    y_t[ii][jj] += x[k][0];
                    y_2t[ii][jj] += x[k][1];
                }
            }

            for(size_t k = num_checked; k < num_shares; k++){
                y_t[ii][jj] += x[k][0];
                y_2t[ii][jj] += x[k][1];
            }
        }
    }

    free(rand_bits);

    for(size_t i = 0; i < num_batches; i++){
        vector<GR<K,D>> xs_ein = recon_publ(y_t[i], t);
        vector<GR<K,D>> xs_eout = recon_publ(y_2t[i], 2*t);
        for(size_t j = 0; j < batch_size; j++){
            // check preimage of each xs_ein and xs_eout match
            // this also checks xs_ein is an I-encoding
            if(!SIMDencodings_match<K,D>(xs_ein[j], xs_eout[j])){
                return false;
            }
        }
    }

    // throw out the last lambda double shares that we hardcoded to be 1
    x.erase(x.end() - num_masks, x.end());

#ifdef TESTING
    std::cout << "simd dshares check passed!\n";
#endif
    return true;
}

template<size_t K, size_t D>
auto ShamirProtocol<K,D>::rand_bits_simd(const size_t m) -> vector<share_t>{
    const size_t simd_delta = (D+1)/4+1;

    const size_t batch_size = D*(n-t);
    assert(m % batch_size == 0);
    const size_t num_batches = m / batch_size;

    const size_t share_size = GR<K, D>::byte_size();
    vector<vector<unsigned char>> sbufs(n);
    vector<vector<unsigned char>> rbufs(n);

    for(size_t i = 0; i < n; i++){
        sbufs[i].resize(num_batches * share_size * D);
        rbufs[i].resize(num_batches * share_size * D);
    }

    vector<share_t> s_t(n);
    vector<Z2k<K>> s(simd_delta);
    for(size_t k = 0; k < num_batches; k++){
        for(size_t i = 0; i < D; i++){
            for(size_t j = 0; j < simd_delta; j++){
                s[j] = randomize<Z2k<K>>();
            }

            s_t = share(SIMDEin<K,D>(s), t);

            for(size_t j = 0; j < n; j++){
                s_t[j].pack(sbufs[j].data() + (k*D + i)*share_size);
            }
        }
    }

    _network.exchange_all(sbufs, rbufs);

    vector<array<share_t, D>> s_ein(n);
    vector<share_t> us;
    us.reserve(m);

    for(size_t k = 0; k < num_batches; k++){
        for(size_t i = 0; i < n; i++){
            for(size_t j = 0; j < D; j++){
                s_ein[i][j] = GR<K,D>{ rbufs[i].data() + (k*D + j)*share_size };
            }
        }

        vector<share_t> r_ein = matrix_mult_special(hyper_matrix, s_ein, n-t);

        for(size_t i = 0; i < batch_size; i++){
            us.emplace_back(r_ein[i]);
        }
    }

    const size_t recon_batch_size = n-t;
    const vector<Z2k<K>> ones(simd_delta, Z2k<K>::one);
    GR<K,D> ein_one = SIMDEin(ones);

    vector<share_t> eout_zero_shares = gen_simd_zero_shares(m + stat_sec);

    vector<share_t> as_ein(m);
    vector<vector<share_t>> asquare_eout(m/recon_batch_size);

    for(size_t i = 0; i < m/recon_batch_size; i++){
        for(size_t j = 0; j < recon_batch_size; j++){
            as_ein[i * recon_batch_size + j] = 2*us[i * recon_batch_size + j] + ein_one;
            asquare_eout[i][j] = as_ein[i * recon_batch_size + j] * as_ein[i * recon_batch_size + j] + eout_zero_shares[i * recon_batch_size + j];
        }
    }

    // todo: implement steps opening asquare_eout
    // todo: implement finding smallest root? ohno
    // todo: implement computing d
    // todo: implement final output
    std::cout << "This is currently not fully implemented!\n";
    return eout_zero_shares;
}

template<size_t K, size_t D>
vector<GR<K,D>> ShamirProtocol<K,D>::gen_simd_const_shares(const size_t m, const vector<Z2k<K>> a){
    const size_t simd_delta = (D+1)/4+1;

    GR<K,D> A = SIMDEout_const(a);
    vector<share_t> shares = gen_simd_zero_shares(m);

    // note: the reason this is shares.size() rather than m is because we throw out some shares when we do the check
    for(size_t i = 0; i < shares.size(); i++){
        shares[i] += A;
    }

    return shares;
}

template<size_t K, size_t D>
vector<GR<K,D>> ShamirProtocol<K,D>::gen_simd_zero_shares(const size_t m){
    // generate random shares of zero
    // this is going to be very similar to how we preprocess double shares
    const size_t simd_delta = (D+1)/4+1;
    const vector<Z2k<K>> zeros(simd_delta, Z2k<K>::zero);

    const size_t batch_size = D*(n-t);
    assert(m % batch_size == 0);
    const size_t num_batches = m / batch_size;

    const size_t share_size = GR<K, D>::byte_size();
    vector<vector<unsigned char>> sbufs(n);
    vector<vector<unsigned char>> rbufs(n);

    for(size_t i = 0; i < n; i++){
        sbufs[i].resize(num_batches * share_size * D);
        rbufs[i].resize(num_batches * share_size * D);
    }

    vector<share_t> z_t(n);
    for(size_t i = 0 ; i < num_batches; i++){
        for(size_t j = 0; j < D; j++){
            z_t = share(SIMDEout<K,D>(zeros), t);

            for(size_t k = 0; k < n; k++){
                z_t[k].pack(sbufs[k].data() + (i*D + j)*share_size);
            }
        }
    }

    _network.exchange_all(sbufs, rbufs);

    vector<array<share_t, D>> shares(n);
    vector<share_t> zero_shares(m);

    for(size_t k = 0; k < num_batches; k++){
        for(size_t i = 0; i < n; i++){
            for(size_t j = 0; j < D; j++){
                shares[i][j] = GR<K,D>{ rbufs[i].data() + (k*D + j) * share_size };
            }
        }

        vector<share_t> rand_shares = matrix_mult_special(hyper_matrix, shares, n-t);
        for(size_t i = 0; i < batch_size; i++){
            zero_shares[k * batch_size + i] = rand_shares[i];
        }
    }

    if(!check_simd_zero_shares(zero_shares)){
        throw std::runtime_error("trouble making simd zero shares!");
    }

#ifdef TESTING
    std::cout << "simd zero shares passed check!\n";
#endif
    return zero_shares;
}

template<size_t K, size_t D>
bool ShamirProtocol<K,D>::check_simd_zero_shares(vector<share_t>& x){
    // check these are indeed shares of zero
    const size_t num_shares = x.size();
    const size_t num_masks = stat_sec;
    const size_t num_checked = num_shares - num_masks;
    const size_t num_rand_bits = stat_sec * num_checked;

    unsigned char *rand_bits = (unsigned char *)malloc(num_rand_bits);

    if(!rand_bits)
        throw std::runtime_error("malloc error");

    _prg->next(rand_bits, num_rand_bits);

    // recon_publ expects this many elements
    const size_t batch_size = n-t;
    const size_t num_batches = stat_sec/batch_size;
    vector<vector<share_t>> y_t(num_batches);

    for(size_t i = 0; i < num_batches; i++){
        y_t[i].resize(batch_size);
        for(size_t j = 0; j < batch_size; j++){
            y_t[i][j] = GR<K,D>::zero();

            for(size_t k = 0; k < num_checked; k++){
                unsigned char bit = 1 & rand_bits[(i * batch_size + j)*num_checked + k];
                if(bit==1){
                    y_t[i][j] += x[k];
                }
            }

            for(size_t k = num_checked; k < num_shares; k++){
                y_t[i][j] += x[k];
            }
        }
    }

    free(rand_bits);

    for(size_t i = 0; i < num_batches; i++){
        vector<GR<K,D>> xs_eout = recon_publ(y_t[i], t);
        if(!are_zero_encodings(xs_eout)){
            return false;
        }
    }

    // throw out the last lambda double shares that we hardcoded to be 1
    x.erase(x.end() - num_masks, x.end());

    return true;
}

// interpolate
template<size_t K, size_t D>
GR<K, D> ShamirProtocol<K, D>::reconstruct(const vector<GR<K, D>> &shares, const size_t d, const size_t num_check){
    // num_check says how many shares to check lie on the polynomial
    //  (at least d+1, but not more than n)
    const size_t num_shares = shares.size();

    assert (d < num_shares);
    assert (num_check <= num_shares);

    GR<K, D> secret;

#define X(i) excep_seq[i]

    for(size_t i = 0; i < d+1; i++){
        auto L = GR<K, D>::one();
        for(size_t j = 0; j < d+1; j++){
            if(i==j){
                continue;
            }
            L = L * (-X(j) / (X(i) - X(j)));
        }
        secret = secret + L * shares[i];
    }

    if(d+1 == n || num_check == d+1){
        return secret;
    }

    // Checks that shares are consistent with polynomial
    for(size_t k = d+1; k < num_check; k++){
        GR<K, D> xk;
        for(size_t i = 0; i < d+1; i++){
            auto L = GR<K, D>::one();
            for(size_t j = 0; j < d+1; j++){
                if(i==j)
                    continue;
                L = L * ((X(k) - X(j)) / (X(i) - X(j)));
            }
            xk = xk + L * shares[i];
        }


        if(!(xk == shares[k])){
#ifdef TESTING
        std::cout << "k = " << k << "\nxk \t\t= " << xk << "\nshares[" << k << "] \t= " << shares[k] << "\n";
#endif
            throw std::runtime_error("inconsistent reconstruction!");
        }
    }

    return secret;

#undef X
}

// interpolate
template<size_t K, size_t D>
vector<GR<K, D>> ShamirProtocol<K, D>::recon_coeffs(const vector<GR<K, D>> &shares, const size_t d, const size_t num_check){
    // num_check says how many shares to check lie on the polynomial
    //  (at least d+1, but not more than n)
    const size_t num_shares = shares.size();

    assert (d < num_shares);
    assert (d < n);
    assert (num_check <= num_shares);

    // first we're going to get the coefficients, then we'll check that at least num_check shares match the coefficients
    vector<GR<K,D>> coeffs;
    coeffs.resize(d+1);
    if(d==n-t){
        // if we're given something of this particular size, then
        // use the inverse of the vandermonde matrix that was precomputed
        for(size_t i = 0; i < d+1; i++){
            coeffs[i] = inv_vand_matrix[i][0] * shares[0];
            for(size_t j = 1; j < d+1; j++){
                coeffs[i] = coeffs[i] + inv_vand_matrix[i][j] * shares[j];
            }
        }
    }else{
        // if we're given anything else, then use normal gaussian elim to compute coeffs
        coeffs = solve_lin_eq(vand_matrix, shares, d+1);
    }

#define X(i) excep_seq[i]
    // Check that the other points are consistent with the polynomimal
    if(d+1 != n && num_check != d+1){
        // Checks that shares are consistent with polynomial
        for(size_t k = d+1; k < num_check; k++){
            GR<K, D> xk;
            for(size_t i = 0; i < d+1; i++){
                auto L = GR<K, D>::one();
                for(size_t j = 0; j < d+1; j++){
                    if(i==j)
                        continue;
                    L = L * ((X(k) - X(j)) / (X(i) - X(j)));
                }
                xk = xk + L * shares[i];
            }


            if(!(xk == shares[k])){
    #ifdef TESTING
            std::cout << "k = " << k << "\nxk \t\t= " << xk << "\nshares[" << k << "] \t= " << shares[k] << "\n";
    #endif
                throw std::runtime_error("inconsistent reconstruction!");
            }
        }
    }

    return coeffs;

#undef X
}

// this current assumes everyone inputs the same number of shares
template<size_t K, size_t D>
auto ShamirProtocol<K,D>::input(const vector<int>& inputs) -> vector<share_t>{

    const size_t num_inputs = inputs.size();
    const size_t share_size = GR<K, D>::byte_size();

    vector<vector<unsigned char>> sbufs(n);
    vector<vector<unsigned char>> rbufs(n);

    // prep the buffers to send/receive
    for(size_t i = 0; i < n; i++){
        sbufs[i].resize(num_inputs * share_size);
        rbufs[i].resize(num_inputs * share_size);
    }

    // create shares for each inputs
    vector<share_t> shares(n);
    for(size_t i = 0; i < num_inputs; i++){
        shares = share(inputs[i], t);
        for(size_t j = 0; j < n; j++){
            shares[j].pack(sbufs[j].data() + i*share_size);
        }
    }

    _network.exchange_all(sbufs, rbufs);

    vector<share_t> input_shares(n * num_inputs);

    for(size_t i = 0; i < (size_t)n; i++){
        for(size_t j = 0; j < num_inputs; j++){
            input_shares[i*num_inputs + j] = GR<K, D>{rbufs[i].data() + (j * share_size)};
        }
    }

    return input_shares;

    // optionally: MAC the inputs, but don't think we're doing that in this protocol
}

// does some rearranging, but calls recon_publ to open the values
template<size_t K, size_t D>
vector<GR<K, D>> ShamirProtocol<K,D>::open(const vector<GR<K, D>>& shares, const size_t d){
    const size_t batch_size = n-t;
    const size_t num_shares = shares.size();

    const size_t num_batches = num_shares/batch_size;

    vector<vector<GR<K,D>>> xs(num_batches);
    for(size_t i = 0; i < num_batches; i++){
        xs[i].resize(batch_size);
        for(size_t j = 0; j < batch_size; j++){
            xs[i][j] = shares[i*batch_size + j];
        }
    }

    vector<GR<K,D>> vals;
    vector<GR<K,D>> ss;

    for(size_t i = 0; i < num_batches; i++){
        ss = recon_publ(xs[i], d);
        vals.insert(vals.end(), ss.begin(), ss.end());
    }

    const size_t num_extras = num_shares % batch_size;
    // // If the number of shares is not divisble by n-t, we'll just do them separately
    if(num_extras != 0){
        vector<GR<K,D>> extras(batch_size, GR<K,D>::zero());
        for(size_t i = 0; i < num_extras; i++){
            extras[i] = shares[num_batches * batch_size + i];
        }

        ss = recon_publ(extras, d);
        vals.insert(vals.end(), ss.begin(), ss.begin()+num_extras);
    }

    return vals;
}

template<size_t K, size_t D>
auto ShamirProtocol<K,D>::reduce_degree_regular(const vector<share_2t>& xs) -> vector<share_t>{
    // 1. subtract by the 2t double share
    // 2. open z' = z - 2t share of r
    // 3. parties take share as z' + t share of r

    const size_t num_shares = xs.size();

    if((num_shares + _dcntr) > _dshares.size()){
        throw std::runtime_error("Not enough preprocessed double shares for degree reduction!");
    }

    vector<share_2t> z_2t(num_shares);
    for(size_t i = 0; i < num_shares; i++){
        z_2t[i] = xs[i] - _dshares[_dcntr + i][1];
    }

    auto z = open(z_2t, 2*t);

    vector<share_t> shares(num_shares);
    for(size_t i = 0; i < num_shares; i++){
        shares[i] = z[i] + _dshares[_dcntr + i][0];
    }

    _dcntr += num_shares;

    return shares;
}

template<size_t K, size_t D>
auto ShamirProtocol<K,D>::reduce_degree_simd(const vector<share_2t>& xs) -> vector<share_t>{

    const size_t num_shares = xs.size();

    if((num_shares + _sdcntr) > _simddshares.size()){
        throw std::runtime_error("Not enough preprocessed double shares for degree reduction!");
    }

    vector<share_2t> z_2t(num_shares);
    for(size_t i = 0; i < num_shares; i++){
        z_2t[i] = xs[i] - _simddshares[_sdcntr + i][1];
    }

    auto z_out = open(z_2t, 2*t);
    auto z_in =  SIMDEout2Ein<K,D>(z_out);

    vector<share_t> shares(num_shares);
    for(size_t i = 0; i < num_shares; i++){
        shares[i] = z_in[i] + _simddshares[_sdcntr + i][0];
    }

    _sdcntr += num_shares;

    return shares;
}

template<size_t K, size_t D>
auto ShamirProtocol<K,D>::reduce_degree_flex(const vector<share_2t>& xs) -> vector<share_t>{

    const size_t num_shares = xs.size();
    const size_t batch_size = D;

    if((num_shares/batch_size + _fdcntr) > _flexdshares.size()){
        throw std::runtime_error("Not enough preprocessed double shares for degree reduction!");
    }

    // only handle full batches for now
    assert(num_shares % batch_size == 0);
    const size_t num_batches = num_shares / batch_size;

    vector<share_2t> ws(num_batches);
    for(size_t k = 0; k < num_batches; k++){
        // <z> = sum <xs_i> \cdot \xi^i
        share_2t z = xs[k*batch_size];
        size_t powi = 1;
        for(size_t i = 1; i < D; i++){
            powi *= 2;
            z += xs[k*batch_size + i] * excep_seq[powi - 1];
        }

        // <w>_2t = <z>_2t - <r>_2t
        ws[k] = z - _flexdshares[_fdcntr + k][D];
    }

    vector<GR<K,D>> w = open(ws, 2*t);

#ifdef TESTING
    std::cout << "w = " << w[0] << "\n";
#endif

    vector<share_t> shares(num_shares);
    for(size_t i = 0; i < num_batches; i++){
        for(size_t j = 0; j < batch_size; j++){
            shares[i*batch_size + j] = w[i][j] + _flexdshares[_fdcntr + i][j];
        }
    }

    _fdcntr+=num_batches;

    return shares;
}

template<size_t K, size_t D>
auto ShamirProtocol<K,D>::reduce_degree_innerprod(const vector<share_2t>& xs) -> vector<share_t>{
    // todo
    std::cout << "IMPLEMENT ME\n";
    return xs;
}

template<size_t K, size_t D>
void ShamirProtocol<K,D>::broadcast(const unsigned int sender_id, vector<GR<K,D>>& xs){
    // Non robust broadcast:
    //  1. Pi sends x1,...,x(n-t) to every Pj
    //  2. Every Pj applies HIM to x to get vector y1,...,yn
    //  3. Every Pj sends yk to Pk
    //  4. Pk checks all received yk are equal. Send 1 if so, otherwise send 0
    //  5. Abort if you receive 0 from anyone

    // note this is hard-coded to broadcast n-t values
    const size_t share_size = GR<K, D>::byte_size();
    const size_t num_vals = n-t;

    // Pi sends x1,...,x(n-t) to every other party
    if(sender_id==_network.id()){
        assert(xs.size() == num_vals);

        vector<unsigned char> bcast_sbuf(num_vals * share_size);

        for(size_t i = 0; i < num_vals; i++){
            xs[i].pack(bcast_sbuf.data() + i * share_size);
        }

        _network.broadcast_send(bcast_sbuf);

    } else {
        vector<unsigned char> bcast_rbufs(share_size*num_vals);

        _network.broadcast_recv(sender_id, bcast_rbufs);
        xs.resize(num_vals);
        for(size_t i = 0; i < num_vals; i++){
            xs[i] = GR<K,D>{bcast_rbufs.data() + (i * share_size)};
        }
    }

    // compute y1,...,yn and send yk to Pk
    vector<GR<K,D>> ys = matrix_mult(hyper_matrix, xs, n);
    vector<vector<unsigned char>> y_sbufs(n);
    vector<vector<unsigned char>> y_rbufs(n);
    for(size_t i = 0; i < n; i++){
        y_rbufs[i].resize(share_size);
        y_sbufs[i].resize(share_size);
        ys[i].pack(y_sbufs[i].data());
    }

    _network.exchange_all(y_sbufs, y_rbufs);

    bool is_ok = true;

    // check that all received yk are the same
    GR<K,D> yk = GR<K,D>{y_rbufs[0].data()};
    for(size_t i = 1; i < n; i++){
        is_ok = is_ok && (yk == GR<K,D>{y_rbufs[i].data()});

#ifdef TESTING
        if(yk != GR<K,D>{y_rbufs[i].data()}){
            std::cout << "Inconsistent on " << i << "\n";
        }
#endif
    }

    // send the is_ok bit to everyone
    vector<unsigned char> bit_sbuf(sizeof(bool));
    vector<vector<unsigned char>> bit_rbufs(n);

    memcpy(bit_sbuf.data(), &is_ok, sizeof(bool));

#ifdef TESTING
    std::cout << "sending bit " << is_ok << "\n";
#endif

    _network.broadcast_send(bit_sbuf);
    for(size_t i = 0; i < n; i++){
        bit_rbufs[i].resize(sizeof(bool));
        _network.broadcast_recv(i, bit_rbufs[i]);

        bool r_bit;
        memcpy(&r_bit, bit_rbufs[i].data(), sizeof(bool));

#ifdef TESTING
        std::cout << "bit received from party "<< i <<": " << r_bit << "\n";
#endif

        is_ok = is_ok && r_bit;
    }

    // abort if
    if(!is_ok){
        throw std::runtime_error("inconsistent broadcast");
    }
}

template<size_t K, size_t D>
void ShamirProtocol<K,D>::recon_priv(GR<K,D> share, const size_t d, const unsigned int recv_id, GR<K,D>& secret){
    const size_t share_size = GR<K, D>::byte_size();
    vector<unsigned char> sbufs(share_size);

    share.pack(sbufs.data());

    _network.send_to(recv_id, sbufs);

    // only reconstruct if I'm the person whose supposed to get the output
    if(recv_id == _network.id()){
        vector<vector<unsigned char>> rbufs(n);
        vector<share_t> shares(n);

        for(size_t i = 0; i < (size_t)n; i++){
            rbufs[i].resize(share_size);
            _network.recv_from(i, rbufs[i]);
            shares[i] = GR<K, D>{rbufs[i].data()};
        }

        // Reconstruct and check that t+d+1 points are consistent with polynomial.
        secret = reconstruct(shares, d, t+d+1);
    }
}

template<size_t K, size_t D>
vector<GR<K,D>> ShamirProtocol<K,D>::recon_publ(vector<GR<K,D>> &shares, const size_t d){

    // This is written for exactly this many shares. Dunno what to do if we receive some other number
    assert(shares.size() == (n-t));

    vector<share_t> u_shares(n);
    for(size_t i = 0; i < n; i++){
        GR<K,D> alpha = GR<K,D>::one();
        u_shares[i] = shares[0];
        for(size_t j = 1; j < n-t; j++){
            alpha = alpha * excep_seq[i];
            u_shares[i] = u_shares[i] + shares[j]*alpha;
        }
    }

    // reconstruct [u_i] to P_i
    vector<GR<K,D>> us(n); // what's a better way to do output from recon_priv?
    for(size_t i = 0; i < n; i++){
        recon_priv(u_shares[i], d, i, us[i]);
    }

    // send u_i to every other party
    const size_t share_size = GR<K, D>::byte_size();
    // only going to send one value to everyone (u_i)
    vector<unsigned char> bcast_sbuf(share_size);
    // going to receive one value from everyone
    vector<vector<unsigned char>> bcast_rbufs(n);

    for(size_t i = 0; i < n; i++){
        bcast_rbufs[i].resize(share_size);
    }

    us[_network.id()].pack(bcast_sbuf.data());
    _network.broadcast_send(bcast_sbuf); // send this single value to everyone

    for(size_t i = 0; i < n; i++){
        _network.broadcast_recv(i, bcast_rbufs[i]);
        us[i] = GR<K,D>{bcast_rbufs[i].data()};
    }

    // reconstruct and check that n-t points lie on the polynomial
    // use recon_coeffs to check n-t points lie on the poly and get the coefficients
    vector<GR<K,D>> coeffs = recon_coeffs(us, n-t-1, n-t);

    return coeffs;
}

#ifdef TESTING
template<size_t K, size_t D>
void ShamirProtocol<K,D>::test(){

    std::cout << "n = " << n << "\nt = " << t << "\n";

    vector<int> v = { (int)(_network.id()) };
    vector<share_t> shares = input(v);

    vector<share_t> shares_to_recon(n-t);
    for(size_t i = 0; i < n-t; i++){
        shares_to_recon[i] = shares[i];
    }

    auto r = recon_publ(shares_to_recon, t);

    bool is_good = true;

    for(size_t i = 0; i < n-t; i++){
        is_good = is_good && r[i] == GR<K,D>{(int)i};
    }

    if(!is_good){
        std::cout << "There was an issue in reconstructing the coefficients for recon_publ!\n";
    } else {
        std::cout << "recon_publ was fine!\n";
    }
}

template<size_t K, size_t D>
void ShamirProtocol<K,D>::test_bcast(){

    std::cout << "n = " << n << "\nt = " << t << "\n";

    vector<share_t> x(n-t);
    if(_network.id()==0){
        for(size_t i = 0; i < n-t; i++){
            x[i] = GR<K,D>{ 5 };
        }
    }

    broadcast(0, x);
}

template<size_t K, size_t D>
void ShamirProtocol<K,D>::test_flex_lincombo(){
    vector<int> v = { (int)(_network.id()) };
    vector<share_t> yis = input(v);

    // 0
    share_2t z = yis[0] * excep_seq[0];
    // 1
    z += yis[1] * excep_seq[1];
    // 2
    z += yis[2] * excep_seq[3];
    // 3
    z += yis[3] * excep_seq[7];

    vector<share_2t> zs(1);
    zs[0] = z;

    vector<GR<K,D>> y = open(zs, t);
    std::cout << "y = " << y[0] << "\n";

    auto opened = open(yis, t);

    for(size_t i = 0; i < 4; i++){
        std::cout << "y = " << y[0][i] << "\n";
        std::cout << "yis " << opened[i][0] << "\n";
    }
}

template<size_t K, size_t D>
void ShamirProtocol<K,D>::test_flex_shares(){
    const size_t batch_size = n-t;
    const size_t num_batches = (50000 / D) / batch_size;
    std::cout << "Testing flex_deshare num_batchs=" << num_batches << " total=" << (num_batches*batch_size) << "\n";

    preproc_dshares_flex(num_batches);

}

template<size_t K, size_t D>
auto ShamirProtocol<K,D>::test_make_share() -> vector<share_t> {
    // generate a sharing of v

    vector<share_t> points (t+1);
    vector<share_t> shares (n);
    points[0] = GR<K, D>{ (int)_network.id() };

    // This is polynomial evaluation via. Horners method.
    for (size_t i = 1; i < t+1; i++)
        points[i] = randomize<share_t>();

    for (size_t i = 0; i < n; i++) {
        const auto x = excep_seq[i+n];
        auto px = points[t] * x + points[t-1];

        // note: can't use size_t because size_t doesn't have negative values
        for (ssize_t j = t-2; j >= 0; j--){
            std::cout << "j = " << j << "\n";
            px = px * x + points[j];
        }

        shares[i] = px;
        std::cout << "shares[" << i << "] = " << shares[i] << "\n";
    }

    return shares;
}

template<size_t K, size_t D>
void ShamirProtocol<K,D>::test_deg_reduce(){

    const size_t batch_size = n-t;
    const size_t num_shares = batch_size * D * 2;
    preproc_dshares_regular(num_shares);

    std::cout << "Testing degree reduction now\n";
    Z2k<K> val_to_share0 = randomize<Z2k<K>>();
    Z2k<K> val_to_share1 = randomize<Z2k<K>>();
    
    vector<share_2t> shares_2t = share(val_to_share0, 2*t);
    vector<share_2t> shares2_2t = share(val_to_share1, 2*t);

    const size_t share_size = GR<K,D>::byte_size();
    vector<vector<unsigned char>> sbufs(n);
    vector<vector<unsigned char>> rbufs(n);
    for(size_t i = 0; i < n; i++){
        sbufs[i].resize(2*share_size);
        rbufs[i].resize(2*share_size);

        shares_2t[i].pack(sbufs[i].data());
        shares2_2t[i].pack(sbufs[i].data() + share_size);
    }

    _network.exchange_all(sbufs, rbufs);

    vector<share_2t> received_vals(2*n);
    for(size_t i = 0; i < n; i++){
        received_vals[2*i] = GR<K,D>{ rbufs[i].data() };
        received_vals[2*i + 1] = GR<K,D>{ rbufs[i].data() + share_size };
    }

    // we're only going to reduce the first D
    received_vals.erase(received_vals.begin() + 2*D, received_vals.end());
    assert(received_vals.size() == 2*D);

    std::cout << "Calling reduce_degree_regular\n";
    vector<share_t> new_shares = reduce_degree_regular(received_vals);

    std::cout << "Opening values\n";
    vector<GR<K,D>> opened_vals = open(new_shares, t);
    for(size_t i = 0; i < opened_vals.size(); i++){
        if((_network.id() * 2)==i){
            if(opened_vals[i]==val_to_share0)
            {
                std::cout << "The shared value matches val_to_share0!\n";
            } else{
                std::cout << "Hmmm" << val_to_share0 << "\n";
            }
        } else if((_network.id() * 2 + 1) == i){
            if(opened_vals[i]==val_to_share1){
                std::cout << "The shared value matches val_to_share1!\n";
            } else{
                std::cout << "Hmmm" << val_to_share1 << "\n";
            }
        }
    }
}

template<size_t K, size_t D>
void ShamirProtocol<K,D>::test_simd(){

    const size_t batch_size = n-t;
    const size_t num_shares = batch_size * D * 100;
    preproc_dshares_simd(num_shares);

    vector<share_t> r_t(num_shares);
    vector<share_2t> r_2t(num_shares);

    for(size_t i = 0; i < num_shares; i++){
        r_t[i] = _simddshares[i][0];
        r_2t[i] = _simddshares[i][1];
    }

    auto xs_t = open(r_t, t);
    auto xs_2t = open(r_2t, 2*t);

    for(size_t i = 0; i < num_shares; i++){
        if(!SIMDencodings_match(xs_t[i], xs_2t[i])){
            std::cout << "mismatch at i = " << i << "!\n";
            std::cout << "xs_t = " << xs_t[i] << "\n";
            std::cout << "xs_2t = " << xs_2t[i] << "\n";
        }
    }

    std::cout << "Things were fine?\n";

    std::cout << "Testing degree reduction now\n";
    vector<Z2k<K>> val_to_share0(2);
    vector<Z2k<K>> val_to_share1(2);
    for(size_t i = 0; i < 2; i++){
        val_to_share0[i] = randomize<Z2k<K>>();
        val_to_share1[i] = randomize<Z2k<K>>();
    }

    GR<K,D> eout0 = SIMDEout<K,D>(val_to_share0);
    GR<K,D> eout1 = SIMDEout<K,D>(val_to_share1);
    vector<share_2t> shares_2t = share(eout0, 2*t);
    vector<share_2t> shares2_2t = share(eout1, 2*t);

    const size_t share_size = GR<K,D>::byte_size();
    vector<vector<unsigned char>> sbufs(n);
    vector<vector<unsigned char>> rbufs(n);
    for(size_t i = 0; i < n; i++){
        sbufs[i].resize(2*share_size);
        rbufs[i].resize(2*share_size);

        shares_2t[i].pack(sbufs[i].data());
        shares2_2t[i].pack(sbufs[i].data() + share_size);
    }

    _network.exchange_all(sbufs, rbufs);

    vector<share_2t> received_vals(2*n);
    for(size_t i = 0; i < n; i++){
        received_vals[2*i] = GR<K,D>{ rbufs[i].data() };
        received_vals[2*i + 1] = GR<K,D>{ rbufs[i].data() + share_size };
    }

    // we're only going to reduce the first D
    received_vals.erase(received_vals.begin() + 2*D, received_vals.end());
    assert(received_vals.size() == 2*D);

    std::cout << "Calling reduce_degree_simd\n";
    vector<share_t> new_shares = reduce_degree_simd(received_vals);

    std::cout << "Opening values\n";
    vector<GR<K,D>> opened_vals = open(new_shares, t);
    for(size_t i = 0; i < opened_vals.size(); i++){
        if((_network.id() * 2)==i){
            if(SIMDencodings_match<K,D>(opened_vals[i], eout0))
            {
                std::cout << "The shared value matches eout0!\n";
            } else{
                std::cout << "Hmmm" << eout0 << "\n";
            }
        } else if((_network.id() * 2 + 1) == i){
            if(SIMDencodings_match<K,D>(opened_vals[i], eout1)){
                std::cout << "The shared value matches eout1!\n";
            } else{
                std::cout << "Hmmm" << eout1 << "\n";
            }
        }
    }
}

template<size_t K, size_t D>
void ShamirProtocol<K,D>::test_prg(){
    const size_t num_bytes = 1;
    unsigned char buf[num_bytes];
    _prg->next(buf, num_bytes);
    std::cout << (1 & buf[0]) << "\n";
    std::cout << (1 & buf[0]>>1) << "\n";
    std::cout << (1 & buf[0]>>2) << "\n";
    std::cout << (1 & buf[0]>>3) << "\n";
    std::cout << (1 & buf[0]>>4) << "\n";
    std::cout << (1 & buf[0]>>5) << "\n";
    std::cout << (1 & buf[0]>>6) << "\n";
    std::cout << (1 & buf[0]>>7) << "\n";
}

#endif

} // goat


#endif // _SHAMIR_HPP
