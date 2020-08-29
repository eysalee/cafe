#include "shamir.hpp"

#include "ncomm.hpp"
#include "crypto.hpp"

#include <iostream>

using namespace ncomm;
using namespace std;
using namespace goat;

int main(int argc, char** argv) {

    if (argc < 3) {
    cout << "usage: " << argv[0] << " <party ID> <network_info.txt>\n";
    return -1;
    }

    partyid_t id = stoul(argv[1]);

    Network nw (id, argv[2]);
    nw.connect();

    const size_t K = 64;
    const size_t D = 4;

    auto prot = ShamirProtocol<K, D>(nw);

#ifdef TESTING
    // testing functions go under here
    // prot.preproc_dshares_regular(2);
    // prot.test();
    // prot.test_bcast();
    // prot.test_flex_lincombo();
    // prot.test_flex_shares();
    // prot.test_prg();
    // // prot.sync();

    // size_t n = prot.get_n();
    // size_t t = prot.get_t();
    // size_t batch_size = D * (n-t);

    // prot.preproc_dshares_simd(batch_size * 2);

    prot.test_simd();

    // prot.gen_simd_zero_shares(batch_size * 2);

    // prot.preproc_dshares_regular(batch_size * 2);
    // prot.test_deg_reduce();

    // SIMDin<K, D>();
    // vector<Z2k<K>> xs(5, Z2k<K>{0});
    // SIMDout<K,D>(xs);

    // std::cout << "This is a test\n";
#endif
}