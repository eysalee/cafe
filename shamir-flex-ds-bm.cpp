#include "shamir.hpp"

#include "ncomm.hpp"
#include "crypto.hpp"
#include "timing.hpp"

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

    // we can open D mults with one double share.
    auto num_shares = stoul(argv[3]);

    cerr << "num_shares=" << num_shares << "\n";

    // auto batch_size = prot.get_n() - prot.get_t();
    // auto num_batches = num_shares / batch_size;

    // prot.preproc_dshares_flex(num_batches);

    // code was updated to take in total number of shares you want to generate
    prot.preproc_dshares_flex(num_shares);

    cout << "Created " << prot.get_number_of_flex_dshares() << " shares\n";
}
