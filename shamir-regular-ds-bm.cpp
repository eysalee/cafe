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

    auto num_shares = stoul(argv[3]);
    // auto num_batches = num_shares / (prot.get_n() - 2*prot.get_t());

    prot.preproc_dshares_regular(num_shares);

    cout << "created " << prot.get_number_of_dshares() << " shares\n";
}
