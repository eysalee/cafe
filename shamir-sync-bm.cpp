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

    std::cout << "Running make-share once\n";
    prot.test_make_share();

    std::cout << "Running make-share twice\n";
    prot.test_make_share();

}
