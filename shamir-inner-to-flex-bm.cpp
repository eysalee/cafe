#include "shamir.hpp"
#include "ncomm.hpp"
#include "timing.hpp"

#include <iostream>

using namespace ncomm;
using namespace std;
using namespace goat;

int main(int argc, char** argv) {

    if (argc < 3) {
	cout << "....\n";
	return -1;
    }

    partyid_t id = stoul(argv[1]);
    Network nw (id, argv[2]);
    nw.connect();

    const size_t K = 64;
    const size_t D = 4;

    auto prot = ShamirProtocol<K, D>(nw);

    auto num_shares = 500000;

    prot.preproc_dshares_innerprod(num_shares);

    cout << "Created " << prot.get_number_of_innerprod_dshares() << " shares\n";

}
