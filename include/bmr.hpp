#ifndef _GOAT_BMR_HPP
#define _GOAT_BMR_HPP

#include "ncomm.hpp"

namespace goat {

// TODO
class Circuit;
class GarbledCircuit;

using ncomm::Network;

class BMR {
public:

    BMR(Network& network) :
	_network{network}, n{network.size()}, t{(n-1)/3} {};

    GarbledCircuit garble(Circuit& circuit);

    // TODO
    void run(GarbledCircuit& circuit);

private:
    Network _network;
    size_t n;
    size_t t;
};

}

#endif // _GOAT_BMR_HPP
