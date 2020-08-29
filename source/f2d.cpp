#include "f2d.hpp"

namespace goat {

template<>
elem_t F2d<4>::mtable[16][16] = {0};

template<>
elem_t F2d<4>::itable[16] = {0};

template<>
bool F2d<4>::init_done = false;

template<>
void F2d<4>::init_tables() {
    if (init_done)
	return;

    size_t ic = 0;
    for (size_t i = 0; i < mask + 1; i++) {
	for (size_t j = 0; j < mask + 1; j++) {
	    // hardcode H(X) = X^4 + X + 1
	    const auto c = mul<0x08, 0b10011>(i, j);
	    if (c == 1) {
		itable[i] = j;
		ic++;
	    }
	    mtable[i][j] = c;
	}
    }

    if (ic != mask)
	throw std::runtime_error("error occured during init of f2d");

    init_done = true;
}

template<>
F2d<4>::F2d(const elem_t e) {
    init_tables();
    this->e = e & mask;
}

}
