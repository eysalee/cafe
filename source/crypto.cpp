#include "crypto.hpp"

#include <stdexcept>
#include <cassert>
#include <iostream>

namespace goat {

void prg(unsigned char *buf, size_t size) {
#ifdef AES_PRG
    aes_prg(buf, size);
#else
    randombytes_buf(buf, size);
#endif
}

void prg(unsigned char *buf, size_t size, unsigned char *seed) {
#ifdef AES_PRG
    aes_prg(buf, seed, size);
#else
    randombytes_buf_deterministic(buf, size, seed);
#endif
}


void crypto_init() {
    if (sodium_init() < 0)
	throw std::runtime_error("sodium_init()");
}

template<>
uint64_t randomize() {
    crypto_init();

    uint64_t r;
    prg((unsigned char *)&r, sizeof(r));
    return r;
}

StatefulPRG::StatefulPRG(unsigned char seed[PRG_KEY_SIZE]) {
    memcpy(_seed, seed, PRG_KEY_SIZE);

    std::cout << "buffering prg...\n";

    _spent = 0;
    prg(_buffer, PRG_BUFFER_SIZE, _seed);
    _update_state();
}

void StatefulPRG::next(unsigned char *buf, const size_t length) {

    if (length > PRG_BUFFER_SIZE) {
	prg(buf, length, _seed);
	goto done;
    }

    if (length > (size_t)(PRG_BUFFER_SIZE - _spent)) {
	prg(_buffer, PRG_BUFFER_SIZE, _seed);
	_spent = 0;
    }

    memcpy(buf, _buffer + _spent, length);
    _spent += length;

done:
    _update_state();
}

void StatefulPRG::_update_state() {
    unsigned long long *s = (unsigned long long *)_seed;
    *s += 1;
}

}
