#ifndef _GOAT_CRYPTO_HPP
#define _GOAT_CRYPTO_HPP

#include <vector>
#include <cstring>

#include <sodium.h>

#ifdef AES_PRG
#define PRG_KEY_SIZE  16
#else
#define PRG_KEY_SIZE  randombytes_SEEDBYTES
#endif

#define PRG_BUFFER_SIZE 1024

namespace goat {

void aes_prg(unsigned char *buf, const unsigned char seed[PRG_KEY_SIZE], size_t len);
void aes_prg(unsigned char *buf, size_t len);

void prg(unsigned char *buf, size_t size, unsigned char *seed);
void prg(unsigned char *buf, size_t size);

void crypto_init();

class StatefulPRG {
public:
    StatefulPRG(unsigned char seed[PRG_KEY_SIZE]);
    void next(unsigned char *buf, const size_t length);

private:
    unsigned char _buffer[PRG_BUFFER_SIZE];
    size_t _spent;

    unsigned char _seed[PRG_KEY_SIZE];
    void _update_state();
};

template<typename V>
V randomize() {

    crypto_init();

    unsigned char buf[V::byte_size()];
    prg(buf, V::byte_size());
    // randombytes_buf(buf, V::byte_size());
    return V{buf};
}

template<typename V>
std::vector<V> randomize(const size_t n) {
    crypto_init();

    std::vector<V> r (n);
    if (std::is_fundamental<V>::value) {
	prg(r.data(), n*sizeof(V));
	// randombytes_buf(r.data(), n * sizeof(V));
    } else {
	for (size_t i = 0; i < n; i++)
	    r[i] = randomize<V>();
    }
    return r;
}

template<> uint64_t randomize();

} // goat

#endif // _GOAT_CRYPTO_HPP
