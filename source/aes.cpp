// Below is taken from https://github.com/sebastien-riou/aes-brute-force
#include <string.h>
#include <wmmintrin.h>
#include <iostream>

#include "crypto.hpp"

namespace goat {

#define DO_ENC_BLOCK(m,k)				\
	do{						\
		m = _mm_xor_si128       (m, k[ 0]);	\
		m = _mm_aesenc_si128    (m, k[ 1]);	\
		m = _mm_aesenc_si128    (m, k[ 2]);	\
		m = _mm_aesenc_si128    (m, k[ 3]);	\
		m = _mm_aesenc_si128    (m, k[ 4]);	\
		m = _mm_aesenc_si128    (m, k[ 5]);	\
		m = _mm_aesenc_si128    (m, k[ 6]);	\
		m = _mm_aesenc_si128    (m, k[ 7]);	\
		m = _mm_aesenc_si128    (m, k[ 8]);	\
		m = _mm_aesenc_si128    (m, k[ 9]);	\
		m = _mm_aesenclast_si128(m, k[10]);	\
	}while(0)

#define AES_128_key_exp(k, rcon) aes_128_key_expansion(k, _mm_aeskeygenassist_si128(k, rcon))

inline static __m128i aes_128_key_expansion(__m128i key,
					    __m128i keygened)
{
	keygened = _mm_shuffle_epi32(keygened, _MM_SHUFFLE(3,3,3,3));
	key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
	key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
	key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
	return _mm_xor_si128(key, keygened);
}

inline static void aes128_load_key(unsigned char *enc_key,
				   __m128i *key_schedule)
{
	key_schedule[0] = _mm_loadu_si128((const __m128i*) enc_key);
	key_schedule[1]  = AES_128_key_exp(key_schedule[0], 0x01);
	key_schedule[2]  = AES_128_key_exp(key_schedule[1], 0x02);
	key_schedule[3]  = AES_128_key_exp(key_schedule[2], 0x04);
	key_schedule[4]  = AES_128_key_exp(key_schedule[3], 0x08);
	key_schedule[5]  = AES_128_key_exp(key_schedule[4], 0x10);
	key_schedule[6]  = AES_128_key_exp(key_schedule[5], 0x20);
	key_schedule[7]  = AES_128_key_exp(key_schedule[6], 0x40);
	key_schedule[8]  = AES_128_key_exp(key_schedule[7], 0x80);
	key_schedule[9]  = AES_128_key_exp(key_schedule[8], 0x1B);
	key_schedule[10] = AES_128_key_exp(key_schedule[9], 0x36);
}

inline static void aes128_enc(__m128i *key_schedule,
			      unsigned char *plainText,
			      unsigned char *cipherText)
{
	__m128i m = _mm_loadu_si128((__m128i *) plainText);
	DO_ENC_BLOCK(m, key_schedule);
	_mm_storeu_si128((__m128i *) cipherText, m);
}

void aes_prg(unsigned char* buf, const unsigned char seed[16], size_t size) {

    size_t nblocks = (size/16 + 1);
    unsigned char *out = (unsigned char *)malloc(nblocks * 16);

    if (!out)
	throw -1;

    __m128i ks[11];
    aes128_load_key((unsigned char *)seed, ks);
    unsigned char *p = out;
    unsigned long counter[4] = {0};

    for (size_t i = 0; i < nblocks; i++) {
	aes128_enc(ks, (unsigned char *)counter, p);
	p += 16;
	counter[0]++;
    }

    memcpy(buf, out, size);
    free(out);
}

void aes_prg(unsigned char* buf, size_t size) {
    // obviously this key should be picked at random at some point during
    // init. This suffices for testing though.
    static const unsigned char k[16] = {
	0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
    };
    aes_prg(buf, k, size);
}

}
