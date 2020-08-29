/*
 * Created by Anders Dalskov
 * Date: 31-10-2018
 *
 * Lots of code copied from SPDZ-2 library.
 *
 * word == unsigned long
 */

#include "GF2_128.hpp"

namespace goat {

int GF2_128::n;
int GF2_128::t1;
int GF2_128::t2;
int GF2_128::t3;
int GF2_128::l0;
int GF2_128::l1;
int GF2_128::l2;
int GF2_128::l3;
int GF2_128::nterms;
int128 GF2_128::mask;
int128 GF2_128::lowermask;
int128 GF2_128::uppermask;
bool GF2_128::rewind = false;

std::ostream& operator<<(std::ostream &s, const int128 &a) {
    unsigned long *tmp = (unsigned long*) & a.a;
    s << std::hex;
    s.width(16);
    s.fill('0');
    s << tmp[1];
    s.width(16);
    s << tmp[0] << std::dec;
    return s;
}

inline int128 int128::operator<<(const int& other) const {
    int128 res(_mm_slli_epi64(a, other));
    __m128i mask;
    if (other < 64)
	mask = _mm_srli_epi64(a, 64 - other);
    else
	mask = _mm_slli_epi64(a, other - 64);
    res.a ^= _mm_slli_si128(mask, 8);
    return res;
}

inline int128 int128::operator>>(const int& other) const {
    int128 res(_mm_srli_epi64(a, other));
    __m128i mask;
    if (other < 64)
	mask = _mm_slli_epi64(a, 64 - other);
    else
	mask = _mm_srli_epi64(a, other - 64);
    res.a ^= _mm_srli_si128(mask, 8);
    return res;
}

void GF2_128::initField() {

    nterms = 3;

    // coefficients of irred f.
    n = 128; t1 = 7;  t2 = 2;  t3 = 1;
    l0 = 0;  l1 = t1; l2 = t2; l3 = t3;

    mask = _mm_set_epi64x(-1, -1);
    lowermask = _mm_set_epi64x((1LL << (64 - 7)) - 1, -1);
    uppermask = _mm_set_epi64x(((unsigned long) - 1) << (64 - 7), 0);
}

static inline void avx_memzero(void* dest, size_t length) {
    __m256i* d = (__m256i*)dest;

    __m256i s = _mm256_setzero_si256();
    while (length >= 32) {
	_mm256_storeu_si256(d++, s);
	length -= 32;
    }

    if (length)
	memset((void*)d, 0, length);
}

GF2_128 &GF2_128::mul(const GF2_128 &x, const GF2_128 &y) {
    __m128i res[2];
    avx_memzero(res, sizeof(res));
    mul128(x.a.a, y.a.a, res, res+1);
    // REDUCE_PENT((int128)res[1], (int128)res[0]);
    reducePentanomial(res[1], res[0]);
    return *this;
}

void GF2_128::reduceTrinomial(int128 xh, int128 xl) {

    a = xl;
    a ^= (xh << l0);
    a ^= (xh << l1);

    int128 hi = a >> n;
    while (hi == 0) {
	a &= mask;
	a ^= hi;
	a ^= (hi << t1);
	hi = a >> n;
    }
}

inline void GF2_128::reducePentanomial(int128 xh, int128 xl) {

    int128 upper, lower;
    int128 tmp = 0;

    a = xl;
    upper = xh & uppermask;
    lower = xh & lowermask;

    // Upper part
    tmp ^= (upper >> (n - t1 - l0));
    tmp ^= (upper >> (n - t1 - l1));
    tmp ^= (upper >> (n - t1 - l2));
    tmp ^= (upper >> (n - t1 - l3));
    lower ^= (tmp >> (l1));
    a ^= (tmp << (n - l1));

    // Lower part
    a ^= (lower << l0);
    a ^= (lower << l1);
    a ^= (lower << l2);
    a ^= (lower << l3);
}

bool is_ge(__m128i a, __m128i b) {
    unsigned long aa[2];
    unsigned long bb[2];
    _mm_storeu_si128((__m128i*)aa, a);
    _mm_storeu_si128((__m128i*)bb, b);
    //  cout << hex << "is_ge " << aa[1] << " " << bb[1] << " " << (aa[1] > bb[1]) << " ";
    //  cout << aa[0] << " " << bb[0] << " " << (aa[0] >= bb[0]) << endl;
    return aa[1] == bb[1] ? aa[0] >= bb[0] : aa[1] > bb[1];
}

class int129 {
    private:

    int128 lower;
    bool msb;

    public:

    int129(): lower(_mm_setzero_si128()), msb(false) {}
    int129(int128 lower, bool msb) : lower(lower), msb(msb) {}
    int129(int128 a) : lower(a), msb(false) {}
    int129(unsigned long a) {
	*this = a;
    }

    int128 get_lower() {
	return lower;
    }

    int129& operator=(const __m128i& other) {
	lower = other;
	msb = false;
	return *this;
    }

    int129& operator=(const unsigned long& other) {
	lower = _mm_set_epi64x(0, other);
	msb = false;
	return *this;
    }

    bool operator==(const int129& other) {
	return (lower == other.lower) && (msb == other.msb);
    }

    bool operator!=(const int129& other) {
	return !(*this == other);
    }

    bool operator>=(const int129& other) {
	//cout << ">= " << msb << other.msb <<  (msb > other.msb) << is_ge(lower.a, other.lower.a) << endl;
        return msb == other.msb ? is_ge(lower.a, other.lower.a) : msb > other.msb;
    }

    int129 operator<<(int other) {
	return int129(lower << other, _mm_cvtsi128_si32(((lower >> (128-other)) & 1).a));
    }

    int129& operator>>=(int other) {
	lower >>= other;
	lower |= (int128(msb) << (128-other));
	msb = !other;
	return *this;
    }

    int129 operator^(const int129& other) {
	return int129(lower ^ other.lower, msb ^ other.msb);
    }

    int129& operator^=(const int129& other) {
	lower ^= other.lower;
	msb ^= other.msb;
	return *this;
    }

    int129 operator&(const unsigned long& other) {
	return int129(lower & other, false);
    }

    friend std::ostream& operator<<(std::ostream& s, const int129& a) {
	s << a.msb << a.lower;
	return s;
    }
};

void GF2_128::invert() {

    if (is_one())
	return;

    if (is_zero())
	throw -1;  // division by zero

    int129 u;
    int129 v = a;
    int129 B = 0;
    int129 D = 1;
    int129 mod = 1;

    mod ^= (int129(1) << n);
    mod ^= (int129(1) << t1);

    if (nterms == 3) {
	mod ^= (int129(1) << t2);
	mod ^= (int129(1) << t3);
    }

    u = mod;
    v = a;

    while (u!=0) {

	while ((u&1)==0) {
	    u >>= 1;
	    if ((B & 1)!=0)
		B ^= mod;
	    B >>= 1;
	}

	while ((v & 1) == 0 && v != 0) {
	    v >>= 1;
	    if ((D & 1) != 0)
		D ^= mod;
	    D >>= 1;
	}

	if (u >= v) {
	    u = u ^ v;
	    B = B ^ D;
	} else {
	    v = v ^ u;
	    D = D ^ B;
	}
    }

    a = D.get_lower();
}

}
