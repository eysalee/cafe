/*
 * Created by Anders Dalskov
 * Date: 31-10-2018
 */

#ifndef _HYPERBMR_GF2_128_H
#define _HYPERBMR_GF2_128_H

/*
 * TODO: Add description.
 *
 * Most of this implementation is taken from the MASCOT code. See:
 * https://github.com/bristolcrypto/SPDZ-2/blob/master/Math/gf2nlong.h
 *
 * We only need a subset of the operations:
 *
 *   =, !=, +, -, /, *, *=, +=
 *
 * and some sort of printing.
 */

#include <iostream>
#include <cstring>
#include <immintrin.h>
#include <wmmintrin.h>
#include <xmmintrin.h>
#include <emmintrin.h>

namespace goat {

typedef uint8_t byte;

/*
 * This class was copied without without modification from gf2nlong.h
 */
class int128 {
    public:

    __m128i a;

    int128() : a(_mm_setzero_si128()) { }
    int128(const int128& a) : a(a.a) { }
    int128(const __m128i& a) : a(a) { }
    int128(const unsigned long& a) : a(_mm_cvtsi64_si128(a)) { }
    int128(const unsigned long& upper, const unsigned long& lower)
	: a(_mm_set_epi64x(upper, lower)) { }

    unsigned long getLower() { return (unsigned long)_mm_cvtsi128_si64(a); }

    bool operator==(const int128& other) const {
	return _mm_test_all_zeros(a ^ other.a, a ^ other.a);
    }

    bool operator!=(const int128& other) const {
	return !(*this == other);
    }

    constexpr int128& operator=(const int128& other) = default;

    int128 operator<<(const int& other) const;
    int128 operator>>(const int& other) const;

    int128 operator^(const int128& other) const { return a ^ other.a; }
    int128 operator|(const int128& other) const { return a | other.a; }
    int128 operator&(const int128& other) const { return a & other.a; }

    int128 operator~() const { return ~a; }

    int128& operator<<=(const int& other) { return *this = *this << other; }
    int128& operator>>=(const int& other) { return *this = *this >> other; }

    int128& operator^=(const int128& other) { a ^= other.a; return *this; }
    int128& operator|=(const int128& other) { a |= other.a; return *this; }
    int128& operator&=(const int128& other) { a &= other.a; return *this; }

    friend std::ostream& operator<<(std::ostream& s, const int128& a);
};

class GF2_128 {
    private:

    int128 a;

    // TODO: figure which is needed and which are not.
    static int n, t1, t2, t3, nterms;
    static int l0, l1, l2, l3;
    static int128 mask, lowermask, uppermask;
    static bool rewind;

    void reduceTrinomial(int128 xh, int128 xl);
    void reducePentanomial(int128 xh, int128 xl);

    public:

    void reduce(int128 xh, int128 xl) {
	    reducePentanomial(xh, xl);
    }

    /*
     * We only care about, and support, degree = 128.
     */
    static void initField();

    int degree() const { return 128; }

    int getNterms() const { return nterms; }

    int getT(int i) {
	switch (i) {
	case 0: return t1;
	case 1: return t2;
	case 2: return t3;
	default: return -1;
	}
    }

    const std::string typeString() const { return "GF2_128"; }
    int size() const { return 128; }
    int t() const { return 0; } //???

    int128 get() const { return a; }
    __m128i getm128i() const { return a.a; }

    int is_zero() const { return a == int128(0); }
    int is_one() const { return a == int128(1); }

    /*
     * Assignment functions
     */
    void assign(const GF2_128 &x) { a = x.a; }
    void assignZero() { a = _mm_setzero_si128(); }
    void assignOne() { a = int128(0, 1); }
    void assignX() { a = int128(0, 2); }
    void assign(int128 aa) { a = aa & mask; }
    void assign(const int aa) {
	a = int128(static_cast<unsigned int>(aa)) & mask;
    }
    void assign(const char* buffer) { a = _mm_loadu_si128((__m128i*)buffer); }

    int getBit(int i) const {
	return ((a >> i) & i).getLower();
    }

    void setBit(int i, byte b) {
	if (b == 1)
	    a |= (1UL << i);
	else
	    a &= ~(1UL << i);
    }

    /*
     * Constructors
     */
    GF2_128() { assignZero(); }
    GF2_128(const GF2_128 &x) { assign(x); }
    GF2_128(const int128 &x) { assign(x); }
    GF2_128(const int x) { assign(x); }

    ~GF2_128() = default;

    constexpr GF2_128& operator=(const GF2_128& other) = default;

    /*
     * Operations
     */
    bool isZero() const { return a == int128(0); }
    bool equal(const GF2_128 &other) const { return a == other.a; }
    bool operator==(const GF2_128 &other) const { return equal(other); }
    bool operator!=(const GF2_128 &other) const { return !equal(other); }

    void add(const GF2_128 &x) { a ^= x.a; }

    GF2_128 operator+(const GF2_128 &x) const { return GF2_128(a ^ x.a); }
    GF2_128 operator-(const GF2_128 &x) const { return GF2_128(a ^ x.a); }
    GF2_128 &operator+=(const GF2_128 &x) {
	    add(x);
	    return *this;
    }

    GF2_128 &mul(const GF2_128 &x, const GF2_128 &y);
    GF2_128 operator*(const GF2_128 &x) {
	    GF2_128 r;
	    return r.mul(*this, x);
	    return r;
    }
    GF2_128 operator*=(const GF2_128 &x) {
	GF2_128 y = mul(*this, x);
	assign(y);
	return *this;
    }

    GF2_128 operator/(const GF2_128 &x) {
	    GF2_128 a;
	    a.invert(x);
	    return mul(*this, a);
    }

    void invert(const GF2_128 &a) {
	    *this = a;
	    invert();
    }
    void invert();

    friend std::ostream& operator<<(std::ostream &s,const GF2_128 &x) {
	s << std::hex << x.a << std::dec;
	return s;
    }
};

inline void mul128(__m128i a, __m128i b, __m128i *res1, __m128i *res2) {
    __m128i tmp3, tmp4, tmp5, tmp6;

    tmp3 = _mm_clmulepi64_si128(a, b, 0x00);
    tmp4 = _mm_clmulepi64_si128(a, b, 0x10);
    tmp5 = _mm_clmulepi64_si128(a, b, 0x01);
    tmp6 = _mm_clmulepi64_si128(a, b, 0x11);

    tmp4 = _mm_xor_si128(tmp4, tmp5);
    tmp5 = _mm_slli_si128(tmp4, 8);
    tmp4 = _mm_srli_si128(tmp4, 8);
    tmp3 = _mm_xor_si128(tmp3, tmp5);
    tmp6 = _mm_xor_si128(tmp6, tmp4);
    // initial mul now in tmp3, tmp6
    *res1 = tmp3;
    *res2 = tmp6;
}

}
#endif /* _HYPERBMR_GF2_128_H */
