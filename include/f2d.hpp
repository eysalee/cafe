#ifndef _F2D_HPP
#define _F2D_HPP

#include <sstream>

#include "ztk.hpp"

namespace goat {

using ztk::GR;

constexpr size_t power_of_2(const size_t k) {
    return k == 0 ? 1 : 2 * power_of_2(k - 1);
}

typedef unsigned char elem_t;

template<size_t D>
class F2d {
public:
    F2d(const elem_t e);

    template<size_t K, size_t D_>
    F2d(const GR<K, D_>& x) {
	static_assert (D == D_);
	init_tables();

	e = 0;
	for (ssize_t i = D-1; i >= 0; i--) {
	    e ^= (x[i].is_odd());
	    if (i > 0)
		e <<= 1;
	}
    }

    F2d<D> operator+(const F2d<D>& x) const {
	return F2d{(elem_t)(this->e ^ x.e)};
    };

    F2d<D> operator-(const F2d<D>& x) const {
	return F2d{(elem_t)(this->e ^ x.e)};
    };

    F2d<D> operator*(const F2d<D>& x) const {
	return F2d<D>{mtable[this->e][x.e]};
    };

    F2d<D> operator/(const F2d<D>& x) const {
	auto c = F2d<D>{itable[x.e]};
	return (*this) * c;
    };

    F2d<D> operator=(const F2d<D>& x) {
	this->e = x.e;
    };

    bool operator==(const F2d<D>& x) const {
	return this->e == x.e;
    };

    bool operator==(const int x) const {
	return this->e == ((elem_t)x);
    };

    bool operator!=(const F2d<D>& x) const {
	return !(*this == x);
    };

    bool operator!=(const int x) const {
	return this->e != ((elem_t)x);
    };

    F2d<D> invert() const {
	return F2d<D>{itable[this->e]};
    };

    elem_t get() const {
	return this->e;
    };

    template<size_t L>
    friend std::ostream& operator<<(std::ostream &os, const F2d<L>& x);

private:
    const static elem_t mask = power_of_2(D) - 1;
    elem_t e;

    static elem_t mtable[mask + 1][mask + 1];
    static elem_t itable[mask + 1];
    static bool init_done;
    static void init_tables();
};

template<size_t U, size_t P>
inline elem_t mul(const elem_t x, const elem_t y) {
    elem_t a = x;
    elem_t b = y;
    elem_t c = 0;

    while (a && b) {
	if (b & 1)
	    c ^= a;
	if (a & U)
	    a = (a << 1) ^ P;
	else
	    a <<= 1;
	b >>= 1;
    }
    return c;
}

template<size_t L>
std::ostream& operator<<(std::ostream &os, const F2d<L> &x) {
    os << "{" << (int)(x.e) << "}";
    return os;
}

}

#endif // _F2D_HPP
