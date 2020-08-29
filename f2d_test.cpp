#include <iostream>
#include "f2d.hpp"
#include "crypto.hpp"
#include "ztk.hpp"

using namespace goat;
using namespace std;
using namespace ztk;

int main() {

    for (size_t i = 0; i < 16; i++) {
	for (size_t j = 0; j < 16; j++) {

	    const F2d<4> a {(elem_t)i};
	    const F2d<4> b {(elem_t)j};

	    cout << a << ", " << b << "\n";
	    cout << "a+b=" << a+b << ", a-b=" << a-b << "\n";
	    cout << "a*b=" << a*b << "\n";

	    if (b != 0)
	    	cout << "a/b=" << a/b << "\n";
	}
    }

    cout << "\n\n\n";

    for (size_t i = 0; i < 10; i++) {
	GR<64, 4> x = randomize<GR<64, 4>>();
	cout << x << "\n";
	F2d<4> y {x};
	cout << y << "\n\n";
    }

}
