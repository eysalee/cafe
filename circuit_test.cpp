#include "circuit.hpp"
#include "crypto.hpp"

using namespace goat;
using namespace std;

vector<bool> i2b(const uint64_t x) {
    vector<bool> r (64);
    auto y = x;
    for (size_t i = 0; i < 64; i++) {
	r[i] = (bool)(y & 1);
	y >>= 1;
    }
    return r;
}

typedef uint64_t u64;

size_t argmax(const vector<u64> &vals) {
    u64 max = vals[0];
    size_t idx = 0;

    for (size_t i = 1; i < vals.size(); i++) {
	if (vals[i] > max) {
	    max = vals[i];
	    idx = i;
	}
    }
    return idx;
}

void test_BitSliceCircuit(const size_t n) {

    BitSliceCircuit<64> c (n);

    vector<u64> vals (n);
    vector<vector<bool>> inputs (n);
    for (size_t i = 0; i < n; i++) {
	auto v = randomize<u64>();
	inputs[i] = i2b(v);
	vals[i] = v;
	for (size_t j = 0; j < 64; j++)
	    cout << inputs[i][j];
	cout << " v=" << v << "\n";
    }
    size_t computed_argmax = argmax(vals);
    cout << "argmax=" << computed_argmax << "\n";

    auto outputs = c.eval(inputs);

    bool all_good = false;
    for (size_t i = 0; i < outputs.size(); i++) {
	if (outputs[i] == 1 && all_good) {
	    cout << "conflicting index: " << i << "\n";
	}
	if (outputs[i] == 1 && computed_argmax == i) {
	    all_good = true;
	}
    }
    c.print_circuit_info(SHOW_CIRC | SHOW_OUTPUTS);
    cout << "test good? " << all_good << "\n";

}

void test_Mod2kAdderCircuit(const size_t n) {

    vector<u64> xs (n);
    vector<u64> ys (n);
    vector<u64> zs (n);

    vector<vector<bool>> inputs (2*n);

    for (size_t i = 0; i < n; i++) {
	auto x = randomize<u64>();
	auto y = randomize<u64>();
	xs[i] = x;
	ys[i] = y;
	zs[i] = x + y;
	inputs[i] = i2b(x);
	inputs[i + n] = i2b(y);
    }

    Mod2kAdderCircuit<64> c (n);
    c.print_circuit_info(SHOW_CIRC);

    auto outputs = c.eval(inputs);

    bool all_good = true;
    for (size_t i = 0; i < n; i++) {
	u64 v = 0;
	for (ssize_t j = 63; j >= 0; j--) {
	    const auto bit = outputs[i*64 + j];
	    v <<= 1;
	    v |= bit;
	}
	if (zs[i] != v) {
	    cout << "uh no:\n"
		 << "v=" << v << "\n"
		 << "zs[" << i << "]=" << zs[i] << "\n\n";
	    all_good = false;
	}
    }

    if (all_good)
	cout << "test passsed :-)\n";
}

void test_Soldered() {

}

int main(int argc, char **argv) {

    string s;

    if (argc > 1)
	s = string(argv[1]);
    else
	s = "";

    if (s == "bitslice") {
	test_BitSliceCircuit(1000);
    } else if (s == "mod2kadder") {
	test_Mod2kAdderCircuit(2000);
    } else {
	test_Soldered();
    }
}
