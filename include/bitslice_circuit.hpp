#include "circuit.hpp"

namespace goat {

template<Bits>
void BitSliceCircuit::BitSliceCircuit(const size_t num_inputs) {

    const auto nin = num_inputs;
    const auto nvals = nin / Bits;

    // Assign inputs
    this->_inputs.resize(nin);
    for (size_t i = 0; i < nin; i++)
	this->_inputs[i] = INPUT(i);

    // Compute a bitslice circuit according to ...

    vector<vector<size_t>> v (Bits);
    for (size_t i = 0; i < Bits; i++) {
	v[i].resize(nvals);
	for (size_t j = 0; j < n; j++) {
	    v[i][j] = this->_inputs[i + Bits*j];
	}
    }

    vector<size_t> s = v[Bits-1];
    size_t bmax = fold_OR(V[Bits-1]);

    // will contain a 1 at the position of the maximum input.
    vector<size_t> w (nvals);
    for (size_t i = 0; i < n; i++)
	w[i] = INV(AND(bmax, INV(s[i])));

    for (ssize_t i = Bits - 2; i >= 0; i--) {
	s = vec_AND(w, v[i]);
	bmax = fold_OR(s);
	w = select(bmax, s, w);
    }

    // assign outputs, which at this point is the bits in w.
    this->_outputs.resize(nvals);
    for (size_t i = )

}

template<Bits>
vector<size_t> BitSliceCircuit::select(
    const size_t b,
    const vector<size_t> xs,
    const vector<size_t> ys)
{

    const auto n = xs.size();

    // verify that both arguments are of the same size
    assert (n == ys.size());

    vector<size_t> ws (n);
    for (size_t i = 0; i < n; i++)
	ws[i] = XOR(AND(b, XOR(xs[i], ys[i])), ys[i]);

    return ws;
}

template<Bits>
size_t BitSliceCircuit::fold_OR(const vector<size_t> xs) {
    const auto n = xs.size();

    if (n == 1)
	return xs[0];

    size_t w = OR(xs[0], xs[1]);
    for (size_t i = 1; i < n; i++)
	w = OR(w, xs[i]);
    return w;
}

template<Bits>
vector<size_t> BitSliceCircuit::vec_AND(
    const vector<size_t> us,
    const vector<size_t> vs)
{
    const auto n = us.size();

    assert (n == vs.size());

    vector<size_t> ws (n);
    for (size_t i = 0; i < n; i++)
	ws[i] = AND(us[i], vs[i]);
    return ws;
}

}
