#ifndef _GOAT_CIRCUIT_HPP
#define _GOAT_CIRCUIT_HPP

#include <vector>
#include <string>
#include <sstream>
#include <cassert>

#include <iostream>

// used by the printing function inside Circuit to control the amount of
// information to print.
#define SHOW_CIRC    1
#define SHOW_INPUTS  2
#define SHOW_OUTPUTS 4
#define SHOW_GATES   8

namespace goat {

using std::cout;

using std::vector;
using std::string;

enum class GateType {
    AND,
    XOR,
    INV,
    IN,
    OUT,
};

static inline string GateType_as_string(GateType gt) {
    switch (gt) {
    case GateType::AND:
	return "AND";
    case GateType::XOR:
	return "XOR";
    case GateType::INV:
	return "INV";
    case GateType::IN:
	return "IN";
    case GateType::OUT:
	return "OUT";
    }
    // Never reached, Just here to make GCC shut up :-)
    throw -1;
}

// Since we only care about AND, XOR and INV gates, each gate has at most fan-in
// 2 and always fan-out 1.
struct Gate {
    size_t u;  // left input index
    size_t v;  // right input index
    size_t w;  // output index

    GateType type;

    // used during evaluation
    int val = -1;

    // string representation of the gate
    string to_string() const {
	std::stringstream ss;
	ss << GateType_as_string(type);
	ss << " w=" << w << ", u=" << u;
	if (type == GateType::AND || type == GateType::XOR)
	    ss << ", v=" << v;
	ss << ", val=" << val;
	return ss.str();
    };

};

// Circuit defines a circuit of some number of inputs and outputs.
class Circuit {
public:

    Circuit() {};

    Circuit(const vector<Gate>& gates);

    size_t IN_count() const {
	return _number_of_inputs;
    };

    size_t OUT_count() const {
	return _number_of_outputs;
    };

    // Returns number of AND gates in the circuit
    size_t AND_count() const {
	return _AND_gate_count;
    };

    // ... XOR gates
    size_t XOR_count() const {
	return _XOR_gate_count;
    };

    // ... and INV gates
    size_t INV_count() const {
	return _INV_gate_count;
    };

    // Gates are inserted in topological order. See eval for how evaluation
    // would proceed.
    const vector<Gate>& gates() const {
	return _gates;
    };

    const vector<size_t>& get_outputs() const {
	return _outputs;
    };

    const vector<size_t>& get_inputs() const {
	return _inputs;
    };

    // Helper/debugger. Prints a bunch of info about the circuit. Arguments can
    // be used to control what is printed and what is not. By default only meta
    // information (e.g., number of AND gates) is displayed.
    void print_circuit_info(int opts) const;

    // Evaluate the circuit
    vector<bool> eval(const vector<bool>& inputs) const;

protected:

    // each of the functions below inserts a gate with input labels set
    // according to its arguments and with _current_max_gate_idx as the output
    // wire label. _current_max_gate_idx is incremented afterwards such that
    // calling
    //
    //   auto w = ADD(u, v)
    //
    // means that w is the wire label for the output of the just added AND gate
    // and _current_max_gate_idx now points to the next "free" label.

    size_t AND(const size_t u, const size_t v);
    size_t XOR(const size_t u, const size_t v);
    size_t INV(const size_t u);
    size_t INPUT(const size_t u);
    size_t OUTPUT(size_t u);

    inline size_t OR(const size_t u, const size_t v) {
	// Since we do not have OR as a primitive, we express it in terms of AND
	// and INV's using de-morgan's law:
	//
	//   X \or Y == \not (\not X \and \not Y)
	return INV(AND(INV(u), INV(v)));
    };

    size_t _current_max_gate_idx = 0;

    vector<Gate> _gates;
    vector<size_t> _inputs;
    vector<size_t> _outputs;

    size_t _AND_gate_count = 0;
    size_t _XOR_gate_count = 0;
    size_t _INV_gate_count = 0;

    size_t _number_of_inputs;
    size_t _number_of_outputs;
};



// Circuit which computes an argmax over it's inputs and gives the output to all
// parties.
// Each input has Bits bits.
template<size_t Bits>
class BitSliceCircuit : public Circuit {
public:

    BitSliceCircuit(const size_t num_inputs);

protected:

    void wire_internals();

private:

    using Circuit::AND;
    using Circuit::OR;
    using Circuit::INV;
    using Circuit::XOR;
    using Circuit::INPUT;
    using Circuit::OUTPUT;

    // ws[i] = xs[i] if b == 1 and ys[i] otherwise
    vector<size_t> select(const size_t b, const vector<size_t> xs, const vector<size_t> ys);

    // w = xs[0] \or xs[1] \or ... \or xs[n-1]
    size_t fold_OR(const vector<size_t> xs);

    // ws[i] = us[i] \and vs[i]
    vector<size_t> vec_AND(const vector<size_t> us, const vector<size_t> vs);
};

template<size_t Bits>
void BitSliceCircuit<Bits>::wire_internals() {
    const size_t m = this->_number_of_inputs;
    const size_t n = m / Bits;

    this->_inputs.resize(m);
    for (size_t i = 0; i < m; i++) {
	this->_inputs[i] = i;
	INPUT(i);
    }

    vector<vector<size_t>> V(Bits);
    for (size_t i = 0; i < Bits; i++) {
	V[i].resize(n);
	for (size_t j = 0; j < n; j++)
	    V[i][j] = this->_inputs[i + Bits*j];
    }

    vector<size_t> S = V[Bits-1];
    size_t bmax = fold_OR(V[Bits-1]);

    vector<size_t> W (n);
    for (size_t i = 0; i < n; i++)
	W[i] = INV(AND(bmax, INV(S[i])));

    for (ssize_t i = Bits - 2; i >= 0; i--) {
	S = vec_AND(W, V[i]);
	bmax = fold_OR(S);
	W = select(bmax, S, W);
    }

    this->_outputs.resize(n);
    for (size_t i = 0; i < n; i++)
	this->_outputs[i] = OUTPUT(W[i]);
}


template<size_t Bits>
class Mod2kAdderCircuit : public Circuit {
public:

    Mod2kAdderCircuit(const size_t num_inputs) :
	Circuit(num_inputs * 2 * Bits, num_inputs * Bits)
	{
	    wire_internals();
	};

#ifdef TESTING
    vector<bool> eval(const vector<vector<bool>>& inputs) const;
#endif

protected:
    void wire_internals();

private:

    using Circuit::AND;
    using Circuit::OR;
    using Circuit::INV;
    using Circuit::XOR;
    using Circuit::INPUT;
    using Circuit::OUTPUT;

};

template<size_t Bits>
void Mod2kAdderCircuit<Bits>::wire_internals() {
    // number of inputs as k-bit integers
    const size_t m = this->_number_of_inputs;

    // number of outputs as k-bit integers.
    const size_t n = this->_number_of_outputs;

    // number of k-bit integers we'll end up with
    const size_t kn = n/Bits;

    assert (m/2 == n);

    this->_inputs.resize(m);
    for (size_t i = 0; i < m; i++) {
	this->_inputs[i] = INPUT(i);
    }

    vector<vector<size_t>> xs (kn);
    vector<vector<size_t>> ys (kn);
    vector<vector<size_t>> zs (kn);
    vector<size_t> carries (kn);

    for (size_t i = 0; i < kn; i++) {
	xs[i].resize(Bits);
	ys[i].resize(Bits);
	zs[i].resize(Bits);
	for (size_t j = 0; j < Bits; j++) {
	    xs[i][j] = this->_inputs[j + i*Bits];
	    ys[i][j] = this->_inputs[j + i*Bits + n];
	}
    }

    // half-adder for first bit.
    for (size_t i = 0; i < kn; i++)
	carries[i] = AND(xs[i][0], ys[i][0]);
    for (size_t i = 0; i < kn; i++)
	zs[i][0] = XOR(xs[i][0], ys[i][0]);

    // full-adder for the remaining k-1 bits.
    for (size_t i = 0; i < kn; i++) {
	for (size_t j = 1; j < Bits; j++) {
	    const auto x = xs[i][j];
	    const auto y = ys[i][j];
	    const auto c = carries[i];

	    auto s = XOR(x, y);

	    zs[i][j] = XOR(s, c);
	    carries[i] = OR(AND(x, y), AND(s, c));
	}
    }

    this->_outputs.resize(n);
    for (size_t i = 0; i < kn; i++)
	for (size_t j = 0; j < Bits; j++)
	    this->_outputs[i*Bits + j] = OUTPUT(zs[i][j]);

}

#ifdef TESTING

template<size_t Bits>
vector<bool> Mod2kAdderCircuit<Bits>::eval(const vector<vector<bool>>& inputs) const {
    const auto n = this->_number_of_inputs;
    const auto m = n / Bits;
    const auto kn = m / 2;

    assert (inputs.size() == m);

    cout << "n=" << n << ", m=" << m << ", kn=" << kn << "\n";

    vector<bool> inputs_flat (n);
    size_t k = 0;
    for (size_t i = 0; i < kn; i++) {
	for (size_t j = 0; j < Bits; j++) {
	    inputs_flat[k] = inputs[i][j];
	    inputs_flat[k + kn*Bits] = inputs[i + kn][j];
	    k++;
	}
    }

    cout << "evaluating circuit ...\n";

    auto gates = this->gates();
    for (size_t i = 0; i < gates.size(); i++) {
	auto &gate = gates[i];
	switch (gate.type) {
    	case GateType::IN:
    	    gate.val = inputs_flat[gate.u];
	    break;
    	case GateType::OUT:
	    gate.val = gates[gate.u].val;
	    break;
	case GateType::AND:
	    gate.val = ((bool)(gates[gate.u].val)) and ((bool)(gates[gate.v].val));
	    break;
	case GateType::XOR:
	    gate.val = ((bool)(gates[gate.u].val)) xor ((bool)(gates[gate.v].val));
	    break;
	case GateType::INV:
	    gate.val = not ((bool)(gates[gate.u].val));
	    break;
	default:
	    throw -1;
	}
    }

    cout << "setting outputs ...\n";

    vector<bool> outputs (this->_number_of_outputs);
    for (size_t i = 0; i < this->_number_of_outputs; i++) {
	outputs[i] = gates[this->_outputs[i]].val;
    }

    return outputs;
}

template<size_t Bits>
vector<bool> BitSliceCircuit<Bits>::eval(const vector<vector<bool>> &inputs) const {
    const auto n = this->_number_of_inputs;
    const auto m = n / Bits;

    assert (inputs.size() == m);

    cout << "n=" << n << ", m=" << m << "\n";

    vector<bool> inputs_flat (n);
    for (size_t i = 0; i < m; i++) {
	for (size_t j = 0; j < Bits; j++) {
	    inputs_flat[i*Bits + j] = inputs[i][j];
	}
    }

    cout << "evaluating circuit ...\n";

    auto gates = this->gates();
    for (size_t i = 0; i < gates.size(); i++) {
	auto &gate = gates[i];

	switch (gate.type) {
    	case GateType::IN:
    	    gate.val = inputs_flat[gate.u];
	    break;
    	case GateType::OUT:
	    gate.val = gates[gate.u].val;
	    break;
	case GateType::AND:
	    gate.val = ((bool)(gates[gate.u].val)) and ((bool)(gates[gate.v].val));
	    break;
	case GateType::XOR:
	    gate.val = ((bool)(gates[gate.u].val)) xor ((bool)(gates[gate.v].val));
	    break;
	case GateType::INV:
	    gate.val = not ((bool)(gates[gate.u].val));
	    break;
	default:
	    throw -1;
	}
    }

    cout << "setting outputs ...\n";

    vector<bool> outputs (m);
    for (size_t i = 0; i < m; i++)
	outputs[i] = gates[this->_outputs[i]].val;

    return outputs;
}

#endif


}

#endif // _GOAT_CIRCUIT_HPP
