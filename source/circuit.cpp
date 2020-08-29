#include "circuit.hpp"

namespace goat {

Circuit::Circuit(const vector<Gate>& gates)
    : _gates{gates}
{
    // figure out inputs, outputs and how many of each gates we have.
    for (auto &gate: _gates) {
	_current_max_gate_idx++;

	switch (gate.type) {
	case GateType::AND:
	    _AND_gate_count++;
	    break;
	case GateType::XOR:
	    _XOR_gate_count++;
	    break;
	case GateType::INV:
	    _INV_gate_count++;
	    break;
	case GateType::IN:
	    _inputs.emplace_back(gate.u);
	    _number_of_inputs++;
	    break;
	case GateType::OUT:
	    _outputs.emplace_back(gate.u);
	    _number_of_outputs++;
	    break;
	}
    }
}

void Circuit::print_circuit_info(int opts) const {

    if (opts & SHOW_CIRC) {
	cout << "number of inputs: " << _number_of_inputs << "\n"
	     << "number of outputs: " << _number_of_outputs << "\n"
	     << "#AND gates: " << _AND_gate_count << "\n"
	     << "#XOR gates: " << _XOR_gate_count << "\n"
	     << "#INV gates: " << _INV_gate_count << "\n";
    }

    if (opts & SHOW_INPUTS) {
	cout << "-----inputs----\n";

	for (size_t i = 0; i < _number_of_inputs; i++) {
	    cout << _inputs[i] << " ";
	}
	cout << "\n";
    }

    if (opts & SHOW_OUTPUTS) {
	cout << "-----outputs----\n";

	for (size_t i = 0; i < _number_of_outputs; i++) {
	    cout << _outputs[i] << " ";
	}
	cout << "\n";
    }

    if (opts & SHOW_GATES) {

	cout << "-----gates----\n";
	size_t c = 0;
	for (auto &gate: _gates)
	    cout << c++ << " " << gate.to_string() << "\n";
    }
}

vector<bool> Circuit::eval(const vector<bool>& inputs) const {

    assert (inputs.size() == IN_count());

    vector<bool> outputs (OUT_count());

#define CHECK_GATE(gate) do {					\
	if (!((gate).val == 0 || (gate).val == 1)) {		\
	    cout << "error @ " << (gate).to_string() << "\n";	\
	    throw -1;						\
	}							\
    } while (0)



    auto gates = _gates;
    Gate gu, gv;

    for (size_t i = 0; i < gates.size(); i++) {
	Gate& gate = gates[i];

	switch (gate.type) {

	case GateType::IN:
	    assert (gate.u < IN_count());
	    gate.val = inputs[gate.u];
	    break;

	case GateType::OUT:
	    gate.val = gates[gate.u].val;
	    assert (gate.w < OUT_count());
	    outputs[gate.w] = gate.val;
	    break;

	case GateType::AND:
	    gu = gates[gate.u];
	    gv = gates[gate.v];
	    CHECK_GATE(gu);
	    CHECK_GATE(gv);
	    gate.val = ((bool)gu.val) and ((bool)gv.val);
	    break;

	case GateType::XOR:
	    gu = gates[gate.u];
	    gv = gates[gate.v];
	    CHECK_GATE(gu);
	    CHECK_GATE(gv);
	    gate.val = ((bool)gu.val) xor ((bool)gv.val);
	    break;

	case GateType::INV:
	    gu = gates[gate.u];
	    CHECK_GATE(gu);
	    gate.val = not ((bool)gu.val);
	    break;
	}
    }

#undef CHECK_GATE

    return outputs;
}

// protected methods

size_t Circuit::AND(const size_t u, const size_t v) {
    auto w = _current_max_gate_idx;
    _gates.emplace_back(Gate{u, v, w, GateType::AND});
    _current_max_gate_idx++;
    _AND_gate_count++;
    return w;
}

size_t Circuit::XOR(const size_t u, const size_t v) {
    size_t w = _current_max_gate_idx;
    _gates.emplace_back(Gate{u, v, w, GateType::XOR});
    _current_max_gate_idx++;
    _XOR_gate_count++;
    return w;
}

size_t Circuit::INV(const size_t u) {
    size_t w = _current_max_gate_idx;
    _gates.emplace_back(Gate{u, 0, w, GateType::INV});
    _current_max_gate_idx++;
    _INV_gate_count++;
    return w;
}

size_t Circuit::INPUT(const size_t u) {
    size_t w = _current_max_gate_idx;
    _gates.emplace_back(Gate{u, 0, w, GateType::IN});
    _current_max_gate_idx++;
    return w;
}

size_t Circuit::OUTPUT(const size_t u) {
    size_t w = _current_max_gate_idx;
    _gates.emplace_back(Gate{u, 0, w, GateType::OUT});
    _current_max_gate_idx++;
    return w;
}

} // goat
