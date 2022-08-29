#define _USE_MATH_DEFINES
#include "parametric_circuit.hpp"

#include <cppsim/exception.hpp>
#include <cppsim/gate_factory.hpp>
#include <cppsim/gate_matrix.hpp>
#include <cppsim/gate_merge.hpp>
#include <cppsim/state.hpp>
#include <cppsim/type.hpp>
#include <iostream>

#include "parametric_gate.hpp"
#include "parametric_gate_factory.hpp"

ParametricQuantumCircuit::ParametricQuantumCircuit(
    UINT qubit_count_, std::string style)
    : QuantumCircuit(qubit_count_) {
    if (style == "undefined") {
        _specified_old = false;
        _specified_new = false;
    } else if (style == "old") {
        _specified_old = true;
        _specified_new = false;
    } else if (style == "new") {
        _specified_old = false;
        _specified_new = true;
    } else if (style ==
               "I_understand_that_mixing_ParameterType_is_very_dangerous_but_I_"
               "still_want_to_do_that") {
        _specified_old = true;
        _specified_new = true;
    } else {
        throw InvalidParametricQuantumCircuitStyleOptionException(
            "Error: "
            "ParametricQuantumCircuit::ParametricQuantumCircuit(UINT, "
            "std::string): "
            "unknown style option: " +
            style);
    }
}
ParametricQuantumCircuit::ParametricQuantumCircuit(
    UINT qubit_count_, const ParameterSet& parameter_set, std::string style)
    : QuantumCircuit(qubit_count_), _parameter_set(parameter_set) {
    if (style == "undefined") {
        _specified_old = false;
        _specified_new = false;
    } else if (style == "old") {
        _specified_old = true;
        _specified_new = false;
    } else if (style == "new") {
        _specified_old = false;
        _specified_new = true;
    } else if (style ==
               "I_understand_that_mixing_ParameterType_is_very_dangerous_but_I_"
               "still_want_to_do_that") {
        _specified_old = true;
        _specified_new = true;
    } else {
        throw InvalidParametricQuantumCircuitStyleOptionException(
            "Error: "
            "ParametricQuantumCircuit::ParametricQuantumCircuit(UINT, "
            "std::string): "
            "unknown style option: " +
            style);
    }
}

bool ParametricQuantumCircuit::is_old_style() const {
    if (_specified_old) return true;
    if (_specified_new) return false;
    if (!_parameter_set.empty()) return false;
    if (_parametric_gate_list.size() == 0) return true;
    if (_parametric_gate_list[0]->get_parameter_type() == "gate") return true;
    return false;
}
bool ParametricQuantumCircuit::is_new_style() const {
    if (_specified_new) return true;
    if (_specified_old) return false;
    if (!_parameter_set.empty()) return true;
    if (_parametric_gate_list.size() == 0) return true;
    if (_parametric_gate_list[0]->get_parameter_type() == "id") return true;
    return false;
}

ParametricQuantumCircuit* ParametricQuantumCircuit::copy() const {
    ParametricQuantumCircuit* new_circuit =
        new ParametricQuantumCircuit(this->qubit_count, this->_parameter_set);
    for (UINT gate_pos = 0; gate_pos < this->gate_list.size(); gate_pos++) {
        new_circuit->add_gate(this->gate_list[gate_pos]->copy());
    }
    for (UINT gate_pos : this->_parametric_gate_position) {
        new_circuit->_parametric_gate_list.push_back(
            dynamic_cast<QuantumGate_SingleParameter*>(
                new_circuit->gate_list[gate_pos]));
    }
    return new_circuit;
}

void ParametricQuantumCircuit::create_parameter(
    const ParameterId& parameter_id, double initial_parameter) {
    if (!this->is_new_style()) {
        throw NotImplementedException(
            "Error: "
            "ParametricQuantumCircuit::create_parameter(const "
            "ParameterId&, double) : "
            "this is a new-style function. If you want to use this, do not add "
            "gate with an old-style parameteric_gate whose parameter is "
            "inside.");
    }
    if (_parameter_set.count(parameter_id) > 0) {
        throw ParameterIdDuplicatedException(
            "Error: "
            "ParametricQuantumCircuit::create_parameter(const ParameterKey&, "
            "double): "
            "parameter_id \"" +
            parameter_id + "\" already exists");
    }
    _parameter_set[parameter_id] = initial_parameter;
}
void ParametricQuantumCircuit::remove_parameter(
    const ParameterId& parameter_id) {
    if (!this->is_new_style()) {
        throw NotImplementedException(
            "Error: "
            "ParametricQuantumCircuit::remove_parameter(const ParameterId&): "
            "this is a new-style function. If you want to use this, do not add "
            "gate with an old-style parameteric_gate whose parameter is "
            "inside.");
    }
    if (_parameter_set.count(parameter_id) == 0) {
        throw ParameterIdNotFoundException(
            "Error: ParametricQuantumCircuit::create_parameter(const "
            "ParameterId&): parameter_id \"" +
            parameter_id + "\" is not found");
    }
    _parameter_set.erase(parameter_id);
}
bool ParametricQuantumCircuit::contains_parameter(
    const ParameterId& parameter_id) const {
    if (!this->is_new_style()) {
        throw NotImplementedException(
            "Error: "
            "ParametricQuantumCircuit::contains_paramter(const ParameterId&) "
            "const: "
            "this is a new-style function. If you want to use "
            "this, do not add gate with an old-style parameteric_gate whose "
            "parameter is inside.");
    }
    return _parameter_set.count(parameter_id) > 0;
}
UINT ParametricQuantumCircuit::get_parameter_count() const {
    if (!this->is_old_style()) {
        throw NotImplementedException(
            "Error: "
            "ParametricQuantumCircuit::get_paramter_count() const: "
            "this is an old-style function. If you want to use this, do not "
            "add gate with a new-style parametric_gate which has string "
            "parameter_id.");
    }
    return (UINT)_parametric_gate_list.size();
}
UINT ParametricQuantumCircuit::get_parametric_gate_count() const {
    if (!this->is_new_style()) {
        throw NotImplementedException(
            "Error: "
            "ParametricQuantumCircuit::get_parametric_gate_count() const: "
            "this is a new-style function. If you want to use "
            "this, do not add gate with an old-style parameteric_gate whose "
            "parameter is inside.");
    }
    return (UINT)_parametric_gate_list.size();
}
UINT ParametricQuantumCircuit::get_parameter_id_count() const {
    if (!this->is_new_style()) {
        throw NotImplementedException(
            "Error: "
            "ParametricQuantumCircuit::get_parameter_id_count() const: "
            "this is a new-style function. If you want to use "
            "this, do not add gate with an old-style parameteric_gate whose "
            "parameter is inside.");
    }
    return (UINT)_parameter_set.size();
}
double ParametricQuantumCircuit::get_parameter(UINT index) const {
    if (!this->is_old_style()) {
        throw NotImplementedException(
            "Error: "
            "ParametricQuantumCircuit::get_parameter(UINT) const: "
            "this is an old-style function. If you want to use this, do not "
            "add gate with a new-style parametric_gate which has string "
            "parameter_id.");
    }
    if (index >= this->_parametric_gate_list.size()) {
        throw ParameterIndexOutOfRangeException(
            "Error: "
            "ParametricQuantumCircuit::get_parameter(UINT): "
            "parameter index is out of range");
    }
    return _parametric_gate_list[index]->get_parameter_value();
}
double ParametricQuantumCircuit::get_parameter(
    const ParameterId& parameter_id) const {
    if (!this->is_new_style()) {
        throw NotImplementedException(
            "Error: "
            "ParametricQuantumCircuit::get_paramter(const ParameterId&) const: "
            "this is a new-style function. If you want to use "
            "this, do not add gate with an old-style parameteric_gate whose "
            "parameter is inside.");
    }
    auto it = _parameter_set.find(parameter_id);
    if (it == _parameter_set.end()) {
        throw ParameterIdNotFoundException(
            "Error: "
            "ParametricQuantumCircuit::create_parameter(const "
            "ParameterId&): parameter_id \"" +
            parameter_id + "\" is not found");
    }
    return it->second;
}
double ParametricQuantumCircuit::get_angle(UINT index) const {
    if (!this->is_new_style()) {
        throw NotImplementedException(
            "Error: "
            "ParametricQuantumCircuit::get_angle(UINT index) const: "
            "this is a new-style function. If you want to use "
            "this, do not add gate with an old-style parameteric_gate whose "
            "parameter is inside.\n"
            "Info: use get_parameter(UINT parameter_index) instead.");
    }
    if (index >= this->_gate_list.size()) {
        throw GateIndexOutOfRangeException(
            "ParametricQuantumCircuit::get_angle(UINT): "
            "gate index is out of range");
    }
    if (!this->_gate_list[index]->is_parametric()) {
        throw NotImplementedException(
            "Error: "
            "ParametricQuantumCircuit::get_angle(UINT index) const: "
            "specified gate is not parametric.");
    }
    auto pgate =
        dynamic_cast<QuantumGate_SingleParameter*>(this->_gate_list[index]);
    return pgate->get_angle(this->_parameter_set);
}
void ParametricQuantumCircuit::set_parameter(UINT index, double value) {
    if (!this->is_old_style()) {
        throw NotImplementedException(
            "Error: "
            "ParametricQuantumCircuit::set_parameter(UINT, double): "
            "this is an old-style function. If you want to use this, do not "
            "add gate with a new-style parametric_gate which has string "
            "parameter_id.");
    }
    if (index >= this->_parametric_gate_list.size()) {
        throw ParameterIndexOutOfRangeException(
            "Error: ParametricQuantumCircuit::set_parameter(UINT,double): "
            "parameter index is out of range");
    }
    _parametric_gate_list[index]->set_parameter_value(value);
}
void ParametricQuantumCircuit::set_parameter(
    const ParameterId& parameter_id, double value) {
    if (!this->is_new_style()) {
        throw NotImplementedException(
            "Error: "
            "ParametricQuantumCircuit::set_parameter(const ParameterId&, "
            "double): "
            "this is a new-style function. If you want to use "
            "this, do not add gate with an old-style parameteric_gate whose "
            "parameter is inside.");
    }
    if (!_parameter_set.count(parameter_id)) {
        throw ParameterIdNotFoundException(
            "Error: "
            "ParametricQuantumCircuit::create_parameter(const "
            "ParameterId&): parameter_id \"" +
            parameter_id + "\" is not found");
    }
    _parameter_set[parameter_id] = value;
}
ParameterSet ParametricQuantumCircuit::get_parameter_set() const {
    if (!this->is_new_style()) {
        throw NotImplementedException(
            "Error: "
            "ParametricQuantumCircuit::get_parameter_set(): "
            "this is a new-style function. If you want to use "
            "this, do not add gate with an old-style parameteric_gate whose "
            "parameter is inside.");
    }
    return _parameter_set;
}
void ParametricQuantumCircuit::set_parameter_set(
    const ParameterSet& parameter_set) {
    if (!this->is_new_style()) {
        throw NotImplementedException(
            "Error: "
            "ParametricQuantumCircuit::set_parameter_set(): "
            "this is a new-style function. If you want to use "
            "this, do not add gate with an old-style parameteric_gate whose "
            "parameter is inside.");
    }
    _parameter_set = parameter_set;
}

std::string ParametricQuantumCircuit::to_string() const {
    std::stringstream os;
    os << QuantumCircuit::to_string();
    os << "*** Parameter Info ***" << std::endl;
    if (this->is_old_style()) {
        os << "# of parameter: " << this->get_parameter_count() << std::endl;
    } else {
        os << "# of parameter: " << this->get_parameter_id_count() << std::endl;
    }
    return os.str();
}

std::ostream& operator<<(
    std::ostream& stream, const ParametricQuantumCircuit& circuit) {
    stream << circuit.to_string();
    return stream;
}
std::ostream& operator<<(
    std::ostream& stream, const ParametricQuantumCircuit* gate) {
    stream << *gate;
    return stream;
}

UINT ParametricQuantumCircuit::get_parametric_gate_position(UINT index) const {
    if (index >= this->_parametric_gate_list.size()) {
        throw ParameterIndexOutOfRangeException(
            "Error: "
            "ParametricQuantumCircuit::get_parametric_gate_position(UINT): "
            "parameter index is out of range");
    }

    return _parametric_gate_position[index];
}
void ParametricQuantumCircuit::add_gate(QuantumGateBase* gate) {
    QuantumCircuit::add_gate(gate);
}
void ParametricQuantumCircuit::add_gate(QuantumGateBase* gate, UINT index) {
    QuantumCircuit::add_gate(gate, index);
    for (auto& val : _parametric_gate_position)
        if (val >= index) val++;
}
void ParametricQuantumCircuit::add_gate_copy(const QuantumGateBase* gate) {
    QuantumCircuit::add_gate(gate->copy());
}
void ParametricQuantumCircuit::add_gate_copy(
    const QuantumGateBase* gate, UINT index) {
    QuantumCircuit::add_gate(gate->copy(), index);
    for (auto& val : _parametric_gate_position)
        if (val >= index) val++;
}

void ParametricQuantumCircuit::remove_gate(UINT index) {
    auto ite = std::find(_parametric_gate_position.begin(),
        _parametric_gate_position.end(), (unsigned int)index);
    if (ite != _parametric_gate_position.end()) {
        UINT dist = (UINT)std::distance(_parametric_gate_position.begin(), ite);
        _parametric_gate_position.erase(
            _parametric_gate_position.begin() + dist);
        _parametric_gate_list.erase(_parametric_gate_list.begin() + dist);
    }
    QuantumCircuit::remove_gate(index);
    for (auto& val : _parametric_gate_position)
        if (val >= index) val--;
}
void ParametricQuantumCircuit::merge_circuit(
    const ParametricQuantumCircuit* circuit) {
    if ((!this->is_old_style() && !circuit->is_new_style()) ||
        (!this->is_new_style() && !circuit->is_old_style())) {
        throw NotImplementedException(
            "Error: "
            "ParametricQuantumCircuit::merge_circuit(const "
            "ParametricQuantumCircuit*): "
            "cannot merge circuit which have different style");
    }
    UINT gate_count = this->gate_list.size();
    for (auto gate : circuit->gate_list) {
        this->add_gate_copy(gate);
    }
    for (auto gate_position : circuit->_parametric_gate_position) {
        UINT new_gate_position = gate_position + gate_count;
        this->_parametric_gate_position.push_back(new_gate_position);
        this->_parametric_gate_list.push_back(
            dynamic_cast<QuantumGate_SingleParameter*>(
                this->gate_list[new_gate_position]));
    }
    for (auto& p : circuit->_parameter_set) {
        this->_parameter_set[p.first] = p.second;
    }
    return;
}

void ParametricQuantumCircuit::add_parametric_gate(
    QuantumGate_SingleParameter* gate) {
    if (!this->is_old_style() && gate->get_parameter_type() == "gate") {
        throw NotImplementedException(
            "Error: "
            "ParametricQuantumCircuit::add_parametric_gate(QuantumGate_"
            "SingleParameter*): "
            "cannot add gate-type parameter on a new-style parametric circuit");
    }
    if (!this->is_new_style() && gate->get_parameter_type() == "id") {
        throw NotImplementedException(
            "Error: "
            "ParametricQuantumCircuit::add_parametric_gate(QuantumGate_"
            "SingleParameter*): "
            "cannot add id-type parameter on an old-style parametric circuit");
    }
    this->add_gate(gate);
    _parametric_gate_position.push_back((UINT)gate_list.size() - 1);
    _parametric_gate_list.push_back(gate);
};
void ParametricQuantumCircuit::add_parametric_gate(
    QuantumGate_SingleParameter* gate, UINT index) {
    if (!this->is_old_style() && gate->get_parameter_type() == "gate") {
        throw NotImplementedException(
            "Error: "
            "ParametricQuantumCircuit::add_parametric_gate(QuantumGate_"
            "SingleParameter*, UINT): "
            "cannot add gate-type parameter on a new-style parametric circuit");
    }
    if (!this->is_new_style() && gate->get_parameter_type() == "id") {
        throw NotImplementedException(
            "Error: "
            "ParametricQuantumCircuit::add_parametric_gate(QuantumGate_"
            "SingleParameter*, UINT): "
            "cannot add id-type parameter on an old-style parametric circuit");
    }
    this->add_gate(gate, index);
    _parametric_gate_position.push_back(index);
    _parametric_gate_list.push_back(gate);
}
void ParametricQuantumCircuit::add_parametric_gate_copy(
    QuantumGate_SingleParameter* gate) {
    if (!this->is_old_style() && gate->get_parameter_type() == "gate") {
        throw NotImplementedException(
            "Error: "
            "ParametricQuantumCircuit::add_parametric_gate_copy(QuantumGate_"
            "SingleParameter*): "
            "cannot add gate-type parameter on a new-style parametric circuit");
    }
    if (!this->is_new_style() && gate->get_parameter_type() == "id") {
        throw NotImplementedException(
            "Error: "
            "ParametricQuantumCircuit::add_parametric_gate_copy(QuantumGate_"
            "SingleParameter*): "
            "cannot add id-type parameter on an old-style parametric circuit");
    }
    QuantumGate_SingleParameter* copied_gate = gate->copy();
    this->add_gate(copied_gate);
    _parametric_gate_position.push_back((UINT)gate_list.size() - 1);
    _parametric_gate_list.push_back(copied_gate);
};
void ParametricQuantumCircuit::add_parametric_gate_copy(
    QuantumGate_SingleParameter* gate, UINT index) {
    if (!this->is_old_style() && gate->get_parameter_type() == "gate") {
        throw NotImplementedException(
            "Error: "
            "ParametricQuantumCircuit::add_parametric_gate_copy(QuantumGate_"
            "SingleParameter*, UINT): "
            "cannot add gate-type parameter on a new-style parametric circuit");
    }
    if (!this->is_new_style() && gate->get_parameter_type() == "id") {
        throw NotImplementedException(
            "Error: "
            "ParametricQuantumCircuit::add_parametric_gate_copy(QuantumGate_"
            "SingleParameter*, UINT): "
            "cannot add id-type parameter on an old-style parametric circuit");
    }
    QuantumGate_SingleParameter* copied_gate = gate->copy();
    this->add_gate(copied_gate, index);
    _parametric_gate_position.push_back(index);
    _parametric_gate_list.push_back(copied_gate);
}

void ParametricQuantumCircuit::add_parametric_RX_gate(
    UINT target_index, double initial_angle) {
    this->add_parametric_gate(gate::ParametricRX(target_index, initial_angle));
}
void ParametricQuantumCircuit::add_parametric_RX_gate(
    UINT target_index, const ParameterId& parameter_id, double parameter_coef) {
    if (!this->contains_parameter(parameter_id)) {
        throw ParameterIdNotFoundException(
            "Error: "
            "ParametricQuantumCircuit::add_parametric_RX_gate(UINT, const "
            "ParameterId&, double): "
            "parameter_id \"" +
            parameter_id + "\" does not exists.");
    }
    this->add_parametric_gate(
        gate::ParametricRX(target_index, parameter_id, parameter_coef));
}
void ParametricQuantumCircuit::add_parametric_RX_gate_new_parameter(
    UINT target_index, const ParameterId& parameter_id, double value,
    double parameter_coef) {
    if (this->contains_parameter(parameter_id)) {
        throw ParameterIdDuplicatedException(
            "Error: "
            "ParametericQuantumCircuit::add_parametric_RX_gate_new_parameter("
            "UINT, const ParameterId&, double, double): "
            "parameter_id \"" +
            parameter_id + "\" already exists");
    }
    this->create_parameter(parameter_id, value);
    this->add_parametric_gate(
        gate::ParametricRX(target_index, parameter_id, parameter_coef));
}
void ParametricQuantumCircuit::add_parametric_RY_gate(
    UINT target_index, double initial_angle) {
    this->add_parametric_gate(gate::ParametricRY(target_index, initial_angle));
}
void ParametricQuantumCircuit::add_parametric_RY_gate(
    UINT target_index, const ParameterId& parameter_id, double parameter_coef) {
    if (!this->contains_parameter(parameter_id)) {
        throw ParameterIdNotFoundException(
            "Error: "
            "ParametricQuantumCircuit::add_parametric_RY_gate(UINT, const "
            "ParameterId&, double): "
            "parameter_id \"" +
            parameter_id + "\" does not exists.");
    }
    this->add_parametric_gate(
        gate::ParametricRY(target_index, parameter_id, parameter_coef));
}
void ParametricQuantumCircuit::add_parametric_RY_gate_new_parameter(
    UINT target_index, const ParameterId& parameter_id, double value,
    double parameter_coef) {
    if (this->contains_parameter(parameter_id)) {
        throw ParameterIdDuplicatedException(
            "Error: "
            "ParametericQuantumCircuit::add_parametric_RY_gate_new_parameter("
            "UINT, const ParameterId&, double, double): "
            "parameter_id \"" +
            parameter_id + "\" already exists");
    }
    this->create_parameter(parameter_id, value);
    this->add_parametric_gate(
        gate::ParametricRY(target_index, parameter_id, parameter_coef));
}
void ParametricQuantumCircuit::add_parametric_RZ_gate(
    UINT target_index, double initial_angle) {
    this->add_parametric_gate(gate::ParametricRZ(target_index, initial_angle));
}
void ParametricQuantumCircuit::add_parametric_RZ_gate(
    UINT target_index, const ParameterId& parameter_id, double parameter_coef) {
    if (!this->contains_parameter(parameter_id)) {
        throw ParameterIdNotFoundException(
            "Error: "
            "ParametricQuantumCircuit::add_parametric_RZ_gate(UINT, const "
            "ParameterId&, double): "
            "parameter_id \"" +
            parameter_id + "\" does not exists.");
    }
    this->add_parametric_gate(
        gate::ParametricRZ(target_index, parameter_id, parameter_coef));
}
void ParametricQuantumCircuit::add_parametric_RZ_gate_new_parameter(
    UINT target_index, const ParameterId& parameter_id, double value,
    double parameter_coef) {
    if (this->contains_parameter(parameter_id)) {
        throw ParameterIdDuplicatedException(
            "Error: "
            "ParametericQuantumCircuit::add_parametric_RZ_gate_new_parameter("
            "UINT, const ParameterId&, double, double): "
            "parameter_id \"" +
            parameter_id + "\" already exists");
    }
    this->create_parameter(parameter_id, value);
    this->add_parametric_gate(
        gate::ParametricRZ(target_index, parameter_id, parameter_coef));
}

void ParametricQuantumCircuit::add_parametric_multi_Pauli_rotation_gate(
    std::vector<UINT> target, std::vector<UINT> pauli_id,
    double initial_angle) {
    this->add_parametric_gate(
        gate::ParametricPauliRotation(target, pauli_id, initial_angle));
}
void ParametricQuantumCircuit::add_parametric_multi_Pauli_rotation_gate(
    std::vector<UINT> target, std::vector<UINT> pauli_id,
    const ParameterId& parameter_id, double parameter_coef) {
    if (!this->contains_parameter(parameter_id)) {
        throw ParameterIdNotFoundException(
            "Error: "
            "ParametricQuantumCircuit::add_parametric_multi_Pauli_rotation_"
            "gate(std::"
            "vector<UINT>, std::vector<UINT>, const ParameterId&, double): "
            "parameter_id \"" +
            parameter_id + "\" does not exists.");
    }
    this->add_parametric_gate(gate::ParametricPauliRotation(
        target, pauli_id, parameter_id, parameter_coef));
}
void ParametricQuantumCircuit::
    add_parametric_multi_Pauli_rotation_gate_new_parameter(
        std::vector<UINT> target, std::vector<UINT> pauli_id,
        const ParameterId& parameter_id, double value, double parameter_coef) {
    if (this->contains_parameter(parameter_id)) {
        throw ParameterIdDuplicatedException(
            "Error: "
            "ParametericQuantumCircuit::add_parametric_multi_Pauli_rotation_"
            "gate_new_"
            "parameter("
            "std::vector<UINT>, std::vector<UINT>, const ParameterId&, double, "
            "double): "
            "parameter_id \"" +
            parameter_id + "\" already exists");
    }
    this->create_parameter(parameter_id, value);
    this->add_parametric_gate(gate::ParametricPauliRotation(
        target, pauli_id, parameter_id, parameter_coef));
}

void ParametricQuantumCircuit::update_quantum_state(QuantumStateBase* state) {
    if (state->qubit_count != this->qubit_count) {
        throw InvalidQubitCountException(
            "Error: "
            "QuantumCircuit::update_quantum_state(QuantumStateBase) : "
            "invalid qubit count");
    }

    for (const auto& gate : this->_gate_list) {
        if (gate->is_parametric()) {
            dynamic_cast<QuantumGate_SingleParameter*>(gate)
                ->update_quantum_state(state, _parameter_set);
        } else {
            gate->update_quantum_state(state);
        }
    }
}

void ParametricQuantumCircuit::update_quantum_state(
    QuantumStateBase* state, UINT start, UINT end) {
    if (state->qubit_count != this->qubit_count) {
        throw InvalidQubitCountException(
            "Error: "
            "QuantumCircuit::update_quantum_state(QuantumStateBase,UINT,"
            "UINT) : invalid qubit count");
    }
    if (start > end) {
        throw GateIndexOutOfRangeException(
            "Error: "
            "QuantumCircuit::update_quantum_state(QuantumStateBase,UINT,"
            "UINT) : start must be smaller than or equal to end");
    }
    if (end > this->_gate_list.size()) {
        throw GateIndexOutOfRangeException(
            "Error: "
            "QuantumCircuit::update_quantum_state(QuantumStateBase,UINT,"
            "UINT) : end must be smaller than or equal to gate_count");
    }
    for (UINT cursor = start; cursor < end; ++cursor) {
        auto gate = this->_gate_list[cursor];
        if (gate->is_parametric()) {
            dynamic_cast<QuantumGate_SingleParameter*>(gate)
                ->update_quantum_state(state, _parameter_set);
        } else {
            gate->update_quantum_state(state);
        }
    }
}

std::vector<double> ParametricQuantumCircuit::backprop_inner_product(
    QuantumState* bistate) {
    // circuitを実行した状態とbistateの、inner_productを取った結果を「値」として、それを逆誤差伝搬します
    // bistateはノルムが1のやつでなくてもよい
    int n = this->qubit_count;
    QuantumState* state = new QuantumState(n);
    //これは、ゲートを前から適用したときの状態を示す
    state->set_zero_state();
    this->update_quantum_state(state);  //一度最後までする

    int num_gates = this->gate_list.size();
    std::vector<int> inverse_parametric_gate_position(num_gates, -1);
    for (UINT i = 0; i < this->get_parameter_count(); i++) {
        inverse_parametric_gate_position[this->_parametric_gate_position[i]] =
            i;
    }
    std::vector<double> ans(this->get_parameter_count());

    /*
    現在、2番のゲートを見ているとする
    ゲート 0 1 2 3 4 5
         state | bistate
    前から2番までのゲートを適用した状態がstate
    最後の微分値から逆算して3番までgateの逆行列を掛けたのがbistate

    1番まで掛けて、YのΘ微分した行列を掛けたやつと、bistateの内積の実数部分をとれば答えが出ることが知られている(知られてないかも)

    ParametricR? の微分値を計算した行列は、Θに180°を足した行列/2 と等しい

    だから、2番まで掛けて、 R?(π) を掛けたやつと、bistateの内積を取る

    さらに、見るゲートは逆順である。
    だから、最初にstateを最後までやって、ゲートを進めるたびにstateに逆行列を掛けている
    さらに、bistateが複素共役になっていることを忘れると、bistateに転置行列を掛ける必要がある。
    しかしこのプログラムではbistateはずっと複素共役なので、転置して共役な行列を掛ける必要がある。
    ユニタリ性より、転置して共役な行列 = 逆行列
    なので、両社にadjoint_gateを掛けている
    */
    QuantumState* Astate = new QuantumState(n);  //一時的なやつ
    for (int i = num_gates - 1; i >= 0; i--) {
        QuantumGateBase* gate_now = this->gate_list[i];  // sono gate
        if (inverse_parametric_gate_position[i] != -1) {
            Astate->load(state);
            QuantumGateBase* RcPI;
            if (gate_now->get_name() == "ParametricRX") {
                RcPI = gate::RX(gate_now->get_target_index_list()[0], M_PI);
            } else if (gate_now->get_name() == "ParametricRY") {
                RcPI = gate::RY(gate_now->get_target_index_list()[0], M_PI);
            } else if (gate_now->get_name() == "ParametricRZ") {
                RcPI = gate::RZ(gate_now->get_target_index_list()[0],
                    M_PI);  // 本当はここで2で割りたいけど、行列を割るのは実装が面倒
            } else if (gate_now->get_name() == "ParametricPauliRotation") {
                ClsParametricPauliRotationGate* pauli_gate_now =
                    (ClsParametricPauliRotationGate*)gate_now;
                RcPI =
                    gate::PauliRotation(pauli_gate_now->get_target_index_list(),
                        pauli_gate_now->get_pauli()->get_pauli_id_list(), M_PI);
            } else {
                std::stringstream error_message_stream;
                error_message_stream
                    << "Error: " << gate_now->get_name()
                    << " does not support backprop in parametric";
                throw NotImplementedException(error_message_stream.str());
            }
            RcPI->update_quantum_state(Astate);
            ans[inverse_parametric_gate_position[i]] =
                state::inner_product(bistate, Astate).real() /
                2.0;  //だからここで2で割る
            delete RcPI;
        }
        auto Agate = gate::get_adjoint_gate(gate_now);
        Agate->update_quantum_state(bistate);
        Agate->update_quantum_state(state);
        delete Agate;
    }
    delete Astate;
    delete state;
    return ans;
}  // CPP

std::vector<double> ParametricQuantumCircuit::backprop(
    GeneralQuantumOperator* obs) {
    //オブザーバブルから、最終段階での微分値を求めて、backprop_from_stateに流す関数
    //上側から来た変動量 * 下側の対応する微分値 =
    //最終的な変動量になるようにする。

    int n = this->qubit_count;
    QuantumState* state = new QuantumState(n);
    state->set_zero_state();
    this->update_quantum_state(state);  //一度最後までする
    QuantumState* bistate = new QuantumState(n);
    QuantumState* Astate = new QuantumState(n);  //一時的なやつ

    obs->apply_to_state(Astate, *state, bistate);
    bistate->multiply_coef(2);
    /*一度stateを最後まで求めてから、さらにapply_to_state している。
    なぜなら、量子のオブザーバブルは普通の機械学習と違って、二乗した値の絶対値が観測値になる。
    二乗の絶対値を微分したやつと、値の複素共役*2は等しい

    オブザーバブルよくわからないけど、テストしたらできてた
    */

    //ニューラルネットワークのbackpropにおける、後ろからの微分値的な役目を果たす
    auto ans = backprop_inner_product(bistate);
    delete bistate;
    delete state;
    delete Astate;
    return ans;

}  // CPP
