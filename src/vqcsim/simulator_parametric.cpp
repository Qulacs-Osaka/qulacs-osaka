#include "parametric_circuit.hpp"
#include "parametric_simulator.hpp"

ParametricQuantumCircuitSimulator::ParametricQuantumCircuitSimulator(
    ParametricQuantumCircuit* circuit, QuantumStateBase* state)
    : QuantumCircuitSimulator(circuit, state), _parametric_circuit(circuit) {}

bool ParametricQuantumCircuitSimulator::is_old_style() const {
    return _parametric_circuit->is_old_style();
}
bool ParametricQuantumCircuitSimulator::is_new_style() const {
    return _parametric_circuit->is_new_style();
}
double ParametricQuantumCircuitSimulator::get_parameter(UINT index) const {
    return _parametric_circuit->get_parameter(index);
}
double ParametricQuantumCircuitSimulator::get_parameter_new_style(
    const ParameterId& parameter_id) const {
    return _parametric_circuit->get_parameter_new_style(parameter_id);
}
void ParametricQuantumCircuitSimulator::add_parameter_value(
    UINT index, double value) {
    _parametric_circuit->set_parameter(
        index, _parametric_circuit->get_parameter(index) + value);
}
void ParametricQuantumCircuitSimulator::add_parameter_value_new_style(
    const ParameterId& parameter_id, double value) {
    _parametric_circuit->set_parameter_new_style(parameter_id,
        _parametric_circuit->get_parameter_new_style(parameter_id) + value);
}
void ParametricQuantumCircuitSimulator::set_parameter_value(
    UINT index, double value) {
    _parametric_circuit->set_parameter(index, value);
}
void ParametricQuantumCircuitSimulator::set_parameter_value_new_style(
    const ParameterId& parameter_id, double value) {
    _parametric_circuit->set_parameter(parameter_id, value);
}
UINT ParametricQuantumCircuitSimulator::get_parametric_gate_count() {
    return _parametric_circuit->get_parameter_count();
}
UINT ParametricQuantumCircuitSimulator::get_parameter_id_count() {
    return _parametric_circuit->get_parameter_id_count();
}
UINT ParametricQuantumCircuitSimulator::get_parametric_gate_position(
    UINT index) {
    return _parametric_circuit->get_parametric_gate_position(index);
}
