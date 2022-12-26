#pragma once

#include <cppsim/simulator.hpp>

class DllExport ParametricQuantumCircuitSimulator
    : public QuantumCircuitSimulator {
private:
    ParametricQuantumCircuit* _parametric_circuit;

public:
    ParametricQuantumCircuitSimulator(
        ParametricQuantumCircuit* circuit, QuantumStateBase* state = NULL);
    bool is_old_style() const;
    bool is_new_style() const;
    double get_parameter(UINT index) const;
    double get_parameter_new_style(const ParameterId& parameter_id) const;
    void add_parameter_value(UINT index, double value);
    void add_parameter_value_new_style(
        const ParameterId& parameter_id, double value);
    void set_parameter_value(UINT index, double value);
    void set_parameter_value_new_style(
        const ParameterId& parameter_id, double value);
    UINT get_parametric_gate_count();
    UINT get_parameter_id_count();
    UINT get_parametric_gate_position(UINT index);
};
