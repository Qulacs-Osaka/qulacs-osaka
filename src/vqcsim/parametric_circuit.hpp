#pragma once

#include <cppsim/circuit.hpp>
#include <cppsim/observable.hpp>
#include <cppsim/state.hpp>

#include "parametric_gate.hpp"
#include "utility.hpp"

class DllExport ParametricQuantumCircuit : public QuantumCircuit {
private:
    std::vector<QuantumGate_SingleParameter*> _parametric_gate_list;
    std::vector<UINT> _parametric_gate_position;
    std::vector<SingleParameter> _parameter_list;

public:
    ParametricQuantumCircuit(UINT qubit_count);

    ParametricQuantumCircuit* copy() const;

    virtual UINT get_parameter_count() const;
    virtual double get_parameter(UINT index) const;
    virtual void set_parameter(UINT index, double value);

    virtual UINT get_parametric_gate_position(UINT index) const;
    virtual void add_gate(QuantumGateBase* gate) override;
    virtual void add_gate(QuantumGateBase* gate, UINT index) override;
    virtual void add_gate_copy(const QuantumGateBase* gate) override;
    virtual void add_gate_copy(
        const QuantumGateBase* gate, UINT index) override;
    virtual void remove_gate(UINT index) override;

    virtual std::string to_string() const;
    friend DllExport std::ostream& operator<<(
        std::ostream& os, const ParametricQuantumCircuit&);
    friend DllExport std::ostream& operator<<(
        std::ostream& os, const ParametricQuantumCircuit* circuit);

    virtual void add_parametric_RX_gate_new_parameter(UINT target_index,
        double initial_parameter, double (*angle_func)(double) = identity_map);
    virtual void add_parametric_RX_gate_share_parameter(UINT target_index,
        UINT parameter_id, double (*angle_func)(double) = identity_map);
    virtual void add_parametric_RY_gate_new_parameter(UINT target_index,
        double initial_parameter, double (*angle_func)(double) = identity_map);
    virtual void add_parametric_RY_gate_share_parameter(UINT target_index,
        UINT parameter_id, double (*angle_func)(double) = identity_map);
    virtual void add_parametric_RZ_gate_new_parameter(UINT target_index,
        double initial_parameter, double (*angle_func)(double) = identity_map);
    virtual void add_parametric_RZ_gate_share_parameter(UINT target_index,
        UINT parameter_id, double (*angle_func)(double) = identity_map);
    virtual void add_parametric_multi_Pauli_rotation_gate_new_parameter(
        std::vector<UINT> target, std::vector<UINT> pauli_id,
        double initial_parameter, double (*angle_func)(double) = identity_map);
    virtual void add_parametric_multi_Pauli_rotation_gate_share_parameter(
        std::vector<UINT> target, std::vector<UINT> pauli_id, UINT parameter_id,
        double (*angle_func)(double) = identity_map);
    virtual std::vector<double> backprop(GeneralQuantumOperator* obs);
    virtual std::vector<double> backprop_inner_product(QuantumState* bistate);
};
