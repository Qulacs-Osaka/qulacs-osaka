#pragma once

#include <cppsim/circuit.hpp>
#include <cppsim/observable.hpp>
#include <cppsim/state.hpp>

#include "parametric_gate.hpp"

class QuantumGate_SingleParameter;

class DllExport ParametricQuantumCircuit : public QuantumCircuit {
private:
    std::vector<QuantumGate_SingleParameter*> _parametric_gate_list;
    std::vector<UINT> _parametric_gate_position;
    ParameterSet _parameter_set;
    UINT _next_parameter_index = 0;

    virtual ParameterKey _generate_old_parameter_id(UINT index) const {
        return "index_" + std::to_string(index);
    }

public:
    ParametricQuantumCircuit(UINT qubit_count);

    ParametricQuantumCircuit* copy() const;

    virtual void add_parametric_gate(QuantumGate_SingleParameter* gate);
    virtual void add_parametric_gate(
        QuantumGate_SingleParameter* gate, UINT index);
    virtual void add_parametric_gate_copy(QuantumGate_SingleParameter* gate);
    virtual void add_parametric_gate_copy(
        QuantumGate_SingleParameter* gate, UINT index);
    virtual UINT get_parameter_count() const;
    virtual std::vector<ParameterKey> get_parameter_id_list() const;
    virtual void create_parameter(
        const ParameterKey& parameter_id, double initial_parameter);
    virtual void remove_parameter(const ParameterKey& parameter_id);
    bool contains_parameter(const ParameterKey& parameter_id) const;
    virtual double get_parameter(const ParameterKey& parameter_id) const;
    virtual void set_parameter(const ParameterKey& parameter_id, double value);

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
        ParameterKey& parameter_id, double initial_parameter);
    virtual void add_parametric_RX_gate_share_parameter(
        UINT target_index, ParameterKey& parameter_id);
    virtual void add_parametric_RY_gate_new_parameter(UINT target_index,
        ParameterKey& parameter_id, double initial_parameter);
    virtual void add_parametric_RY_gate_share_parameter(
        UINT target_index, ParameterKey& parameter_id);
    virtual void add_parametric_RZ_gate_new_parameter(UINT target_index,
        ParameterKey& parameter_id, double initial_parameter);
    virtual void add_parametric_RZ_gate_share_parameter(
        UINT target_index, ParameterKey& parameter_id);
    virtual void add_parametric_multi_Pauli_rotation_gate_new_parameter(
        std::vector<UINT> target, std::vector<UINT> pauli_id,
        ParameterKey& parameter_id, double initial_parameter);
    virtual void add_parametric_multi_Pauli_rotation_gate_share_parameter(
        std::vector<UINT> target, std::vector<UINT> pauli_id,
        ParameterKey& parameter_id);
    virtual std::vector<double> backprop(GeneralQuantumOperator* obs);
    virtual std::vector<double> backprop_inner_product(QuantumState* bistate);

    // 互換性のために残す旧メソッド(パラメータのIDとして整数値を利用)
    virtual double get_parameter(UINT index) const;
    virtual void set_parameter(UINT index, double value);
    virtual void add_parametric_RX_gate(
        UINT target_index, double initial_angle);
    virtual void add_parametric_RY_gate(
        UINT target_index, double initial_angle);
    virtual void add_parametric_RZ_gate(
        UINT target_index, double initial_angle);
    virtual void add_parametric_multi_Pauli_rotation_gate(
        std::vector<UINT> target, std::vector<UINT> pauli_id,
        double initial_angle);
};