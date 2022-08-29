#pragma once

#include <cppsim/circuit.hpp>
#include <cppsim/observable.hpp>
#include <cppsim/state.hpp>

#include "parametric_gate.hpp"

class DllExport ParametricQuantumCircuit : public QuantumCircuit {
private:
    std::vector<QuantumGate_SingleParameter*> _parametric_gate_list;
    std::vector<UINT> _parametric_gate_position;
    ParameterSet _parameter_set;
    bool _specified_old;
    bool _specified_new;

public:
    ParametricQuantumCircuit(UINT qubit_count, std::string style = "undefined");
    ParametricQuantumCircuit(UINT qubit_count,
        const ParameterSet& parameter_set, std::string style = "undefined");
    virtual bool is_old_style() const;
    virtual bool is_new_style() const;

    ParametricQuantumCircuit* copy() const;

    virtual void create_parameter(
        const ParameterId& parameter_id, double initial_parameter);
    virtual void remove_parameter(const ParameterId& parameter_id);
    bool contains_parameter(const ParameterId& parameter_id) const;
    virtual UINT get_parameter_count() const;
    virtual UINT get_parametric_gate_count() const;
    virtual UINT get_parameter_id_count() const;
    virtual double get_parameter(UINT index) const;
    virtual double get_parameter(const ParameterId& parameter_id) const;
    virtual double get_angle(UINT index) const;
    virtual void set_parameter(UINT index, double value);
    virtual void set_parameter(const ParameterId& parameter_id, double value);
    virtual ParameterSet get_parameter_set() const;
    virtual void set_parameter_set(const ParameterSet& parameter_set);

    virtual UINT get_parametric_gate_position(UINT index) const;
    virtual void add_gate(QuantumGateBase* gate) override;
    virtual void add_gate(QuantumGateBase* gate, UINT index) override;
    virtual void add_gate_copy(const QuantumGateBase* gate) override;
    virtual void add_gate_copy(
        const QuantumGateBase* gate, UINT index) override;
    virtual void remove_gate(UINT index) override;
    /**
     *  \~japanese-en 量子回路をマージする。
     *
     * 引数で与えた量子回路のゲートを後ろに追加していく。
     * マージされた側の量子回路に変更を加えてもマージした側の量子回路には変更は加わらないことに注意する。
     * パラメータゲートに対応するため、ParametricQuantumCircuit にも
     * merge_circuit() を追加する circuit1.add_circuit(circuit2)
     * circuit2.add_gate(gate) # これをしても、circuit1にgateは追加されない
     *
     * @param[in] circuit マージする量子回路
     */
    virtual void merge_circuit(const ParametricQuantumCircuit* circuit);
    virtual std::string to_string() const;
    friend DllExport std::ostream& operator<<(
        std::ostream& os, const ParametricQuantumCircuit&);
    friend DllExport std::ostream& operator<<(
        std::ostream& os, const ParametricQuantumCircuit* circuit);

    virtual void add_parametric_gate(QuantumGate_SingleParameter* gate);
    virtual void add_parametric_gate(
        QuantumGate_SingleParameter* gate, UINT index);
    virtual void add_parametric_gate_copy(QuantumGate_SingleParameter* gate);
    virtual void add_parametric_gate_copy(
        QuantumGate_SingleParameter* gate, UINT index);

    virtual void add_parametric_RX_gate(
        UINT target_index, double initial_angle);
    virtual void add_parametric_RX_gate(UINT target_index,
        const ParameterId& parameter_id, double parameter_coef = 1.);
    virtual void add_parametric_RX_gate_new_parameter(UINT target_index,
        const ParameterId& parameter_id, double value,
        double parameter_coef = 1.);
    virtual void add_parametric_RY_gate(
        UINT target_index, double initial_angle);
    virtual void add_parametric_RY_gate(UINT target_index,
        const ParameterId& parameter_id, double parameter_coef = 1.);
    virtual void add_parametric_RY_gate_new_parameter(UINT target_index,
        const ParameterId& parameter_id, double value,
        double parameter_coef = 1.);
    virtual void add_parametric_RZ_gate(
        UINT target_index, double initial_angle);
    virtual void add_parametric_RZ_gate(UINT target_index,
        const ParameterId& parameter_id, double parameter_coef = 1.);
    virtual void add_parametric_RZ_gate_new_parameter(UINT target_index,
        const ParameterId& parameter_id, double value,
        double parameter_coef = 1.);
    virtual void add_parametric_multi_Pauli_rotation_gate(
        std::vector<UINT> target, std::vector<UINT> pauli_id,
        double initial_angle);
    virtual void add_parametric_multi_Pauli_rotation_gate(
        std::vector<UINT> target, std::vector<UINT> pauli_id,
        const ParameterId& parameter_id, double parameter_coef = 1.);
    virtual void add_parametric_multi_Pauli_rotation_gate_new_parameter(
        std::vector<UINT> target, std::vector<UINT> pauli_id,
        const ParameterId& parameter_id, double value,
        double parameter_coef = 1.);

    virtual void update_quantum_state(QuantumStateBase* state);
    virtual void update_quantum_state(
        QuantumStateBase* state, UINT start_index, UINT end_index);

    virtual std::vector<double> backprop(GeneralQuantumOperator* obs);
    virtual std::vector<double> backprop_inner_product(QuantumState* bistate);
};
