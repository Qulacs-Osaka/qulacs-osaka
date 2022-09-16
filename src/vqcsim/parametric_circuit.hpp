#pragma once

#include <cppsim/circuit.hpp>
#include <cppsim/observable.hpp>
#include <cppsim/state.hpp>

#include "parametric_gate.hpp"

class DllExport ParametricQuantumCircuit : public QuantumCircuit {
private:
    std::vector<QuantumGate_SingleParameter*> _parametric_gate_list;
    std::vector<UINT> _parametric_gate_position;
    std::vector<double> _parameter_list;
    bool _specified_old;
    bool _specified_new;

public:
    ParametricQuantumCircuit(UINT qubit_count, std::string style = "undefined");
    ParametricQuantumCircuit(UINT qubit_count,
        const std::vector<double>&
            parameter_list);            // automatically specify new-style
    virtual bool is_old_style() const;  // both-style
    virtual bool is_new_style() const;  // both-style

    virtual ParametricQuantumCircuit* copy() const;  // both-style

    virtual ParameterId create_parameter(
        double initial_parameter);  // new-style
    virtual bool contains_parameter(
        const ParameterId& parameter_id) const;                // new-style
    virtual UINT get_parameter_count() const;                  // old-style
    virtual UINT get_parametric_gate_count() const;            // new-style
    virtual UINT get_parameter_id_count() const;               // new-style
    virtual double get_parameter(UINT parameter_index) const;  // old-style
    virtual double get_parameter_new_style(const ParameterId& parameter_id)
        const;  // new-style (いつかget_parameter()をこっちにしたい)
    virtual void set_parameter(
        UINT parameter_index, double value);  // old-style
    virtual void set_parameter_new_style(const ParameterId& parameter_id,
        double value);  // new-style (いつかset_parameter()をこっちにしたい)
    virtual double get_angle(UINT gate_index) const;         // both-style
    virtual std::vector<double> get_parameter_list() const;  // new-style
    virtual void set_parameter_list(
        const std::vector<double>& parameter_list);  // new-style

    virtual UINT get_parametric_gate_position(
        UINT parameter_index) const;                        // both-style
    virtual void add_gate(QuantumGateBase* gate) override;  // both-style
    virtual void add_gate(
        QuantumGateBase* gate, UINT index) override;  // both-style
    virtual void add_gate_copy(
        const QuantumGateBase* gate) override;  // both-style
    virtual void add_gate_copy(
        const QuantumGateBase* gate, UINT index) override;  // both-style
    virtual void remove_gate(UINT index) override;          // both-style
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
     * @param[in] share_parameter_id
     * マージ元と共有するパラメータIDの変換表(マージする側→マージされる側)
     */
    virtual void merge_circuit(const ParametricQuantumCircuit* circuit,
        const std::map<ParameterId, ParameterId>& share_parameter_id =
            {});  // both-style(share_parameter_idを渡すのはnewのみ)
    virtual std::string to_string() const;  // both-style
    friend DllExport std::ostream& operator<<(std::ostream& os,
        const ParametricQuantumCircuit& circuit);  // both-style
    friend DllExport std::ostream& operator<<(std::ostream& os,
        const ParametricQuantumCircuit* circuit);  // both-style

    virtual ParameterId add_parametric_gate(QuantumGate_SingleParameter*
            gate);  // both-style
                    // (new-styleで新規パラメータにしたいならcreate_parameter()を先に呼ぶ)
    virtual ParameterId add_parametric_gate(QuantumGate_SingleParameter* gate,
        UINT
            index);  // both-style
                     // (new-styleで新規パラメータにしたいならcreate_parameter()を先に呼ぶ)
    virtual ParameterId add_parametric_gate_copy(
        const QuantumGate_SingleParameter*
            gate);  // both-style
                    // (newで新規パラメータにしたいならcreate_parameter()を先に呼ぶ)
    virtual ParameterId add_parametric_gate_copy(
        const QuantumGate_SingleParameter* gate,
        UINT
            index);  // both-style
                     // (new-styleで新規パラメータにしたいならcreate_parameter()を先に呼ぶ)

    virtual ParameterId add_parametric_RX_gate(
        UINT target_index, double initial_angle);  // old-style
    virtual ParameterId add_parametric_RX_gate_existing_parameter(
        UINT target_index, const ParameterId& parameter_id,
        double parameter_coef = 1.);  // new-style
    virtual ParameterId add_parametric_RX_gate_new_parameter(UINT target_index,
        double value,
        double parameter_coef = 1.);  // new-style
    virtual ParameterId add_parametric_RY_gate(
        UINT target_index, double initial_angle);  // old-style
    virtual ParameterId add_parametric_RY_gate_existing_parameter(
        UINT target_index, const ParameterId& parameter_id,
        double parameter_coef = 1.);  // new-style
    virtual ParameterId add_parametric_RY_gate_new_parameter(UINT target_index,
        double value,
        double parameter_coef = 1.);  // new-style
    virtual ParameterId add_parametric_RZ_gate(
        UINT target_index, double initial_angle);  // old-style
    virtual ParameterId add_parametric_RZ_gate_existing_parameter(
        UINT target_index, const ParameterId& parameter_id,
        double parameter_coef = 1.);  // new-style
    virtual ParameterId add_parametric_RZ_gate_new_parameter(UINT target_index,
        double value,
        double parameter_coef = 1.);  // new-style
    virtual ParameterId add_parametric_multi_Pauli_rotation_gate(
        std::vector<UINT> target, std::vector<UINT> pauli_id,
        double initial_angle);  // old-style
    virtual ParameterId
    add_parametric_multi_Pauli_rotation_gate_existing_parameter(
        std::vector<UINT> target, std::vector<UINT> pauli_id,
        const ParameterId& parameter_id,
        double parameter_coef = 1.);  // new-style
    virtual ParameterId add_parametric_multi_Pauli_rotation_gate_new_parameter(
        std::vector<UINT> target, std::vector<UINT> pauli_id, double value,
        double parameter_coef = 1.);  // new-style

    virtual void update_quantum_state(QuantumStateBase* state);  // both-style
    virtual void update_quantum_state(QuantumStateBase* state, UINT start_index,
        UINT end_index);  // both-style

    virtual std::vector<double> backprop(
        GeneralQuantumOperator* obs);  // old-style
    virtual std::vector<double> backprop_inner_product(
        QuantumState* bistate);  // old-style
};
