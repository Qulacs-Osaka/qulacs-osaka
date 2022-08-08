#pragma once

#include <cppsim/exception.hpp>
#include <cppsim/gate.hpp>
#include <cppsim/pauli_operator.hpp>
#include <cppsim/state.hpp>
#include <cppsim/utility.hpp>
#include <csim/update_ops.hpp>
#include <csim/update_ops_dm.hpp>

#ifdef _USE_GPU
#include <gpusim/update_ops_cuda.h>
#endif

#include "parameter.hpp"

class QuantumGate_SingleParameter : public QuantumGateBase {
protected:
    ParameterKey _parameter_key;
    double _parameter_coef;

public:
    QuantumGate_SingleParameter(
        const ParameterKey& parameter_key, double parameter_coef = 1.)
        : _parameter_key(parameter_key), _parameter_coef(parameter_coef) {
        _gate_property |= FLAG_PARAMETRIC;
    }
    ParameterKey get_parameter_key() { return _parameter_key; }
    ParameterId get_parameter_id() {
        return parameter::get_parameter_id(_parameter_key);
    }
    double get_parameter_coef() { return _parameter_coef; }
    void set_parameter_value(
        double parameter, ParameterSet& parameter_set) const {
        parameter::set_parameter_value(
            _parameter_key, parameter / _parameter_coef, parameter_set);
    }
    void set_parameter_value(double parameter) const {
        parameter::set_parameter_value(
            _parameter_key, parameter / _parameter_coef);
    }
    double get_parameter_value(const ParameterSet& parameter_set = {}) const {
        return parameter::get_parameter_value(_parameter_key, parameter_set) *
               _parameter_coef;
    }
    virtual QuantumGate_SingleParameter* copy() const override = 0;
};

class QuantumGate_SingleParameterOneQubitRotation
    : public QuantumGate_SingleParameter {
protected:
    typedef void(T_UPDATE_FUNC)(UINT, double, CTYPE*, ITYPE);
    typedef void(T_GPU_UPDATE_FUNC)(UINT, double, void*, ITYPE, void*, UINT);
    T_UPDATE_FUNC* _update_func = NULL;
    T_UPDATE_FUNC* _update_func_dm = NULL;
    T_GPU_UPDATE_FUNC* _update_func_gpu = NULL;

public:
    QuantumGate_SingleParameterOneQubitRotation(
        const ParameterKey& parameter_key, double parameter_coef = 1.)
        : QuantumGate_SingleParameter(parameter_key, parameter_coef) {}
    virtual void update_quantum_state(
        QuantumStateBase* state, const ParameterSet& parameter_set) {
        double angle = this->get_parameter_value(parameter_set);
        if (state->is_state_vector()) {
#ifdef _USE_GPU
            if (state->get_device_name() == "gpu") {
                if (_update_func_gpu == NULL) {
                    throw UndefinedUpdateFuncException(
                        "Error: "
                        "QuantumGate_SingleParameterOneQubitRotation::update_"
                        "quantum_state(QuantumStateBase) : update function is "
                        "undefined");
                }
                _update_func_gpu(this->_target_qubit_list[0].index(), angle,
                    state->data(), state->dim, state->get_cuda_stream(),
                    state->device_number);
                return;
            }
#endif
            if (_update_func == NULL) {
                throw UndefinedUpdateFuncException(
                    "Error: "
                    "QuantumGate_SingleParameterOneQubitRotation::update_"
                    "quantum_state(QuantumStateBase) : update function is "
                    "undefined");
            }
            _update_func(this->_target_qubit_list[0].index(), angle,
                state->data_c(), state->dim);
        } else {
            if (_update_func_dm == NULL) {
                throw UndefinedUpdateFuncException(
                    "Error: "
                    "QuantumGate_SingleParameterOneQubitRotation::update_"
                    "quantum_state(QuantumStateBase) : update function is "
                    "undefined");
            }
            _update_func_dm(this->_target_qubit_list[0].index(), angle,
                state->data_c(), state->dim);
        }
    }
    virtual void update_quantum_state(QuantumStateBase* state) {
        update_quantum_state(state, {});
    }
};

class ClsParametricRXGate : public QuantumGate_SingleParameterOneQubitRotation {
public:
    ClsParametricRXGate(UINT target_qubit_index,
        const ParameterKey& parameter_key, double parameter_coef = 1.)
        : QuantumGate_SingleParameterOneQubitRotation(
              parameter_key, parameter_coef) {
        this->_name = "ParametricRX";
        this->_update_func = RX_gate;
        this->_update_func_dm = dm_RX_gate;
#ifdef _USE_GPU
        this->_update_func_gpu = RX_gate_host;
#endif
        this->_target_qubit_list.push_back(
            TargetQubitInfo(target_qubit_index, FLAG_X_COMMUTE));
    }
    virtual void set_matrix(
        ComplexMatrix& matrix, const ParameterSet& parameter_set) const {
        double angle = this->get_parameter_value(parameter_set);
        matrix = ComplexMatrix::Zero(2, 2);
        matrix << cos(angle / 2), sin(angle / 2) * 1.i, sin(angle / 2) * 1.i,
            cos(angle / 2);
    }
    virtual void set_matrix(ComplexMatrix& matrix) const {
        set_matrix(matrix, {});
    }
    virtual QuantumGate_SingleParameter* copy() const override {
        auto* copied_gate = new ClsParametricRXGate(*this);
        ParameterKey key = this->_parameter_key;
        if (parameter::get_parameter_type(key) == "global") {
            copied_gate->_parameter_key = parameter::create_parameter(
                parameter::get_parameter_value(key, {}));
        }
        return copied_gate;
    };
};

class ClsParametricRYGate : public QuantumGate_SingleParameterOneQubitRotation {
public:
    ClsParametricRYGate(UINT target_qubit_index,
        const ParameterKey& parameter_key, double parameter_coef = 1.)
        : QuantumGate_SingleParameterOneQubitRotation(
              parameter_key, parameter_coef) {
        this->_name = "ParametricRY";
        this->_update_func = RY_gate;
        this->_update_func_dm = dm_RY_gate;
#ifdef _USE_GPU
        this->_update_func_gpu = RY_gate_host;
#endif
        this->_target_qubit_list.push_back(
            TargetQubitInfo(target_qubit_index, FLAG_Y_COMMUTE));
    }
    virtual void set_matrix(
        ComplexMatrix& matrix, const ParameterSet& parameter_set) const {
        double angle = this->get_parameter_value(parameter_set);
        matrix = ComplexMatrix::Zero(2, 2);
        matrix << cos(angle / 2), sin(angle / 2), -sin(angle / 2),
            cos(angle / 2);
    }
    virtual void set_matrix(ComplexMatrix& matrix) const {
        set_matrix(matrix, {});
    }
    virtual QuantumGate_SingleParameter* copy() const override {
        auto* copied_gate = new ClsParametricRYGate(*this);
        ParameterKey key = this->_parameter_key;
        if (parameter::get_parameter_type(key) == "global") {
            copied_gate->_parameter_key = parameter::create_parameter(
                parameter::get_parameter_value(key, {}));
        }
        return copied_gate;
    };
};

class ClsParametricRZGate : public QuantumGate_SingleParameterOneQubitRotation {
public:
    ClsParametricRZGate(UINT target_qubit_index,
        const ParameterKey& parameter_key, double parameter_coef = 1.)
        : QuantumGate_SingleParameterOneQubitRotation(
              parameter_key, parameter_coef) {
        this->_name = "ParametricRZ";
        this->_update_func = RZ_gate;
        this->_update_func_dm = dm_RZ_gate;
#ifdef _USE_GPU
        this->_update_func_gpu = RZ_gate_host;
#endif
        this->_target_qubit_list.push_back(
            TargetQubitInfo(target_qubit_index, FLAG_Z_COMMUTE));
    }
    virtual void set_matrix(
        ComplexMatrix& matrix, const ParameterSet& parameter_set) const {
        double angle = this->get_parameter_value(parameter_set);
        matrix = ComplexMatrix::Zero(2, 2);
        matrix << cos(angle / 2) + 1.i * sin(angle / 2), 0, 0,
            cos(angle / 2) - 1.i * sin(angle / 2);
    }
    virtual void set_matrix(ComplexMatrix& matrix) const {
        set_matrix(matrix, {});
    }
    virtual QuantumGate_SingleParameter* copy() const override {
        auto* copied_gate = new ClsParametricRZGate(*this);
        ParameterKey key = this->_parameter_key;
        if (parameter::get_parameter_type(key) == "global") {
            copied_gate->_parameter_key = parameter::create_parameter(
                parameter::get_parameter_value(key, {}));
        }
        return copied_gate;
    };
};

class ClsParametricPauliRotationGate : public QuantumGate_SingleParameter {
protected:
    PauliOperator* _pauli;

public:
    ClsParametricPauliRotationGate(PauliOperator* pauli,
        const ParameterKey& parameter_key, double parameter_coef = 1.)
        : QuantumGate_SingleParameter(parameter_key, parameter_coef) {
        _pauli = pauli;
        this->_name = "ParametricPauliRotation";
        auto target_index_list = _pauli->get_index_list();
        auto pauli_id_list = _pauli->get_pauli_id_list();
        for (UINT index = 0; index < target_index_list.size(); ++index) {
            UINT commutation_relation = 0;
            if (pauli_id_list[index] == 1)
                commutation_relation = FLAG_X_COMMUTE;
            else if (pauli_id_list[index] == 2)
                commutation_relation = FLAG_Y_COMMUTE;
            else if (pauli_id_list[index] == 3)
                commutation_relation = FLAG_Z_COMMUTE;
            this->_target_qubit_list.push_back(TargetQubitInfo(
                target_index_list[index], commutation_relation));
        }
    }
    virtual ~ClsParametricPauliRotationGate() { delete _pauli; }
    virtual void update_quantum_state(
        QuantumStateBase* state, const ParameterSet& parameter_set) {
        auto target_index_list = _pauli->get_index_list();
        auto pauli_id_list = _pauli->get_pauli_id_list();
        double angle = this->get_parameter_value(parameter_set);
        if (state->is_state_vector()) {
#ifdef _USE_GPU
            if (state->get_device_name() == "gpu") {
                multi_qubit_Pauli_rotation_gate_partial_list_host(
                    target_index_list.data(), pauli_id_list.data(),
                    (UINT)target_index_list.size(), angle, state->data(),
                    state->dim, state->get_cuda_stream(), state->device_number);
            } else {
                multi_qubit_Pauli_rotation_gate_partial_list(
                    target_index_list.data(), pauli_id_list.data(),
                    (UINT)target_index_list.size(), angle, state->data_c(),
                    state->dim);
            }
#else
            multi_qubit_Pauli_rotation_gate_partial_list(
                target_index_list.data(), pauli_id_list.data(),
                (UINT)target_index_list.size(), angle, state->data_c(),
                state->dim);
#endif
        } else {
            dm_multi_qubit_Pauli_rotation_gate_partial_list(
                target_index_list.data(), pauli_id_list.data(),
                (UINT)target_index_list.size(), angle, state->data_c(),
                state->dim);
        }
    };
    virtual void update_quantum_state(QuantumStateBase* state) {
        update_quantum_state(state, {});
    }
    virtual QuantumGate_SingleParameter* copy() const override {
        auto* copied_gate =
            new ClsParametricPauliRotationGate(_pauli->copy(), _parameter_key);
        ParameterKey key = this->_parameter_key;
        if (parameter::get_parameter_type(key) == "global") {
            copied_gate->_parameter_key = parameter::create_parameter(
                parameter::get_parameter_value(key, {}));
        }
        return copied_gate;
    };
    virtual void set_matrix(
        ComplexMatrix& matrix, const ParameterSet& parameter_set) const {
        double angle = this->get_parameter_value(parameter_set);
        get_Pauli_matrix(matrix, _pauli->get_pauli_id_list());
        matrix = cos(angle / 2) *
                     ComplexMatrix::Identity(matrix.rows(), matrix.cols()) +
                 1.i * sin(angle / 2) * matrix;
    };
    virtual void set_matrix(ComplexMatrix& matrix) const {
        set_matrix(matrix, {});
    }
    virtual PauliOperator* get_pauli() const { return _pauli; };
};
