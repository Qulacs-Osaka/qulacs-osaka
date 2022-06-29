
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

#include "utility.hpp"

class SingleParameter {
protected:
    double _value;
    UINT _id;
    std::vector<UINT> _gate_indices;

public:
    SingleParameter(double value, UINT id) : _value(value), _id(id) {}
    void set_parameter_value(double value) { _value = value; }
    double get_parameter_value() const { return _value; }
    UINT get_parameter_id() const { return _id; }
    std::vector<UINT> get_gate_indices() const { return _gate_indices; }
    void push_gate_index(UINT index) { _gate_indices.push_back(index); }
    void pop_gate_index() {
        if (_gate_indices.empty()) {
            throw GateIndexOutOfRangeException(
                "SingleParameter::pop_gate_index(): gate_indices is empty");
        }
        _gate_indices.pop_back();
    }
};

class QuantumGate_SingleParameter : public QuantumGateBase {
protected:
    SingleParameter* _parameter;

public:
    QuantumGate_SingleParameter(SingleParameter* parameter)
        : _parameter(parameter) {
        _gate_property |= FLAG_PARAMETRIC;
    }
    virtual QuantumGate_SingleParameter* copy() const override = 0;
};

class QuantumGate_SingleParameterOneQubitRotation
    : public QuantumGate_SingleParameter {
protected:
    typedef double(T_ANGLE_FUNC)(double);
    typedef void(T_UPDATE_FUNC)(UINT, double, CTYPE*, ITYPE);
    typedef void(T_GPU_UPDATE_FUNC)(UINT, double, void*, ITYPE, void*, UINT);
    T_ANGLE_FUNC* _angle_func;
    T_UPDATE_FUNC* _update_func = NULL;
    T_UPDATE_FUNC* _update_func_dm = NULL;
    T_GPU_UPDATE_FUNC* _update_func_gpu = NULL;

    QuantumGate_SingleParameterOneQubitRotation(
        SingleParameter* parameter, T_ANGLE_FUNC* angle_func = identity_map)
        : QuantumGate_SingleParameter(parameter) {
        _angle_func = angle_func;
    }

public:
    virtual void update_quantum_state(QuantumStateBase* state) override {
        double angle = _angle_func(_parameter->get_parameter_value());
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
};

class ClsParametricRXGate : public QuantumGate_SingleParameterOneQubitRotation {
public:
    ClsParametricRXGate(UINT target_qubit_index, SingleParameter* parameter,
        typeof(_angle_func) angle_func = identity_map)
        : QuantumGate_SingleParameterOneQubitRotation(parameter, angle_func) {
        this->_name = "ParametricRX";
        this->_update_func = RX_gate;
        this->_update_func_dm = dm_RX_gate;
#ifdef _USE_GPU
        this->_update_func_gpu = RX_gate_host;
#endif
        this->_target_qubit_list.push_back(
            TargetQubitInfo(target_qubit_index, FLAG_X_COMMUTE));
    }
    virtual void set_matrix(ComplexMatrix& matrix) const override {
        double angle = _angle_func(_parameter->get_parameter_value());
        matrix = ComplexMatrix::Zero(2, 2);
        matrix << cos(angle / 2), sin(angle / 2) * 1.i, sin(angle / 2) * 1.i,
            cos(angle / 2);
    }
    virtual QuantumGate_SingleParameter* copy() const override {
        throw NotImplementedException(
            "Error: ClsParametricRXGate::copy(): ParametricGate "
            "cannot be copied");
    };
};

class ClsParametricRYGate : public QuantumGate_SingleParameterOneQubitRotation {
public:
    ClsParametricRYGate(UINT target_qubit_index, SingleParameter* parameter,
        typeof(_angle_func) angle_func = identity_map)
        : QuantumGate_SingleParameterOneQubitRotation(parameter, angle_func) {
        this->_name = "ParametricRY";
        this->_update_func = RY_gate;
        this->_update_func_dm = dm_RY_gate;
#ifdef _USE_GPU
        this->_update_func_gpu = RY_gate_host;
#endif
        this->_target_qubit_list.push_back(
            TargetQubitInfo(target_qubit_index, FLAG_Y_COMMUTE));
    }
    virtual void set_matrix(ComplexMatrix& matrix) const override {
        double angle = _angle_func(_parameter->get_parameter_value());
        matrix = ComplexMatrix::Zero(2, 2);
        matrix << cos(angle / 2), sin(angle / 2), -sin(angle / 2),
            cos(angle / 2);
    }
    virtual QuantumGate_SingleParameter* copy() const override {
        throw NotImplementedException(
            "Error: ClsParametricRYGate::copy(): ParametricGate "
            "cannot be copied");
    };
};

class ClsParametricRZGate : public QuantumGate_SingleParameterOneQubitRotation {
public:
    ClsParametricRZGate(UINT target_qubit_index, SingleParameter* parameter,
        typeof(_angle_func) angle_func = identity_map)
        : QuantumGate_SingleParameterOneQubitRotation(parameter, angle_func) {
        this->_name = "ParametricRZ";
        this->_update_func = RZ_gate;
        this->_update_func_dm = dm_RZ_gate;
#ifdef _USE_GPU
        this->_update_func_gpu = RZ_gate_host;
#endif
        this->_target_qubit_list.push_back(
            TargetQubitInfo(target_qubit_index, FLAG_Z_COMMUTE));
    }
    virtual void set_matrix(ComplexMatrix& matrix) const override {
        double angle = _angle_func(_parameter->get_parameter_value());
        matrix = ComplexMatrix::Zero(2, 2);
        matrix << cos(angle / 2) + 1.i * sin(angle / 2), 0, 0,
            cos(angle / 2) - 1.i * sin(angle / 2);
    }
    virtual QuantumGate_SingleParameter* copy() const override {
        throw NotImplementedException(
            "Error: ClsParametricRZGate::copy(): ParametricGate "
            "cannot be copied");
    };
};

class ClsParametricPauliRotationGate : public QuantumGate_SingleParameter {
protected:
    PauliOperator* _pauli;
    typedef double(T_ANGLE_FUNC)(double);
    T_ANGLE_FUNC* _angle_func;

public:
    ClsParametricPauliRotationGate(SingleParameter* parameter,
        PauliOperator* pauli, T_ANGLE_FUNC* angle_func = identity_map)
        : QuantumGate_SingleParameter(parameter) {
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
        _angle_func = angle_func;
    }
    virtual ~ClsParametricPauliRotationGate() { delete _pauli; }
    virtual void update_quantum_state(QuantumStateBase* state) override {
        auto target_index_list = _pauli->get_index_list();
        auto pauli_id_list = _pauli->get_pauli_id_list();
        double angle = _angle_func(_parameter->get_parameter_value());
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
    virtual QuantumGate_SingleParameter* copy() const override {
        throw NotImplementedException(
            "Error: ClsParametricPauliRotationGate::copy(): ParametricGate "
            "cannot be copied");
    };
    virtual void set_matrix(ComplexMatrix& matrix) const override {
        double angle = _angle_func(_parameter->get_parameter_value());
        get_Pauli_matrix(matrix, _pauli->get_pauli_id_list());
        matrix = cos(angle / 2) *
                     ComplexMatrix::Identity(matrix.rows(), matrix.cols()) +
                 1.i * sin(angle / 2) * matrix;
    };
    virtual PauliOperator* get_pauli() const { return _pauli; };
};
