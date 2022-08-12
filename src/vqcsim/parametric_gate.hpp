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
    ParameterType _parameter_type;
    ParameterId _parameter_id;
    double _parameter_coef;
    double _angle;

public:
    QuantumGate_SingleParameter(double angle) : _angle(angle) {
        _parameter_type = "gate";
        _parameter_coef = 1.;
        _gate_property |= FLAG_PARAMETRIC;
    }
    QuantumGate_SingleParameter(
        const ParameterId& parameter_id, double parameter_coef = 1.)
        : _parameter_id(parameter_id), _parameter_coef(parameter_coef) {
        _parameter_type = "id";
        _gate_property |= FLAG_PARAMETRIC;
    }
    ParameterType get_parameter_type() const { return _parameter_type; }
    ParameterId get_parameter_id() const {
        if (this->get_parameter_type() != "id") {
            throw NotImplementedException(
                "Error: "
                "QuantumGate_SingleParameter::get_parameter_id() const: "
                "this is a new-style function. You cannot use this function "
                "with an old-style parametric_gate which contains a parameter "
                "value");
        }
        return _parameter_id;
    }
    double get_parameter_coef() const { return _parameter_coef; }
    void set_parameter_value(double value) {
        if (this->get_parameter_type() != "gate") {
            throw NotImplementedException(
                "Error: "
                "QuantumGate_SingleParameter::set_parameter_value(double): "
                "this is an old-style function. You cannot use this function "
                "with a new-style parametric_gate which has string "
                "parameter_id.");
        }
        _angle = value;
    }
    double get_parameter_value() const {
        if (this->get_parameter_type() != "gate") {
            throw NotImplementedException(
                "Error: "
                "QuantumGate_SingleParameter::get_parameter_value() const: "
                "this is an old-style function. You cannot use this function "
                "with a new-style parametric_gate which has string "
                "parameter_id.");
        }
        return _angle;
    }
    double get_angle(const ParameterSet& parameter_set = {}) const {
        if (this->get_parameter_type() == "gate") return _angle;
        if (this->get_parameter_type() == "id") {
            auto it = parameter_set.find(get_parameter_id());
            if (it != parameter_set.end()) return it->second * _parameter_coef;
        }
        throw ParameterIdNotFoundException(
            "Error: QuantumGate_SingleParameter::get_angle(const "
            "ParameterSet&) const: parameter_id \"" +
            get_parameter_id() + "\" is not found");
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

    QuantumGate_SingleParameterOneQubitRotation(double angle)
        : QuantumGate_SingleParameter(angle) {}
    QuantumGate_SingleParameterOneQubitRotation(
        const ParameterId& parameter_id, double parameter_coef = 1.)
        : QuantumGate_SingleParameter(parameter_id, parameter_coef) {}

public:
    virtual void update_quantum_state(
        QuantumStateBase* state, const ParameterSet& parameter_set) {
        double angle = this->get_angle(parameter_set);
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
        if (this->get_parameter_type() != "gate") {
            throw NotImplementedException(
                "Error: "
                "QuantumGate_SingleParameter::update_quantum_state("
                "QuantumStateBase*): "
                "this is an old-style function. You cannot use this function "
                "with a new-style parametric_gate which has string "
                "parameter_id.");
        }
        update_quantum_state(state, {});
    }
};

class ClsParametricRXGate : public QuantumGate_SingleParameterOneQubitRotation {
public:
    ClsParametricRXGate(UINT target_qubit_index, double angle)
        : QuantumGate_SingleParameterOneQubitRotation(angle) {
        this->_name = "ParametricRX";
        this->_update_func = RX_gate;
        this->_update_func_dm = dm_RX_gate;
#ifdef _USE_GPU
        this->_update_func_gpu = RX_gate_host;
#endif
        this->_target_qubit_list.push_back(
            TargetQubitInfo(target_qubit_index, FLAG_X_COMMUTE));
    }
    ClsParametricRXGate(UINT target_qubit_index,
        const ParameterId& parameter_id, double parameter_coef = 1.)
        : QuantumGate_SingleParameterOneQubitRotation(
              parameter_id, parameter_coef) {
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
        double angle = this->get_angle(parameter_set);
        matrix = ComplexMatrix::Zero(2, 2);
        matrix << cos(angle / 2), sin(angle / 2) * 1.i, sin(angle / 2) * 1.i,
            cos(angle / 2);
    }
    virtual void set_matrix(ComplexMatrix& matrix) const {
        if (this->get_parameter_type() != "gate") {
            throw NotImplementedException(
                "Error: "
                "ClsParametricRXGate::set_matrix("
                "ComplexMatrix&): "
                "this is an old-style function. You cannot use this function "
                "with a new-style parametric_gate which has string "
                "parameter_id.");
        }
        set_matrix(matrix, {});
    }
    virtual QuantumGate_SingleParameter* copy() const override {
        return new ClsParametricRXGate(*this);
    };
};

class ClsParametricRYGate : public QuantumGate_SingleParameterOneQubitRotation {
public:
    ClsParametricRYGate(UINT target_qubit_index, double angle)
        : QuantumGate_SingleParameterOneQubitRotation(angle) {
        this->_name = "ParametricRY";
        this->_update_func = RY_gate;
        this->_update_func_dm = dm_RY_gate;
#ifdef _USE_GPU
        this->_update_func_gpu = RY_gate_host;
#endif
        this->_target_qubit_list.push_back(
            TargetQubitInfo(target_qubit_index, FLAG_Y_COMMUTE));
    }
    ClsParametricRYGate(UINT target_qubit_index,
        const ParameterId& parameter_id, double parameter_coef = 1.)
        : QuantumGate_SingleParameterOneQubitRotation(
              parameter_id, parameter_coef) {
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
        double angle = this->get_angle(parameter_set);
        matrix = ComplexMatrix::Zero(2, 2);
        matrix << cos(angle / 2), sin(angle / 2), -sin(angle / 2),
            cos(angle / 2);
    }
    virtual void set_matrix(ComplexMatrix& matrix) const {
        if (this->get_parameter_type() != "gate") {
            throw NotImplementedException(
                "Error: "
                "ClsParametricRYGate::set_matrix("
                "ComplexMatrix&): "
                "this is an old-style function. You cannot use this function "
                "with a new-style parametric_gate which has string "
                "parameter_id.");
        }
        set_matrix(matrix, {});
    }
    virtual QuantumGate_SingleParameter* copy() const override {
        return new ClsParametricRYGate(*this);
    };
};

class ClsParametricRZGate : public QuantumGate_SingleParameterOneQubitRotation {
public:
    ClsParametricRZGate(UINT target_qubit_index, double angle)
        : QuantumGate_SingleParameterOneQubitRotation(angle) {
        this->_name = "ParametricRZ";
        this->_update_func = RZ_gate;
        this->_update_func_dm = dm_RZ_gate;
#ifdef _USE_GPU
        this->_update_func_gpu = RZ_gate_host;
#endif
        this->_target_qubit_list.push_back(
            TargetQubitInfo(target_qubit_index, FLAG_Z_COMMUTE));
    }
    ClsParametricRZGate(UINT target_qubit_index,
        const ParameterId& parameter_id, double parameter_coef = 1.)
        : QuantumGate_SingleParameterOneQubitRotation(
              parameter_id, parameter_coef) {
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
        double angle = this->get_angle(parameter_set);
        matrix = ComplexMatrix::Zero(2, 2);
        matrix << cos(angle / 2) + 1.i * sin(angle / 2), 0, 0,
            cos(angle / 2) - 1.i * sin(angle / 2);
    }
    virtual void set_matrix(ComplexMatrix& matrix) const {
        if (this->get_parameter_type() != "gate") {
            throw NotImplementedException(
                "Error: "
                "ClsParametricRZGate::set_matrix("
                "ComplexMatrix&): "
                "this is an old-style function. You cannot use this function "
                "with a new-style parametric_gate which has string "
                "parameter_id.");
        }
        set_matrix(matrix, {});
    }
    virtual QuantumGate_SingleParameter* copy() const override {
        return new ClsParametricRZGate(*this);
    };
};

class ClsParametricPauliRotationGate : public QuantumGate_SingleParameter {
protected:
    PauliOperator* _pauli;

public:
    ClsParametricPauliRotationGate(double angle, PauliOperator* pauli)
        : QuantumGate_SingleParameter(angle) {
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
    ClsParametricPauliRotationGate(PauliOperator* pauli,
        const ParameterId& parameter_id, double parameter_coef = 1.)
        : QuantumGate_SingleParameter(parameter_id, parameter_coef) {
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
        double angle = this->get_angle(parameter_set);
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
        if (this->get_parameter_type() != "gate") {
            throw NotImplementedException(
                "Error: "
                "ClsParametricPauliRotationGate::update_quantum_state("
                "QuantumStateBase*): "
                "this is an old-style function. You cannot use this function "
                "with a new-style parametric_gate which has string "
                "parameter_id.");
        }
        update_quantum_state(state, {});
    }
    virtual QuantumGate_SingleParameter* copy() const override {
        if (this->get_parameter_type() == "gate") {
            return new ClsParametricPauliRotationGate(_angle, _pauli->copy());
        }
        if (this->get_parameter_type() == "id") {
            return new ClsParametricPauliRotationGate(
                _pauli->copy(), this->get_parameter_id());
        }
        assert(false);
        return nullptr;
    }
    virtual void set_matrix(
        ComplexMatrix& matrix, const ParameterSet& parameter_set) const {
        double angle = this->get_angle(parameter_set);
        get_Pauli_matrix(matrix, _pauli->get_pauli_id_list());
        matrix = cos(angle / 2) *
                     ComplexMatrix::Identity(matrix.rows(), matrix.cols()) +
                 1.i * sin(angle / 2) * matrix;
    };
    virtual void set_matrix(ComplexMatrix& matrix) const {
        if (this->get_parameter_type() != "gate") {
            throw NotImplementedException(
                "Error: "
                "ClsParametricPauliRotationGate::set_matrix(ComplexMatrix&) "
                "const: "
                "this is an old-style function used for parameter_id 'gate:'");
        }
        if (this->get_parameter_type() != "gate") {
            throw NotImplementedException(
                "Error: "
                "ClsParametricPauliRotationGate::set_matrix("
                "ComplexMatrix&): "
                "this is an old-style function. You cannot use this function "
                "with a new-style parametric_gate which has string "
                "parameter_id.");
        }
        set_matrix(matrix, {});
    }
    virtual PauliOperator* get_pauli() const { return _pauli; };
};