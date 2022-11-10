#include "parametric_gate.hpp"

QuantumGate_SingleParameter::QuantumGate_SingleParameter(double angle)
    : _angle(angle) {
    _parameter_id = OLD_STYLE_PARAMETER;
    _parameter_coef = 1.;
    _gate_property |= FLAG_PARAMETRIC;
}
QuantumGate_SingleParameter::QuantumGate_SingleParameter(
    const ParameterId& parameter_id, double parameter_coef)
    : _parameter_id(parameter_id), _parameter_coef(parameter_coef) {
    _gate_property |= FLAG_PARAMETRIC;
}

bool QuantumGate_SingleParameter::is_old_style() const {
    return _parameter_id == OLD_STYLE_PARAMETER;
}
bool QuantumGate_SingleParameter::is_new_style() const {
    return _parameter_id != OLD_STYLE_PARAMETER;
}

ParameterId QuantumGate_SingleParameter::get_parameter_id() const {
    return _parameter_id;
}
void QuantumGate_SingleParameter::set_parameter_id(
    const ParameterId& parameter_id) {
    if (!is_new_style()) {
        throw NotImplementedException(
            "Error: "
            "QuantumGate_SingleParameter::set_parameter_id(const "
            "ParameterId&): "
            "this is a new-style function. You cannot use this function "
            "with an old-style parametric_gate whose parameter_id is -1");
    }
    _parameter_id = parameter_id;
}
double QuantumGate_SingleParameter::get_parameter_coef() const {
    return _parameter_coef;
}
double QuantumGate_SingleParameter::get_parameter_value() const {
    if (!this->is_old_style()) {
        throw NotImplementedException(
            "Error: "
            "QuantumGate_SingleParameter::set_parameter_value(double): "
            "this is an old-style function. You cannot use this function "
            "with a new-style parametric_gate whose parameter_id is not -1");
    }
    return _angle;
}
void QuantumGate_SingleParameter::set_parameter_value(double value) {
    if (!this->is_old_style()) {
        throw NotImplementedException(
            "Error: "
            "QuantumGate_SingleParameter::set_parameter_value(double): "
            "this is an old-style function. You cannot use this function "
            "with a new-style parametric_gate whose parameter_id is not -1");
    }
    _angle = value;
}
double QuantumGate_SingleParameter::get_angle() const {
    if (!this->is_old_style()) {
        throw NotImplementedException(
            "Error: "
            "QuantumGate_SingleParameter::get_angle(): "
            "this is an old-style function. You should use the function "
            "with passing parameter_list");
    }
    return _angle;
}
double QuantumGate_SingleParameter::get_angle(
    const std::vector<double> parameter_list) const {
    if (!this->is_new_style()) {
        throw NotImplementedException(
            "Error: "
            "QuantumGate_SingleParameter::get_angle(const "
            "std::vector<double>&): "
            "this is a new-style function. You should use the function "
            "without passing parameter_list");
    }
    if (this->_parameter_id >= parameter_list.size()) {
        throw ParameterIndexOutOfRangeException(
            "Error: QuantumGate_SingleParameter::get_angle(const "
            "std::vector<double>&) const: parameter_list index out of range");
    }
    return parameter_list[_parameter_id] * _parameter_coef;
}

QuantumGate_SingleParameterOneQubitRotation::
    QuantumGate_SingleParameterOneQubitRotation(double angle)
    : QuantumGate_SingleParameter(angle) {}
QuantumGate_SingleParameterOneQubitRotation::
    QuantumGate_SingleParameterOneQubitRotation(
        const ParameterId& parameter_id, double parameter_coef)
    : QuantumGate_SingleParameter(parameter_id, parameter_coef) {}

void QuantumGate_SingleParameterOneQubitRotation::_update_quantum_state(
    QuantumStateBase* state, double angle) {
    if (state->is_state_vector()) {
#ifdef _USE_GPU
        if (state->get_device_name() == "gpu") {
            if (_update_func_gpu == nullptr) {
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
        if (_update_func == nullptr) {
            throw UndefinedUpdateFuncException(
                "Error: "
                "QuantumGate_SingleParameterOneQubitRotation::update_"
                "quantum_state(QuantumStateBase) : update function is "
                "undefined");
        }
        _update_func(this->_target_qubit_list[0].index(), angle,
            state->data_c(), state->dim);
    } else {
        if (_update_func_dm == nullptr) {
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
void QuantumGate_SingleParameterOneQubitRotation::update_quantum_state(
    QuantumStateBase* state) {
    if (!this->is_old_style()) {
        throw NotImplementedException(
            "Error: "
            "QuantumGate_SingleParameter::update_qunautm_state("
            "QuantumStateBase&, const std::vector<double>&): "
            "this is an old-style function. You should use the function "
            "with passing parameter_list");
    }
    this->_update_quantum_state(state, this->get_angle());
}
void QuantumGate_SingleParameterOneQubitRotation::update_quantum_state(
    QuantumStateBase* state, const std::vector<double>& parameter_list) {
    if (!this->is_new_style()) {
        throw NotImplementedException(
            "Error: "
            "QuantumGate_SingleParameter::update_qunautm_state("
            "QuantumStateBase&, const std::vector<double>&): "
            "this is a new-style function. You should use the function "
            "without passing parameter_list");
    }
    this->_update_quantum_state(state, this->get_angle(parameter_list));
}

ClsParametricRXGate::ClsParametricRXGate(UINT target_qubit_index, double angle)
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
ClsParametricRXGate::ClsParametricRXGate(UINT target_qubit_index,
    const ParameterId& parameter_id, double parameter_coef)
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

void ClsParametricRXGate::_set_matrix(
    ComplexMatrix& matrix, double angle) const {
    matrix = ComplexMatrix::Zero(2, 2);
    matrix << cos(angle / 2), sin(angle / 2) * 1.i, sin(angle / 2) * 1.i,
        cos(angle / 2);
}
void ClsParametricRXGate::set_matrix(ComplexMatrix& matrix) const {
    if (!this->is_old_style()) {
        throw NotImplementedException(
            "Error: "
            "ClsParametricRXGate::set_matrix(ComplexMatrix& matrix): "
            "this is an old-style function. You should use the function "
            "with passing parameter_list");
    }
    _set_matrix(matrix, this->get_angle());
}
void ClsParametricRXGate::set_matrix(
    ComplexMatrix& matrix, const std::vector<double>& parameter_list) const {
    if (!this->is_new_style()) {
        throw NotImplementedException(
            "Error: "
            "ClsParametricRXGate::set_matrix(ComplexMatrix& matrix): "
            "this is a new-style function. You should use the function "
            "without passing parameter_list");
    }

    _set_matrix(matrix, this->get_angle(parameter_list));
}
QuantumGate_SingleParameter* ClsParametricRXGate::copy() const {
    return new ClsParametricRXGate(*this);
};

ClsParametricRYGate::ClsParametricRYGate(UINT target_qubit_index, double angle)
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
ClsParametricRYGate::ClsParametricRYGate(UINT target_qubit_index,
    const ParameterId& parameter_id, double parameter_coef)
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

void ClsParametricRYGate::_set_matrix(
    ComplexMatrix& matrix, double angle) const {
    matrix = ComplexMatrix::Zero(2, 2);
    matrix << cos(angle / 2), sin(angle / 2), -sin(angle / 2), cos(angle / 2);
}
void ClsParametricRYGate::set_matrix(ComplexMatrix& matrix) const {
    if (!this->is_old_style()) {
        throw NotImplementedException(
            "Error: "
            "ClsParametricRYGate::set_matrix(ComplexMatrix& matrix): "
            "this is an old-style function. You should use the function "
            "with passing parameter_list");
    }
    _set_matrix(matrix, this->get_angle());
}
void ClsParametricRYGate::set_matrix(
    ComplexMatrix& matrix, const std::vector<double>& parameter_list) const {
    if (!this->is_new_style()) {
        throw NotImplementedException(
            "Error: "
            "ClsParametricRYGate::set_matrix(ComplexMatrix& matrix): "
            "this is a new-style function. You should use the function "
            "without passing parameter_list");
    }

    _set_matrix(matrix, this->get_angle(parameter_list));
}
QuantumGate_SingleParameter* ClsParametricRYGate::copy() const {
    return new ClsParametricRYGate(*this);
};

ClsParametricRZGate::ClsParametricRZGate(UINT target_qubit_index, double angle)
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
ClsParametricRZGate::ClsParametricRZGate(UINT target_qubit_index,
    const ParameterId& parameter_id, double parameter_coef)
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

void ClsParametricRZGate::_set_matrix(
    ComplexMatrix& matrix, double angle) const {
    matrix = ComplexMatrix::Zero(2, 2);
    matrix << cos(angle / 2) + 1.i * sin(angle / 2), 0, 0,
        cos(angle / 2) - 1.i * sin(angle / 2);
}
void ClsParametricRZGate::set_matrix(ComplexMatrix& matrix) const {
    if (!this->is_old_style()) {
        throw NotImplementedException(
            "Error: "
            "ClsParametricRZGate::set_matrix(ComplexMatrix& matrix): "
            "this is an old-style function. You should use the function "
            "with passing parameter_list");
    }
    this->_set_matrix(matrix, this->get_angle());
}
void ClsParametricRZGate::set_matrix(
    ComplexMatrix& matrix, const std::vector<double>& parameter_list) const {
    if (!this->is_new_style()) {
        throw NotImplementedException(
            "Error: "
            "ClsParametricRZGate::set_matrix(ComplexMatrix& matrix): "
            "this is a new-style function. You should use the function "
            "without passing parameter_list");
    }

    this->_set_matrix(matrix, this->get_angle(parameter_list));
}
QuantumGate_SingleParameter* ClsParametricRZGate::copy() const {
    return new ClsParametricRZGate(*this);
};

ClsParametricPauliRotationGate::ClsParametricPauliRotationGate(
    double angle, PauliOperator* pauli)
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
        this->_target_qubit_list.push_back(
            TargetQubitInfo(target_index_list[index], commutation_relation));
    }
}
ClsParametricPauliRotationGate::ClsParametricPauliRotationGate(
    PauliOperator* pauli, const ParameterId& parameter_id,
    double parameter_coef)
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
        this->_target_qubit_list.push_back(
            TargetQubitInfo(target_index_list[index], commutation_relation));
    }
}
ClsParametricPauliRotationGate::~ClsParametricPauliRotationGate() {
    delete _pauli;
}

void ClsParametricPauliRotationGate::_update_quantum_state(
    QuantumStateBase* state, double angle) {
    auto target_index_list = _pauli->get_index_list();
    auto pauli_id_list = _pauli->get_pauli_id_list();
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
        multi_qubit_Pauli_rotation_gate_partial_list(target_index_list.data(),
            pauli_id_list.data(), (UINT)target_index_list.size(), angle,
            state->data_c(), state->dim);
#endif
    } else {
        dm_multi_qubit_Pauli_rotation_gate_partial_list(
            target_index_list.data(), pauli_id_list.data(),
            (UINT)target_index_list.size(), angle, state->data_c(), state->dim);
    }
};
void ClsParametricPauliRotationGate::update_quantum_state(
    QuantumStateBase* state) {
    if (!this->is_old_style()) {
        throw NotImplementedException(
            "Error: "
            "ClsParametricPauliRotationGate::update_qunautm_state("
            "QuantumStateBase&, const std::vector<double>&): "
            "this is an old-style function. You should use the function "
            "with passing parameter_list");
    }
    this->_update_quantum_state(state, this->get_angle());
}
void ClsParametricPauliRotationGate::update_quantum_state(
    QuantumStateBase* state, const std::vector<double>& parameter_list) {
    if (!this->is_new_style()) {
        throw NotImplementedException(
            "Error: "
            "ClsParametricPauliRotationGate::update_qunautm_state("
            "QuantumStateBase&, const std::vector<double>&): "
            "this is a new-style function. You should use the function "
            "without passing parameter_list");
    }
    this->_update_quantum_state(state, this->get_angle(parameter_list));
}

QuantumGate_SingleParameter* ClsParametricPauliRotationGate::copy() const {
    if (this->is_old_style()) {
        return new ClsParametricPauliRotationGate(_angle, _pauli->copy());
    }
    if (this->is_new_style()) {
        return new ClsParametricPauliRotationGate(
            _pauli->copy(), this->get_parameter_id());
    }
    assert(false);
    return nullptr;
}

void ClsParametricPauliRotationGate::_set_matrix(
    ComplexMatrix& matrix, double angle) const {
    get_Pauli_matrix(matrix, _pauli->get_pauli_id_list());
    matrix =
        cos(angle / 2) * ComplexMatrix::Identity(matrix.rows(), matrix.cols()) +
        1.i * sin(angle / 2) * matrix;
}
void ClsParametricPauliRotationGate::set_matrix(ComplexMatrix& matrix) const {
    if (!this->is_old_style()) {
        throw NotImplementedException(
            "Error: "
            "ClsParametricPauliRotationGate::set_matrix(ComplexMatrix& "
            "matrix): "
            "this is an old-style function. You should use the function "
            "with passing parameter_list");
    }
    this->_set_matrix(matrix, this->get_angle());
}
void ClsParametricPauliRotationGate::set_matrix(
    ComplexMatrix& matrix, const std::vector<double>& parameter_list) const {
    if (!this->is_new_style()) {
        throw NotImplementedException(
            "Error: "
            "ClsParametricPauliRotationGate::set_matrix(ComplexMatrix& "
            "matrix): "
            "this is a new-style function. You should use the function "
            "without passing parameter_list");
    }
    this->_set_matrix(matrix, this->get_angle(parameter_list));
}

PauliOperator* ClsParametricPauliRotationGate::get_pauli() const {
    return this->_pauli;
}
