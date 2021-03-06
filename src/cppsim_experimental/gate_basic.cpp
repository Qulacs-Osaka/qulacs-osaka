
#include "gate_basic.hpp"

void QuantumGateBasic::_update_state_vector_cpu_special(
    QuantumStateBase* state) const {
    if (_special_func_type == GateI) {
        // pass
    } else if (_special_func_type == GateX) {
        X_gate(_target_qubit_index[0], state->data_c(), state->dim);
    } else if (_special_func_type == GateY) {
        Y_gate(_target_qubit_index[0], state->data_c(), state->dim);
    } else if (_special_func_type == GateZ) {
        Z_gate(_target_qubit_index[0], state->data_c(), state->dim);
    } else if (_special_func_type == GateSqrtX) {
        sqrtX_gate(_target_qubit_index[0], state->data_c(), state->dim);
    } else if (_special_func_type == GateSqrtXdag) {
        sqrtXdag_gate(_target_qubit_index[0], state->data_c(), state->dim);
    } else if (_special_func_type == GateSqrtY) {
        sqrtY_gate(_target_qubit_index[0], state->data_c(), state->dim);
    } else if (_special_func_type == GateSqrtYdag) {
        sqrtYdag_gate(_target_qubit_index[0], state->data_c(), state->dim);
    } else if (_special_func_type == GateH) {
        H_gate(_target_qubit_index[0], state->data_c(), state->dim);
    } else if (_special_func_type == GateS) {
        S_gate(_target_qubit_index[0], state->data_c(), state->dim);
    } else if (_special_func_type == GateSdag) {
        Sdag_gate(_target_qubit_index[0], state->data_c(), state->dim);
    } else if (_special_func_type == GateT) {
        T_gate(_target_qubit_index[0], state->data_c(), state->dim);
    } else if (_special_func_type == GateTdag) {
        Tdag_gate(_target_qubit_index[0], state->data_c(), state->dim);
    } else if (_special_func_type == GateP0) {
        P0_gate(_target_qubit_index[0], state->data_c(), state->dim);
    } else if (_special_func_type == GateP1) {
        P1_gate(_target_qubit_index[0], state->data_c(), state->dim);
    } else if (_special_func_type == GateRX) {
        // invert
        RX_gate(_target_qubit_index[0], -_rotation_angle, state->data_c(),
            state->dim);
    } else if (_special_func_type == GateRY) {
        // invert
        RY_gate(_target_qubit_index[0], -_rotation_angle, state->data_c(),
            state->dim);
    } else if (_special_func_type == GateRZ) {
        // invert
        RZ_gate(_target_qubit_index[0], -_rotation_angle, state->data_c(),
            state->dim);
    } else if (_special_func_type == GateCX) {
        CNOT_gate(_control_qubit_index[0], _target_qubit_index[0],
            state->data_c(), state->dim);
    } else if (_special_func_type == GateCZ) {
        CZ_gate(_control_qubit_index[0], _target_qubit_index[0],
            state->data_c(), state->dim);
    } else if (_special_func_type == GateSWAP) {
        SWAP_gate(_target_qubit_index[0], _target_qubit_index[1],
            state->data_c(), state->dim);
    } else {
        throw std::invalid_argument("Unsupported special gate");
    }
}

void QuantumGateBasic::_update_state_vector_cpu_general(
    QuantumStateBase* state) const {
    if (_matrix_type == DenseMatrix) {
        const CTYPE* matrix_ptr =
            reinterpret_cast<const CTYPE*>(this->_dense_matrix_element.data());
        // single qubit dense matrix gate
        if (_target_qubit_index.size() == 1) {
            // no control qubit
            if (_control_qubit_index.size() == 0) {
                single_qubit_dense_matrix_gate(_target_qubit_index[0],
                    matrix_ptr, state->data_c(), state->dim);
            }
            // single control qubit
            else if (_control_qubit_index.size() == 1) {
                single_qubit_control_single_qubit_dense_matrix_gate(
                    _control_qubit_index[0], _control_qubit_value[0],
                    _target_qubit_index[0], matrix_ptr, state->data_c(),
                    state->dim);
            }
            // multiple control qubits
            else {
                multi_qubit_control_single_qubit_dense_matrix_gate(
                    _control_qubit_index.data(), _control_qubit_value.data(),
                    (UINT)(_control_qubit_index.size()), _target_qubit_index[0],
                    matrix_ptr, state->data_c(), state->dim);
            }
        }

        // multi qubit dense matrix gate
        else {
            // no control qubit
            if (_control_qubit_index.size() == 0) {
                multi_qubit_dense_matrix_gate(_target_qubit_index.data(),
                    (UINT)(_target_qubit_index.size()), matrix_ptr,
                    state->data_c(), state->dim);
            }
            // single control qubit
            else if (_control_qubit_index.size() == 1) {
                single_qubit_control_multi_qubit_dense_matrix_gate(
                    _control_qubit_index[0], _control_qubit_value[0],
                    _target_qubit_index.data(),
                    (UINT)(_target_qubit_index.size()), matrix_ptr,
                    state->data_c(), state->dim);
            }
            // multiple control qubit
            else {
                multi_qubit_control_multi_qubit_dense_matrix_gate(
                    _control_qubit_index.data(), _control_qubit_value.data(),
                    (UINT)(_control_qubit_index.size()),
                    _target_qubit_index.data(),
                    (UINT)(_target_qubit_index.size()), matrix_ptr,
                    state->data_c(), state->dim);
            }
        }
    } else if (_matrix_type == DiagonalMatrix) {
        const CTYPE* matrix_ptr = reinterpret_cast<const CTYPE*>(
            this->_diagonal_matrix_element.data());
        if (_target_qubit_index.size() == 1)
            single_qubit_diagonal_matrix_gate(_target_qubit_index[0],
                matrix_ptr, state->data_c(), state->dim);
        else
            multi_qubit_diagonal_matrix_gate(_target_qubit_index.data(),
                (UINT)_target_qubit_index.size(), matrix_ptr, state->data_c(),
                state->dim);
    } else if (_matrix_type == SparseMatrix) {
        multi_qubit_sparse_matrix_gate_eigen(_target_qubit_index.data(),
            (UINT)(_target_qubit_index.size()), this->_sparse_matrix_element,
            state->data_c(), state->dim);
    } else if (_matrix_type == PauliMatrix) {
        if (_target_qubit_index.size() == 1) {
            if (fabs(_rotation_angle) < 1e-16) {
                single_qubit_Pauli_gate(_target_qubit_index[0], _pauli_id[0],
                    state->data_c(), state->dim);
            } else {
                // invert
                single_qubit_Pauli_rotation_gate(_target_qubit_index[0],
                    _pauli_id[0], -_rotation_angle, state->data_c(),
                    state->dim);
            }
        } else {
            if (fabs(_rotation_angle) < 1e-16) {
                multi_qubit_Pauli_gate_partial_list(_target_qubit_index.data(),
                    _pauli_id.data(), (UINT)_target_qubit_index.size(),
                    state->data_c(), state->dim);
            } else {
                // invert
                multi_qubit_Pauli_rotation_gate_partial_list(
                    _target_qubit_index.data(), _pauli_id.data(),
                    (UINT)_target_qubit_index.size(), -_rotation_angle,
                    state->data_c(), state->dim);
            }
        }
    }
}

void QuantumGateBasic::_update_density_matrix_cpu_general(
    QuantumStateBase* state) const {
    if (_matrix_type == DenseMatrix) {
        const CTYPE* matrix_ptr =
            reinterpret_cast<const CTYPE*>(this->_dense_matrix_element.data());
        if (_control_qubit_index.size() == 0) {
            if (_target_qubit_index.size() == 1) {
                dm_single_qubit_dense_matrix_gate(_target_qubit_index[0],
                    matrix_ptr, state->data_c(), state->dim);
            } else {
                dm_multi_qubit_dense_matrix_gate(_target_qubit_index.data(),
                    (UINT)_target_qubit_index.size(), matrix_ptr,
                    state->data_c(), state->dim);
            }
        } else {
            if (_target_qubit_index.size() == 1) {
                dm_multi_qubit_control_single_qubit_dense_matrix_gate(
                    _control_qubit_index.data(), _control_qubit_value.data(),
                    (UINT)_control_qubit_index.size(), _target_qubit_index[0],
                    matrix_ptr, state->data_c(), state->dim);
            } else {
                dm_multi_qubit_control_multi_qubit_dense_matrix_gate(
                    _control_qubit_index.data(), _control_qubit_value.data(),
                    (UINT)_control_qubit_index.size(),
                    _target_qubit_index.data(),
                    (UINT)_target_qubit_index.size(), matrix_ptr,
                    state->data_c(), state->dim);
            }
        }
    } else {
        throw std::invalid_argument(
            "Only DenseMatrix gate type is supported for density matrix");
    }
}

#ifdef _USE_GPU
void QuantumGateBasic::_update_state_vector_gpu(QuantumStateBase* state) {
    if (_matrix_type == DenseMatrix) {
        const CTYPE* matrix_ptr =
            reinterpret_cast<const CTYPE*>(this->_dense_matrix_element.data());
        // single qubit dense matrix gate
        if (_target_qubit_index.size() == 1) {
            // no control qubit
            if (_control_qubit_index.size() == 0) {
                single_qubit_dense_matrix_gate_host(_target_qubit_index[0],
                    (const CPPCTYPE*)matrix_ptr, state->data(), state->dim,
                    state->get_cuda_stream(), state->device_number);
            }
            // single control qubit
            else if (_control_qubit_index.size() == 1) {
                single_qubit_control_single_qubit_dense_matrix_gate_host(
                    _control_qubit_index[0], _control_qubit_value[0],
                    _target_qubit_index[0], (const CPPCTYPE*)matrix_ptr,
                    state->data(), state->dim, state->get_cuda_stream(),
                    state->device_number);
            }
            // multiple control qubits
            else {
                multi_qubit_control_multi_qubit_dense_matrix_gate_host(
                    _control_qubit_index.data(), _control_qubit_value.data(),
                    (UINT)(_control_qubit_index.size()),
                    _target_qubit_index.data(),
                    (UINT)(_target_qubit_index.size()),
                    (const CPPCTYPE*)matrix_ptr, state->data(), state->dim,
                    state->get_cuda_stream(), state->device_number);
            }
        }

        // multi qubit dense matrix gate
        else {
            // no control qubit
            if (_control_qubit_index.size() == 0) {
                multi_qubit_dense_matrix_gate_host(_target_qubit_index.data(),
                    (UINT)(_target_qubit_index.size()),
                    (const CPPCTYPE*)matrix_ptr, state->data(), state->dim,
                    state->get_cuda_stream(), state->device_number);
            }
            // single control qubit
            else if (_control_qubit_index.size() == 1) {
                single_qubit_control_multi_qubit_dense_matrix_gate_host(
                    _control_qubit_index[0], _control_qubit_value[0],
                    _target_qubit_index.data(),
                    (UINT)(_target_qubit_index.size()),
                    (const CPPCTYPE*)matrix_ptr, state->data(), state->dim,
                    state->get_cuda_stream(), state->device_number);
            }
            // multiple control qubit
            else {
                multi_qubit_control_multi_qubit_dense_matrix_gate_host(
                    _control_qubit_index.data(), _control_qubit_value.data(),
                    (UINT)(_control_qubit_index.size()),
                    _target_qubit_index.data(),
                    (UINT)(_target_qubit_index.size()),
                    (const CPPCTYPE*)matrix_ptr, state->data(), state->dim,
                    state->get_cuda_stream(), state->device_number);
            }
        }
    } else {
        throw std::invalid_argument(
            "Only DenseMatrix gate type is supported for density matrix");
    }
}
void _update_density_matrix_gpu(QuantumStateBase* state) {
    throw std::runtime_error(
        "Density matrix simulation is not supported on GPU.");
}
#endif

void get_new_qubit_list(const QuantumGateBase* gate_first,
    const QuantumGateBase* gate_second,
    std::vector<UINT>& new_target_index_list,
    std::vector<UINT>& new_target_commutation_list,
    std::vector<UINT>& new_control_index_list,
    std::vector<UINT>& new_control_value_list);
void get_extended_matrix(const QuantumGateBase* gate,
    const std::vector<UINT>& new_target_index_list,
    const std::vector<UINT>& new_control_index_list, ComplexMatrix& matrix);

// Create target_gate_set and control_gate_set after merging
// Any qubit index is classified as 9 cases :  (first_target, first_control,
// not_in_first) * (second_target, second_control, not_in_second) Since each
// target_qubit_list is not sorted and since we cannot sort them without
// corrupsing matrix correspondense, we cannot use stl set functions. Currently,
// all the indices are classified with a dirty way.
void get_new_qubit_list(const QuantumGateBase* gate_first,
    const QuantumGateBase* gate_second,
    std::vector<UINT>& new_target_index_list,
    std::vector<UINT>& new_target_commutation_list,
    std::vector<UINT>& new_control_index_list,
    std::vector<UINT>& new_control_value_list) {
    std::vector<UINT> gate_first_target_index =
        gate_first->get_target_index_list();
    std::vector<UINT> gate_first_target_commutation =
        gate_first->get_target_commutation_list();
    std::vector<UINT> gate_first_control_index =
        gate_first->get_control_index_list();
    std::vector<UINT> gate_first_control_value =
        gate_first->get_control_value_list();
    std::vector<UINT> gate_second_target_index =
        gate_second->get_target_index_list();
    std::vector<UINT> gate_second_target_commutation =
        gate_second->get_target_commutation_list();
    std::vector<UINT> gate_second_control_index =
        gate_second->get_control_index_list();
    std::vector<UINT> gate_second_control_value =
        gate_second->get_control_value_list();
    for (UINT i = 0; i < gate_first_target_index.size(); ++i) {
        // case 0 : qubit belongs to both target_set -> property is merged,
        // goto new_target_set
        UINT gate_first_target_qubit = gate_first_target_index[i];
        auto ite_target = std::find_if(gate_second_target_index.begin(),
            gate_second_target_index.end(),
            [&gate_first_target_qubit](UINT gate_second_target_qubit) {
                return gate_second_target_qubit == gate_first_target_qubit;
            });
        if (ite_target != gate_second_target_index.end()) {
            UINT gate_second_list_index =
                std::distance(gate_second_target_index.begin(), ite_target);
            new_target_index_list.push_back(gate_first_target_qubit);
            new_target_commutation_list.push_back(
                gate_first_target_commutation[i] &
                gate_second_target_commutation[gate_second_list_index]);
            continue;
        }

        // case 1 : qubit belongs to first gate and second control -> first
        // property is merged to Z, goto new_target_set
        auto ite_control = std::find_if(gate_second_control_index.begin(),
            gate_second_control_index.end(),
            [&gate_first_target_qubit](UINT gate_second_control_qubit) {
                return gate_first_target_qubit == gate_second_control_qubit;
            });
        if (ite_control != gate_second_control_index.end()) {
            new_target_index_list.push_back(gate_first_target_qubit);
            new_target_commutation_list.push_back(
                gate_first_target_commutation[i] & FLAG_COMMUTE_Z);
            continue;
        }

        // case 2 : qubit belongs to first gate and not in second -> first
        // property is preserved, goto new_target_set
        else {
            new_target_index_list.push_back(gate_first_target_qubit);
            new_target_commutation_list.push_back(
                gate_first_target_commutation[i]);
        }
    }

    for (UINT i = 0; i < gate_first_control_index.size(); ++i) {
        // case 3 : qubit belongs to first control and second target -> second
        // property is merged with Z, goto new_target_set
        UINT gate_first_control_qubit = gate_first_control_index[i];
        auto ite_target = std::find_if(gate_second_target_index.begin(),
            gate_second_target_index.end(),
            [&gate_first_control_qubit](UINT gate_second_target_qubit) {
                return gate_first_control_qubit == gate_second_target_qubit;
            });
        if (ite_target != gate_second_target_index.end()) {
            UINT list_index =
                std::distance(gate_second_target_index.begin(), ite_target);
            new_target_index_list.push_back(gate_first_control_qubit);
            new_target_commutation_list.push_back(
                gate_second_target_commutation[list_index] & FLAG_COMMUTE_Z);
            continue;
        }

        // case 4 : qubit belongs to first control and second control ->  if
        // control_value is equal, goto new_control_set. If not, goto
        // new_target_set with Z_COMMUTE
        auto ite_control = std::find_if(gate_second_control_index.begin(),
            gate_second_control_index.end(),
            [&gate_first_control_qubit](UINT gate_second_control_qubit) {
                return gate_first_control_qubit == gate_second_control_qubit;
            });
        if (ite_control != gate_second_control_index.end()) {
            UINT gate_first_list_index = i;
            UINT gate_second_list_index =
                std::distance(gate_second_control_index.begin(), ite_control);

            if (gate_first_control_value[gate_first_list_index] ==
                gate_second_control_value[gate_second_list_index]) {
                new_control_index_list.push_back(gate_first_control_qubit);
                new_control_value_list.push_back(
                    gate_first_control_value[gate_first_list_index]);
            } else {
                new_target_index_list.push_back(gate_first_control_qubit);
                new_target_commutation_list.push_back(FLAG_COMMUTE_Z);
            }
            continue;
        }

        // case 5 : qubit belongs to first control and not in second -> goto
        // new_target_set with Z_COMMUTE
        else {
            new_target_index_list.push_back(gate_first_control_qubit);
            new_target_commutation_list.push_back(FLAG_COMMUTE_Z);
        }
    }

    for (UINT i = 0; i < gate_second_target_index.size(); ++i) {
        UINT gate_second_target_qubit = gate_second_target_index[i];
        auto ite_target = std::find_if(gate_first_target_index.begin(),
            gate_first_target_index.end(),
            [&gate_second_target_qubit](UINT gate_first_target_qubit) {
                return gate_first_target_qubit == gate_second_target_qubit;
            });
        if (ite_target != gate_first_target_index.end()) {
            continue;
        }
        auto ite_control = std::find_if(gate_first_control_index.begin(),
            gate_first_control_index.end(),
            [&gate_second_target_qubit](UINT gate_first_control_qubit) {
                return gate_first_control_qubit == gate_second_target_qubit;
            });
        if (ite_control != gate_first_control_index.end()) {
            continue;
        }

        // case 6 : qubit belongs to second target but not in first -> goto
        // new_target_set with second property
        else {
            new_target_index_list.push_back(gate_second_target_qubit);
            new_target_commutation_list.push_back(
                gate_second_target_commutation[i]);
        }
    }
    for (UINT i = 0; i < gate_second_control_index.size(); ++i) {
        UINT gate_second_control_qubit = gate_second_control_index[i];
        auto ite_target = std::find_if(gate_first_target_index.begin(),
            gate_first_target_index.end(),
            [&gate_second_control_qubit](UINT gate_first_target_qubit) {
                return gate_first_target_qubit == gate_second_control_qubit;
            });
        if (ite_target != gate_first_target_index.end()) {
            continue;
        }
        auto ite_control = std::find_if(gate_first_control_index.begin(),
            gate_first_control_index.end(),
            [&gate_second_control_qubit](UINT gate_first_control_qubit) {
                return gate_first_control_qubit == gate_second_control_qubit;
            });
        if (ite_control != gate_first_control_index.end()) {
            continue;
        }

        // case 7 : qubit belongs to second control but not in first -> goto
        // new_target_set with Z_COMMUTE
        else {
            new_target_index_list.push_back(gate_second_control_qubit);
            new_target_commutation_list.push_back(FLAG_COMMUTE_Z);
        }
    }

    // case 8 : qubit belongs to nothing -> do nothing
}

// Join new qubit indices to target_list according to a given new_target_list,
// and set a new matrix to "matrix"
void get_extended_matrix(const QuantumGateBase* gate,
    const std::vector<UINT>& new_target_index_list, const std::vector<UINT>&,
    ComplexMatrix& matrix) {
    // New qubits index may be in either gate_target_index, gate_control_index,
    // or it comes from the other gate. Case 0 : If qubit index is in
    // gate_target_index -> named A = original gate_target_index (Order must not
    // be changed!!!)
    std::vector<UINT> join_from_target = gate->get_target_index_list();

    // Case 1 : If qubit index is in gate_control_index -> named B
    std::vector<UINT> join_from_control;
    ITYPE control_mask = 0;

    std::vector<UINT> gate_target_index_list = gate->get_target_index_list();
    std::vector<UINT> gate_control_index_list = gate->get_control_index_list();
    std::vector<UINT> gate_control_value_list = gate->get_control_value_list();
    for (UINT i = 0; i < new_target_index_list.size(); ++i) {
        UINT val = new_target_index_list[i];
        auto ite = std::find_if(gate_control_index_list.begin(),
            gate_control_index_list.end(),
            [&val](const UINT info) { return info == val; });
        if (ite != gate_control_index_list.end()) {
            int list_index =
                std::distance(gate_control_index_list.begin(), ite);
            join_from_control.push_back(gate_control_index_list[list_index]);

            if (gate_control_value_list[list_index] == 1)
                control_mask ^= (1ULL << (join_from_control.size() - 1));
        }
    }
    // Case 2 : If qubit index is not in both -> named C
    std::vector<UINT> join_from_other_gate;
    for (UINT i = 0; i < new_target_index_list.size(); ++i) {
        UINT val = new_target_index_list[i];
        auto ite1 = std::find_if(gate_target_index_list.begin(),
            gate_target_index_list.end(),
            [&val](const UINT info) { return info == val; });
        auto ite2 = std::find_if(gate_control_index_list.begin(),
            gate_control_index_list.end(),
            [&val](const UINT info) { return info == val; });
        if (ite1 == gate_target_index_list.end() &&
            ite2 == gate_control_index_list.end()) {
            join_from_other_gate.push_back(val);
        }
    }
    // At first, qubit indices are ordered as (A,C,B)
    std::vector<UINT> unsorted_new_target_index_list = join_from_target;
    unsorted_new_target_index_list.insert(unsorted_new_target_index_list.end(),
        join_from_other_gate.begin(), join_from_other_gate.end());
    unsorted_new_target_index_list.insert(unsorted_new_target_index_list.end(),
        join_from_control.begin(), join_from_control.end());

    // *** NOTE ***
    // Order of tensor product is reversed!!!
    // U_0 I_1 = I \tensor U = [[U,0], [0,U]]
    // 0-control-U_0 = |0><0| \tensor U + |1><1| \tensor I = [[U,0],[0,I]]
    // 1-control-U_0 = |0><0| \tensor I + |1><1| \tensor U = [[I,0],[0,U]]

    // *** Algorithm ***
    // The gate matrix corresponding to indices (A,C,B) has 2^|B| blocks of
    // gate matrix with (A,C). The (control_mask)-th block matrix is (A,C),
    // and the others are Identity. The gate matrix with (A,C) has 2^|C|
    // blocks of gate matrix with A, which is equal to the original gate
    // matrix.

    // Thus, the following steps work.
    // 1. Enumerate set B and C. -> Done
    // 2. Generate 2^{|A|+|B|+|C|}-dim identity matrix
    size_t new_matrix_qubit_count = (UINT)new_target_index_list.size();
    size_t new_matrix_dim = 1ULL << new_matrix_qubit_count;
    matrix = ComplexMatrix::Identity(new_matrix_dim, new_matrix_dim);
    // 3. Decide correct 2^{|A|+|C|}-dim block matrix from control values.
    ITYPE start_block_basis =
        (1ULL << (join_from_target.size() + join_from_other_gate.size())) *
        control_mask;

    // 4. Repeat 2^{|C|}-times paste of original gate matrix A .

    ComplexMatrix org_matrix;
    gate->get_matrix(org_matrix);

    size_t org_matrix_dim = 1ULL << gate_target_index_list.size();
    ITYPE repeat_count = 1ULL << join_from_other_gate.size();
    for (ITYPE repeat_index = 0; repeat_index < repeat_count; ++repeat_index) {
        size_t paste_start =
            (size_t)(start_block_basis + repeat_index * org_matrix_dim);
        matrix.block(paste_start, paste_start, org_matrix_dim, org_matrix_dim) =
            org_matrix;
    }
    // 5. Since the order of (C,B,A) is different from that of the other
    // gate, we sort (C,B,A) after generating matrix. We do nothing if it is
    // already sorted
    if (!std::is_sorted(unsorted_new_target_index_list.begin(),
            unsorted_new_target_index_list.end())) {
        // generate ascending index of the INDEX_NUMBER of
        // unsorted_target_qubit_index_list.
        std::vector<std::pair<UINT, UINT>> sorted_element_position;
        for (UINT i = 0; i < unsorted_new_target_index_list.size(); ++i) {
            sorted_element_position.push_back(
                std::make_pair(unsorted_new_target_index_list[i], i));
        }
        std::sort(
            sorted_element_position.begin(), sorted_element_position.end());
        std::vector<UINT> sorted_index(sorted_element_position.size(), -1);
        for (UINT i = 0; i < sorted_index.size(); ++i)
            sorted_index[sorted_element_position[i].second] = i;

        // If target qubit is not in the sorted position, we swap the
        // element to the element in correct position. If not, go next
        // index. This sort will stop with n-time swap in the worst case,
        // which is smaller than the cost of std::sort. We cannot directly
        // sort target qubit list in order to swap matrix rows and columns
        // with respect to qubit ordering.
        UINT ind1 = 0;
        while (ind1 < sorted_index.size()) {
            if (sorted_index[ind1] != ind1) {
                UINT ind2 = sorted_index[ind1];

                // move to correct position
                std::swap(sorted_index[ind1], sorted_index[ind2]);
                std::swap(unsorted_new_target_index_list[ind1],
                    unsorted_new_target_index_list[ind2]);

                // create masks
                const UINT min_index = std::min(ind1, ind2);
                const UINT max_index = std::max(ind1, ind2);
                const ITYPE min_mask = 1ULL << min_index;
                const ITYPE max_mask = 1ULL << max_index;

                const ITYPE loop_dim = new_matrix_dim >> 2;

                for (ITYPE state_index = 0; state_index < loop_dim;
                     ++state_index) {
                    ITYPE basis_00 = state_index;
                    basis_00 = insert_zero_to_basis_index(
                        basis_00, min_mask, min_index);
                    basis_00 = insert_zero_to_basis_index(
                        basis_00, max_mask, max_index);
                    ITYPE basis_01 = basis_00 ^ min_mask;
                    ITYPE basis_10 = basis_00 ^ max_mask;

                    matrix.col((size_t)basis_01)
                        .swap(matrix.col((size_t)basis_10));
                    matrix.row((size_t)basis_01)
                        .swap(matrix.row((size_t)basis_10));
                }
            } else
                ind1++;
        }
    }

    // std::cout << "unsorted " << std::endl;
    // for (auto val : unsorted_target_list) std::cout << val << " ";
    // std::cout
    // << std::endl; std::cout << matrix << std::endl;
    // sort_target_qubit(unsorted_new_target_index_list, matrix);
    // std::cout << "sorted " << std::endl;
    // for (auto val : unsorted_target_list) std::cout << val << " ";
    // std::cout
    // << std::endl; std::cout << matrix << std::endl;
}

namespace gate {
DllExport QuantumGateBasic* merge(
    const QuantumGateBase* gate_first, const QuantumGateBase* gate_second) {
    // obtain updated qubit information
    std::vector<UINT> new_target_index_list, new_target_commutation_list,
        new_control_index_list, new_control_value_list;
    get_new_qubit_list(gate_first, gate_second, new_target_index_list,
        new_target_commutation_list, new_control_index_list,
        new_control_value_list);

    // sort by index
    std::vector<std::pair<UINT, UINT>> new_target_qubit_list,
        new_control_qubit_list;
    for (UINT i = 0; i < new_target_index_list.size(); ++i) {
        new_target_qubit_list.push_back(std::pair<UINT, UINT>(
            new_target_index_list[i], new_target_commutation_list[i]));
    }
    for (UINT i = 0; i < new_control_index_list.size(); ++i) {
        new_control_qubit_list.push_back(std::pair<UINT, UINT>(
            new_control_index_list[i], new_control_value_list[i]));
    }
    std::sort(new_target_qubit_list.begin(), new_target_qubit_list.end());
    std::sort(new_control_qubit_list.begin(), new_control_qubit_list.end());
    for (UINT i = 0; i < new_target_qubit_list.size(); ++i) {
        new_target_index_list[i] = new_target_qubit_list[i].first;
        new_target_commutation_list[i] = new_target_qubit_list[i].second;
    }
    for (UINT i = 0; i < new_control_qubit_list.size(); ++i) {
        new_control_index_list[i] = new_control_qubit_list[i].first;
        new_control_value_list[i] = new_control_qubit_list[i].second;
    }
    // extend gate matrix to whole qubit list
    ComplexMatrix matrix_first, matrix_second;
    get_extended_matrix(gate_first, new_target_index_list,
        new_control_index_list, matrix_first);
    get_extended_matrix(gate_second, new_target_index_list,
        new_control_index_list, matrix_second);

    /*
    ComplexMatrix orgmat1, orgmat2;
    gate_first->get_matrix(orgmat1);
    gate_second->get_matrix(orgmat2);
    std::cout << "first gate is extended from \n"
              << orgmat1 << " \nto\n"
              << matrix_first << "\n\n";
    std::cout << "second gate is extended from \n"
              << orgmat2 << " \nto\n"
              << matrix_second << "\n\n";
    */
    ComplexMatrix new_matrix = matrix_second * matrix_first;

    // generate new matrix gate
    QuantumGateBasic* new_gate = QuantumGateBasic::DenseMatrixGate(
        new_target_index_list, new_matrix, new_target_commutation_list);
    for (UINT i = 0; i < new_control_index_list.size(); ++i) {
        new_gate->add_control_qubit(
            new_control_index_list[i], new_control_value_list[i]);
    }
    new_gate->set_gate_property(
        gate_first->get_property_value() & gate_second->get_property_value());

    // std::cout << "result matrix is " << new_gate << "\n\n";
    return new_gate;
}

// TODO: code is almost common with merge except * or +
DllExport QuantumGateBasic* add(
    const QuantumGateBase* gate_first, const QuantumGateBase* gate_second) {
    // obtain updated qubit information
    std::vector<UINT> new_target_index_list, new_target_commutation_list,
        new_control_index_list, new_control_value_list;
    get_new_qubit_list(gate_first, gate_second, new_target_index_list,
        new_target_commutation_list, new_control_index_list,
        new_control_value_list);

    // sort by index
    std::vector<std::pair<UINT, UINT>> new_target_qubit_list,
        new_control_qubit_list;
    for (UINT i = 0; i < new_target_index_list.size(); ++i) {
        new_target_qubit_list.push_back(std::pair<UINT, UINT>(
            new_target_index_list[i], new_target_commutation_list[i]));
    }
    for (UINT i = 0; i < new_control_index_list.size(); ++i) {
        new_control_qubit_list.push_back(std::pair<UINT, UINT>(
            new_control_index_list[i], new_control_value_list[i]));
    }
    std::sort(new_target_qubit_list.begin(), new_target_qubit_list.end());
    std::sort(new_control_qubit_list.begin(), new_control_qubit_list.end());
    for (UINT i = 0; i < new_target_qubit_list.size(); ++i) {
        new_target_index_list[i] = new_target_qubit_list[i].first;
        new_target_commutation_list[i] = new_target_qubit_list[i].second;
    }
    for (UINT i = 0; i < new_control_qubit_list.size(); ++i) {
        new_control_index_list[i] = new_control_qubit_list[i].first;
        new_control_value_list[i] = new_control_qubit_list[i].second;
    }
    // extend gate matrix to whole qubit list
    ComplexMatrix matrix_first, matrix_second;
    get_extended_matrix(gate_first, new_target_index_list,
        new_control_index_list, matrix_first);
    get_extended_matrix(gate_second, new_target_index_list,
        new_control_index_list, matrix_second);

    /*
    ComplexMatrix orgmat1, orgmat2;
    gate_first->get_matrix(orgmat1);
    gate_second->get_matrix(orgmat2);
    std::cout << "first gate is extended from \n"
              << orgmat1 << " \nto\n"
              << matrix_first << "\n\n";
    std::cout << "second gate is extended from \n"
              << orgmat2 << " \nto\n"
              << matrix_second << "\n\n";
    */
    ComplexMatrix new_matrix = matrix_second + matrix_first;

    // generate new matrix gate
    QuantumGateBasic* new_gate = QuantumGateBasic::DenseMatrixGate(
        new_target_index_list, new_matrix, new_target_commutation_list);
    for (UINT i = 0; i < new_control_index_list.size(); ++i) {
        new_gate->add_control_qubit(
            new_control_index_list[i], new_control_value_list[i]);
    }
    new_gate->set_gate_property(
        gate_first->get_property_value() & gate_second->get_property_value());

    // std::cout << "result matrix is " << new_gate << "\n\n";
    return new_gate;
}

DllExport QuantumGateBasic* Identity(UINT target_qubit) {
    ComplexMatrix mat = ComplexMatrix::Identity(2, 2);
    auto ptr = QuantumGateBasic::DenseMatrixGate({target_qubit}, mat,
        {FLAG_COMMUTE_X | FLAG_COMMUTE_Y | FLAG_COMMUTE_Z});
    ptr->_set_special_func_type(GateI);
    return ptr;
}
DllExport QuantumGateBasic* X(UINT target_qubit) {
    ComplexMatrix mat(2, 2);
    mat << 0, 1, 1, 0;
    auto ptr = QuantumGateBasic::DenseMatrixGate(
        {target_qubit}, mat, {FLAG_COMMUTE_X});
    ptr->_set_special_func_type(GateX);
    return ptr;
}
DllExport QuantumGateBasic* Y(UINT target_qubit) {
    ComplexMatrix mat(2, 2);
    mat << 0, -1.i, 1.i, 0;
    auto ptr = QuantumGateBasic::DenseMatrixGate(
        {target_qubit}, mat, {FLAG_COMMUTE_Y});
    ptr->_set_special_func_type(GateY);
    return ptr;
}
DllExport QuantumGateBasic* Z(UINT target_qubit) {
    ComplexMatrix mat(2, 2);
    mat << 1, 0, 0, -1;
    auto ptr = QuantumGateBasic::DenseMatrixGate(
        {target_qubit}, mat, {FLAG_COMMUTE_Z});
    ptr->_set_special_func_type(GateZ);
    return ptr;
}
DllExport QuantumGateBasic* sqrtX(UINT target_qubit) {
    ComplexMatrix mat(2, 2);
    mat << 0.5 + 0.5i, 0.5 - 0.5i, 0.5 - 0.5i, 0.5 + 0.5i;
    auto ptr = QuantumGateBasic::DenseMatrixGate(
        {target_qubit}, mat, {FLAG_COMMUTE_X});
    ptr->_set_special_func_type(GateSqrtX);
    return ptr;
}
DllExport QuantumGateBasic* sqrtXdag(UINT target_qubit) {
    ComplexMatrix mat(2, 2);
    mat << 0.5 + 0.5i, 0.5 - 0.5i, 0.5 - 0.5i, 0.5 + 0.5i;
    auto ptr = QuantumGateBasic::DenseMatrixGate(
        {target_qubit}, mat.adjoint(), {FLAG_COMMUTE_X});
    ptr->_set_special_func_type(GateSqrtXdag);
    return ptr;
}
DllExport QuantumGateBasic* sqrtY(UINT target_qubit) {
    ComplexMatrix mat(2, 2);
    mat << 0.5 + 0.5i, -0.5 - 0.5i, 0.5 + 0.5i, 0.5 + 0.5i;
    auto ptr = QuantumGateBasic::DenseMatrixGate(
        {target_qubit}, mat, {FLAG_COMMUTE_Y});
    ptr->_set_special_func_type(GateSqrtY);
    return ptr;
}
DllExport QuantumGateBasic* sqrtYdag(UINT target_qubit) {
    ComplexMatrix mat(2, 2);
    mat << 0.5 + 0.5i, -0.5 - 0.5i, 0.5 + 0.5i, 0.5 + 0.5i;
    auto ptr = QuantumGateBasic::DenseMatrixGate(
        {target_qubit}, mat.adjoint(), {FLAG_COMMUTE_Y});
    ptr->_set_special_func_type(GateSqrtYdag);
    return ptr;
}
DllExport QuantumGateBasic* RX(UINT target_qubit, double rotation_angle) {
    auto ptr = QuantumGateBasic::PauliMatrixGate(
        {target_qubit}, {PAULI_ID_X}, rotation_angle);
    ptr->_set_special_func_type(GateRX);
    return ptr;
}
DllExport QuantumGateBasic* RY(UINT target_qubit, double rotation_angle) {
    auto ptr = QuantumGateBasic::PauliMatrixGate(
        {target_qubit}, {PAULI_ID_Y}, rotation_angle);
    ptr->_set_special_func_type(GateRY);
    return ptr;
}
DllExport QuantumGateBasic* RZ(UINT target_qubit, double rotation_angle) {
    auto ptr = QuantumGateBasic::PauliMatrixGate(
        {target_qubit}, {PAULI_ID_Z}, rotation_angle);
    ptr->_set_special_func_type(GateRZ);
    return ptr;
}
DllExport QuantumGateBasic* H(UINT target_qubit) {
    double invsqrt2 = 1. / sqrt(2.);
    ComplexMatrix mat(2, 2);
    mat << invsqrt2, invsqrt2, invsqrt2, -invsqrt2;
    auto ptr = QuantumGateBasic::DenseMatrixGate({target_qubit}, mat, {});
    ptr->_set_special_func_type(GateH);
    return ptr;
}
DllExport QuantumGateBasic* S(UINT target_qubit) {
    ComplexMatrix mat(2, 2);
    mat << 1., 0., 0., 1.i;
    auto ptr = QuantumGateBasic::DenseMatrixGate(
        {target_qubit}, mat, {FLAG_COMMUTE_Z});
    ptr->_set_special_func_type(GateS);
    return ptr;
}
DllExport QuantumGateBasic* HS(UINT target_qubit) {
    double invsqrt2 = 1. / sqrt(2.);
    ComplexMatrix mat(2, 2);
    mat << invsqrt2, invsqrt2 * 1.i, invsqrt2, -invsqrt2 * 1.i;
    auto ptr = QuantumGateBasic::DenseMatrixGate({target_qubit}, mat, {});
    return ptr;
}
DllExport QuantumGateBasic* Sdag(UINT target_qubit) {
    ComplexMatrix mat(2, 2);
    mat << 1., 0., 0., -1.i;
    auto ptr = QuantumGateBasic::DenseMatrixGate(
        {target_qubit}, mat, {FLAG_COMMUTE_Z});
    ptr->_set_special_func_type(GateSdag);
    return ptr;
}
DllExport QuantumGateBasic* T(UINT target_qubit) {
    ComplexMatrix mat(2, 2);
    mat << 1., 0., 0., (1. + 1.i) / sqrt(2.);
    auto ptr = QuantumGateBasic::DenseMatrixGate(
        {target_qubit}, mat, {FLAG_COMMUTE_Z});
    ptr->_set_special_func_type(GateT);
    return ptr;
}
DllExport QuantumGateBasic* Tdag(UINT target_qubit) {
    ComplexMatrix mat(2, 2);
    mat << 1., 0., 0., (1. - 1.i) / sqrt(2.);
    auto ptr = QuantumGateBasic::DenseMatrixGate(
        {target_qubit}, mat, {FLAG_COMMUTE_Z});
    ptr->_set_special_func_type(GateTdag);
    return ptr;
}
DllExport QuantumGateBasic* P0(UINT target_qubit) {
    ComplexMatrix mat(2, 2);
    mat << 1., 0., 0., 0;
    auto ptr = QuantumGateBasic::DenseMatrixGate(
        {target_qubit}, mat, {FLAG_COMMUTE_Z});
    ptr->_set_special_func_type(GateP0);
    return ptr;
}
DllExport QuantumGateBasic* P1(UINT target_qubit) {
    ComplexMatrix mat(2, 2);
    mat << 0., 0., 0., 1;
    auto ptr = QuantumGateBasic::DenseMatrixGate(
        {target_qubit}, mat, {FLAG_COMMUTE_Z});
    ptr->_set_special_func_type(GateP1);
    return ptr;
}
DllExport QuantumGateBasic* CX(UINT control_qubit, UINT target_qubit) {
    ComplexMatrix mat;
    get_Pauli_matrix(mat, {PAULI_ID_X});
    auto ptr = QuantumGateBasic::DenseMatrixGate(
        {target_qubit}, mat, {FLAG_COMMUTE_X});
    ptr->add_control_qubit(control_qubit, 1);
    ptr->_set_special_func_type(GateCX);
    return ptr;
}
DllExport QuantumGateBasic* CNOT(UINT control_qubit, UINT target_qubit) {
    return CX(control_qubit, target_qubit);
}
DllExport QuantumGateBasic* CY(UINT control_qubit, UINT target_qubit) {
    ComplexMatrix mat;
    get_Pauli_matrix(mat, {PAULI_ID_Y});
    auto ptr = QuantumGateBasic::DenseMatrixGate(
        {target_qubit}, mat, {FLAG_COMMUTE_Y});
    ptr->add_control_qubit(control_qubit, 1);
    ptr->_set_special_func_type(GateCY);
    return ptr;
}
DllExport QuantumGateBasic* CZ(UINT control_qubit, UINT target_qubit) {
    ComplexMatrix mat;
    get_Pauli_matrix(mat, {PAULI_ID_Z});
    auto ptr = QuantumGateBasic::DenseMatrixGate(
        {target_qubit}, mat, {FLAG_COMMUTE_Z});
    ptr->add_control_qubit(control_qubit, 1);
    ptr->_set_special_func_type(GateCZ);
    return ptr;
}
DllExport QuantumGateBasic* SWAP(UINT target_qubit1, UINT target_qubit2) {
    ComplexMatrix mat(4, 4);
    mat << 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1;
    auto ptr = QuantumGateBasic::DenseMatrixGate(
        {target_qubit1, target_qubit2}, mat, {});
    ptr->_set_special_func_type(GateSWAP);
    return ptr;
}
DllExport QuantumGateBasic* ISWAP(UINT target_qubit1, UINT target_qubit2) {
    ComplexMatrix mat(4, 4);
    mat << 1, 0, 0, 0, 0, 0, 1.i, 0, 0, 1.i, 0, 0, 0, 0, 0, 1;
    auto ptr = QuantumGateBasic::DenseMatrixGate(
        {target_qubit1, target_qubit2}, mat, {});
    return ptr;
}
DllExport QuantumGateBasic* Toffoli(
    UINT control_qubit1, UINT control_qubit2, UINT target_qubit) {
    ComplexMatrix mat;
    get_Pauli_matrix(mat, {PAULI_ID_X});
    auto ptr = QuantumGateBasic::DenseMatrixGate(
        {target_qubit}, mat, {FLAG_COMMUTE_X});
    ptr->add_control_qubit(control_qubit1, 1);
    ptr->add_control_qubit(control_qubit2, 1);
    return ptr;
}
DllExport QuantumGateBasic* Fredkin(
    UINT control_qubit, UINT target_qubit1, UINT target_qubit2) {
    ComplexMatrix mat(4, 4);
    mat << 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1;
    auto ptr = QuantumGateBasic::DenseMatrixGate(
        {target_qubit1, target_qubit2}, mat, {});
    ptr->add_control_qubit(control_qubit, 1);
    return ptr;
}

DllExport QuantumGateBasic* DenseMatrix(
    UINT target_index, ComplexMatrix matrix) {
    std::vector<UINT> target_list(1, target_index);
    return QuantumGateBasic::DenseMatrixGate(target_list, matrix);
}
DllExport QuantumGateBasic* DenseMatrix(
    std::vector<UINT> target_list, ComplexMatrix matrix) {
    if (!check_is_unique_index_list(target_list)) {
        throw std::invalid_argument("target list contains duplicated values.");
    }
    return QuantumGateBasic::DenseMatrixGate(target_list, matrix);
}
DllExport QuantumGateBasic* SparseMatrix(
    UINT target_index, SparseComplexMatrix matrix) {
    std::vector<UINT> target_list(1, target_index);
    return QuantumGateBasic::SparseMatrixGate(target_list, matrix);
}
DllExport QuantumGateBasic* SparseMatrix(
    std::vector<UINT> target_list, SparseComplexMatrix matrix) {
    if (!check_is_unique_index_list(target_list)) {
        throw std::invalid_argument("target list contains duplicated values.");
    }
    return QuantumGateBasic::SparseMatrixGate(target_list, matrix);
}
}  // namespace gate
