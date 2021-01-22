#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include <boost/dynamic_bitset.hpp>
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "type.hpp"
#include "utility.hpp"

#ifdef _USE_GPU
#include <gpusim/stat_ops.h>
#endif

#ifndef _MSC_VER
extern "C" {
#include <csim/stat_ops.h>
#include <csim/stat_ops_dm.h>
}
#else
#include <csim/stat_ops.h>
#include <csim/stat_ops_dm.h>
#endif
#include "pauli_operator.hpp"
#include "state.hpp"

PauliOperator::PauliOperator(std::string strings, CPPCTYPE coef) {
    _coef = coef;
    std::stringstream ss(strings);
    std::string pauli_str;
    UINT index, pauli_type = 0;
    UINT max_index = 0;
    while (!ss.eof()) {
        ss >> pauli_str >> index;
        max_index = std::max(max_index, index);
        if (pauli_str.length() == 0) break;
        if (pauli_str == "I" || pauli_str == "i")
            pauli_type = 0;
        else if (pauli_str == "X" || pauli_str == "x")
            pauli_type = 1;
        else if (pauli_str == "Y" || pauli_str == "y")
            pauli_type = 2;
        else if (pauli_str == "Z" || pauli_str == "z")
            pauli_type = 3;
        else {
            fprintf(stderr, "invalid Pauli string is given : %s\n ",
                pauli_str.c_str());
            assert(false);
        }
        if (pauli_type != 0) this->add_single_Pauli(index, pauli_type);
    }
    this->set_bits();
}

PauliOperator::PauliOperator(const std::vector<UINT>& target_qubit_list,
    std::string Pauli_operator_type_list, CPPCTYPE coef) {
    _coef = coef;
    UINT term_count = (UINT)(strlen(Pauli_operator_type_list.c_str()));
    UINT pauli_type = 0;
    for (UINT term_index = 0; term_index < term_count; ++term_index) {
        if (Pauli_operator_type_list[term_index] == 'i' ||
            Pauli_operator_type_list[term_index] == 'I') {
            pauli_type = 0;
        } else if (Pauli_operator_type_list[term_index] == 'x' ||
                   Pauli_operator_type_list[term_index] == 'X') {
            pauli_type = 1;
        } else if (Pauli_operator_type_list[term_index] == 'y' ||
                   Pauli_operator_type_list[term_index] == 'Y') {
            pauli_type = 2;
        } else if (Pauli_operator_type_list[term_index] == 'z' ||
                   Pauli_operator_type_list[term_index] == 'Z') {
            pauli_type = 3;
        } else {
            fprintf(stderr, "invalid Pauli string is given\n");
            assert(false);
        }

        if (pauli_type != 0)
            this->add_single_Pauli(target_qubit_list[term_index], pauli_type);
    }
    this->set_bits();
}

PauliOperator::PauliOperator(
    const std::vector<UINT>& pauli_list, CPPCTYPE coef) {
    _coef = coef;
    for (UINT term_index = 0; term_index < pauli_list.size(); ++term_index) {
        if (pauli_list[term_index] != 0)
            this->add_single_Pauli(term_index, pauli_list[term_index]);
    }
    this->set_bits();
}

PauliOperator::PauliOperator(const std::vector<UINT>& target_qubit_index_list,
    const std::vector<UINT>& target_qubit_pauli_list, CPPCTYPE coef) {
    _coef = coef;
    assert(target_qubit_index_list.size() == target_qubit_pauli_list.size());
    for (UINT term_index = 0; term_index < target_qubit_index_list.size();
         ++term_index) {
        this->add_single_Pauli(target_qubit_index_list[term_index],
            target_qubit_pauli_list[term_index]);
    }
    this->set_bits();
}

PauliOperator::PauliOperator(const boost::dynamic_bitset<>& x,
    const boost::dynamic_bitset<>& z, CPPCTYPE coef = 1.) {
    _x = x;
    _z = z;
    _coef = coef;
    for (UINT i = 0; i < _x.size(); i++) {
        UINT pauli_type = 0;
        if (_x[i] && !_z[i]) {
            pauli_type = 1;
        } else if (_x[i] && _z[i]) {
            pauli_type = 2;
        } else if (!_x[i] && _z[i]) {
            pauli_type = 3;
        }
        if (pauli_type != 0) {
            this->add_single_Pauli(i, pauli_type);
        }
    }
}

void PauliOperator::add_single_Pauli(UINT qubit_index, UINT pauli_type) {
    this->_pauli_list.push_back(SinglePauliOperator(qubit_index, pauli_type));
    if (qubit_index > _x.size()) {
        _x.resize(qubit_index + 1);
        _z.resize(qubit_index + 1);
    }
    if (pauli_type == 1) {
        _x.set(qubit_index);
    } else if (pauli_type == 2) {
        _x.set(qubit_index);
        _z.set(qubit_index);
    } else if (pauli_type == 3) {
        _z.set(qubit_index);
    }
}

void PauliOperator::set_bits() {
    UINT max_index = 0;
    for (int i = 0; i < _pauli_list.size(); i++) {
        max_index = std::max(max_index, _pauli_list[i].index());
    }
    _x.resize(max_index + 1);
    _z.resize(max_index + 1);
    for (int i = 0; i < _pauli_list.size(); i++) {
        if (_pauli_list[i].pauli_id() == 1) {
            _x.set(_pauli_list[i].index());
        } else if (_pauli_list[i].pauli_id() == 2) {
            _x.set(_pauli_list[i].index());
            _z.set(_pauli_list[i].index());
        } else if (_pauli_list[i].pauli_id() == 3) {
            _z.set(_pauli_list[i].index());
        }
    }
}

CPPCTYPE PauliOperator::get_expectation_value(
    const QuantumStateBase* state) const {
    if (state->is_state_vector()) {
#ifdef _USE_GPU
        if (state->get_device_name() == "gpu") {
            return _coef *
                   expectation_value_multi_qubit_Pauli_operator_partial_list_host(
                       this->get_index_list().data(),
                       this->get_pauli_id_list().data(),
                       (UINT)this->get_index_list().size(), state->data(),
                       state->dim, state->get_cuda_stream(),
                       state->device_number);
        } else {
            return _coef *
                   expectation_value_multi_qubit_Pauli_operator_partial_list(
                       this->get_index_list().data(),
                       this->get_pauli_id_list().data(),
                       (UINT)this->get_index_list().size(), state->data_c(),
                       state->dim);
        }
#else
        return _coef *
               expectation_value_multi_qubit_Pauli_operator_partial_list(
                   this->get_index_list().data(),
                   this->get_pauli_id_list().data(),
                   (UINT)this->get_index_list().size(), state->data_c(),
                   state->dim);
#endif
    } else {
        return _coef *
               dm_expectation_value_multi_qubit_Pauli_operator_partial_list(
                   this->get_index_list().data(),
                   this->get_pauli_id_list().data(),
                   (UINT)this->get_index_list().size(), state->data_c(),
                   state->dim);
    }
}

CPPCTYPE PauliOperator::get_transition_amplitude(
    const QuantumStateBase* state_bra,
    const QuantumStateBase* state_ket) const {
    if ((!state_bra->is_state_vector()) || (!state_ket->is_state_vector())) {
        std::cerr
            << "get_transition_amplitude for density matrix is not implemented"
            << std::endl;
    }
#ifdef _USE_GPU
    if (state_ket->get_device_name() == "gpu" &&
        state_bra->get_device_name() == "gpu") {
        return _coef *
               (CPPCTYPE)
                   transition_amplitude_multi_qubit_Pauli_operator_partial_list_host(
                       this->get_index_list().data(),
                       this->get_pauli_id_list().data(),
                       (UINT)this->get_index_list().size(), state_bra->data(),
                       state_ket->data(), state_bra->dim,
                       state_ket->get_cuda_stream(), state_ket->device_number);
    } else {
        return _coef *
               (CPPCTYPE)
                   transition_amplitude_multi_qubit_Pauli_operator_partial_list(
                       this->get_index_list().data(),
                       this->get_pauli_id_list().data(),
                       (UINT)this->get_index_list().size(), state_bra->data_c(),
                       state_ket->data_c(), state_bra->dim);
    }
#else
    return _coef *
           (CPPCTYPE)
               transition_amplitude_multi_qubit_Pauli_operator_partial_list(
                   this->get_index_list().data(),
                   this->get_pauli_id_list().data(),
                   (UINT)this->get_index_list().size(), state_bra->data_c(),
                   state_ket->data_c(), state_bra->dim);
#endif
}

PauliOperator* PauliOperator::copy() const {
    auto pauli = new PauliOperator(this->_coef);
    for (auto val : this->_pauli_list) {
        pauli->add_single_Pauli(val.index(), val.pauli_id());
    }
    return pauli;
}

PauliOperator PauliOperator::operator*(PauliOperator& target) {
    CPPCTYPE bits_coef = 1.;
    auto target_x = target.get_x_bits();
    auto target_z = target.get_x_bits();
    for (int i = 0; i < _x.size(); i++) {
        if (_x[i] && !_z[i]) {  // X
            if (!target_x[i] && target_z[i]) {
                bits_coef *= -1i;
            } else if (target_x[i] && target_z[i]) {
                bits_coef *= 1i;
            }
        } else if (!_x[i] && _z[i]) {           // Z
            if (target_x[i] && !target_z[i]) {  // X
                bits_coef *= -1i;
            } else if (target_x[i] && target_z[i]) {  // Y
                bits_coef *= 1i;
            }
        } else if (_x[i] && _z[i]) {            // Y
            if (target_x[i] && !target_z[i]) {  // X
                bits_coef *= 1i;
            } else if (!target_x[i] && target_z[i]) {  // Z
                bits_coef *= 1i;
            }
        }
        PauliOperator res(_x ^ target.get_x_bits(), _z ^ target.get_z_bits(),
            _coef * target.get_coef() * bits_coef);
        return res;
    }
