#pragma once

#include "gate.hpp"
#include "state.hpp"
#ifndef _MSC_VER
extern "C" {
#include <csim/update_ops.h>
}
#else
#include <csim/update_ops.h>
#endif
#include <cmath>
#include <csim/update_ops_cpp.hpp>

#ifdef _USE_GPU
#include <gpusim/update_ops_cuda.h>
#endif
#include <iostream>

/**
 * \~japanese-en �t�ÓT��H�̂�\���N���X
 */
class ClsReversibleBooleanGate : public QuantumGateBase {
private:
    std::function<ITYPE(ITYPE, ITYPE)> function_ptr;

public:
    ClsReversibleBooleanGate(std::vector<UINT> target_qubit_index_list,
        std::function<ITYPE(ITYPE, ITYPE)> _function_ptr)
        : function_ptr(_function_ptr) {
        for (auto val : target_qubit_index_list) {
            this->_target_qubit_list.push_back(TargetQubitInfo(val, 0));
        }
        this->_name = "ReversibleBoolean";
    };

    /**
     * \~japanese-en �ʎq��Ԃ��X�V����
     *
     * @param state �X�V����ʎq���
     */
    virtual void update_quantum_state(QuantumStateBase* state) override {
        std::vector<UINT> target_index;
        std::transform(this->_target_qubit_list.cbegin(),
            this->_target_qubit_list.cend(), std::back_inserter(target_index),
            [](auto value) { return value.index(); });
        if (state->is_state_vector()) {
#ifdef _USE_GPU
            if (state->get_device_name() == "gpu") {
                std::cerr << "Not Implemented" << std::endl;
                exit(0);
                // reversible_boolean_gate_gpu(target_index.data(),
                // target_index.size(), function_ptr, state->data_c(),
                // state->dim);
            } else {
                reversible_boolean_gate(target_index.data(),
                    (UINT)target_index.size(), function_ptr, state->data_c(),
                    state->dim);
            }
#else
            reversible_boolean_gate(target_index.data(),
                (UINT)target_index.size(), function_ptr, state->data_c(),
                state->dim);
#endif
        } else {
            std::cerr << "not implemented" << std::endl;
        }
    };
    /**
     * \~japanese-en
     * ���g�̃f�B�[�v�R�s�[�𐶐�����
     *
     * @return ���g�̃f�B�[�v�R�s�[
     */
    virtual QuantumGateBase* copy() const override {
        return new ClsReversibleBooleanGate(*this);
    };
    /**
     * \~japanese-en ���g�̃Q�[�g�s����Z�b�g����
     *
     * @param matrix �s����Z�b�g����ϐ��̎Q��
     */
    virtual void set_matrix(ComplexMatrix& matrix) const override {
        ITYPE matrix_dim = 1ULL << this->_target_qubit_list.size();
        matrix = ComplexMatrix::Zero(matrix_dim, matrix_dim);
        for (ITYPE index = 0; index < matrix_dim; ++index) {
            ITYPE target_index = function_ptr(index, matrix_dim);
            matrix(target_index, index) = 1;
        }
    }
};
