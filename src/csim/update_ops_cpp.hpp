#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <functional>
#include "type.hpp"

DllExport void double_qubit_dense_matrix_gate(UINT target_qubit_index0,
    UINT target_qubit_index1, const CTYPE matrix[16], CTYPE* state, ITYPE dim);
DllExport void double_qubit_dense_matrix_gate(UINT target_qubit_index0,
    UINT target_qubit_index1, const Eigen::Matrix4cd& eigen_matrix,
    CTYPE* state, ITYPE dim);
DllExport void double_qubit_dense_matrix_gate_eigen(UINT target_qubit_index0,
    UINT target_qubit_index1, const Eigen::Matrix4cd& eigen_matrix,
    CTYPE* state, ITYPE dim);

DllExport void multi_qubit_dense_matrix_gate_eigen(
    const UINT* target_qubit_index_list, UINT target_qubit_index_count,
    const CTYPE* matrix, CTYPE* state, ITYPE dim);
DllExport void multi_qubit_dense_matrix_gate_eigen(
    const UINT* target_qubit_index_list, UINT target_qubit_index_count,
    const Eigen::MatrixXcd& eigen_matrix, CTYPE* state, ITYPE dim);
DllExport void multi_qubit_dense_matrix_gate_eigen(
    const UINT* target_qubit_index_list, UINT target_qubit_index_count,
    const Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic,
        Eigen::RowMajor>& eigen_matrix,
    CTYPE* state, ITYPE dim);

DllExport void multi_qubit_sparse_matrix_gate_eigen(
    const UINT* target_qubit_index_list, UINT target_qubit_index_count,
    const Eigen::SparseMatrix<std::complex<double>>& eigen_matrix, CTYPE* state,
    ITYPE dim);

/**
 * \~english
 * Apply reversible boolean function as a unitary gate.
 *
 * Apply reversible boolean function as a unitary gate. Boolean function is
 * given as a pointer of int -> int function.
 *
 * @param[in] target_qubit_index_list �^�[�Q�b�g�ʎq�r�b�g�̃��X�g
 * @param[in] target_qubit_index_count �^�[�Q�b�g�ʎq�r�b�g�̐�
 * @param[in] matrix �Y��������ёΏۃr�b�g�̎������󂯎��ƓY������Ԃ��֐�
 * @param[in,out] state �ʎq���
 * @param[in] dim ����
 *
 *
 * \~japanese-en
 * �t��H�֐������j�^���Q�[�g�Ƃ��č�p����
 *
 *  �t��H�֐������j�^���Q�[�g�Ƃ��č�p����B�t��H�֐��͓Y������^����ƌ��ʂ̓Y������Ԃ��֐��B
 *
 * @param[in] target_qubit_index_list �^�[�Q�b�g�ʎq�r�b�g�̃��X�g
 * @param[in] target_qubit_index_count �^�[�Q�b�g�ʎq�r�b�g�̐�
 * @param[in] matrix �Y��������ёΏۃr�b�g�̎������󂯎��ƓY������Ԃ��֐�
 * @param[in,out] state �ʎq���
 * @param[in] dim ����
 *
 */
DllExport void reversible_boolean_gate(const UINT* target_qubit_index_list,
    UINT target_qubit_index_count,
    std::function<ITYPE(ITYPE, ITYPE)> function_ptr, CTYPE* state, ITYPE dim);
