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

/**
 * \~japanese-en 1つの実数パラメータを持った量子ゲートのクラス
 *
 * parametricな量子ゲートを管理するクラス。
 * _parameter_idが0以上のとき、パラメータはQuantumParametricCircuitが管理していて、id
 * -> パラメータの変換表を操作時に渡す(new-style)
 * _parameter_idが-1のとき、パラメータは_angleに直接格納される(old-style)
 */
class QuantumGate_SingleParameter : public QuantumGateBase {
protected:
    ParameterId _parameter_id;
    double _parameter_coef;
    double _angle;

public:
    /**
     * \~japanese-en old-styleのコンストラクタ
     *
     * @param[in] angle パラメータの値
     * @return 生成されたインスタンス
     */
    explicit QuantumGate_SingleParameter(double angle);  // old-style
    /**
     * \~japanese-en new-styleのコンストラクタ
     *
     * @param[in] parameter_id パラメータのId(未割り当ての場合は-1)
     * @param[in] parameter_coef パラメータの倍率
     * @return 生成されたインスタンス
     */
    explicit QuantumGate_SingleParameter(const ParameterId& parameter_id,
        double parameter_coef = 1.);  // new-style

    bool is_old_style() const;  // both-style
    bool is_new_style() const;  // both-style

    ParameterId get_parameter_id() const;                    // both-style
    void set_parameter_id(const ParameterId& parameter_id);  // new-style
    double get_parameter_coef() const;                       // both-style
    void set_parameter_value(double value);                  // old-style
    double get_parameter_value() const;                      // old-style
    double get_angle() const;                                // old-style
    double get_angle(
        const std::vector<double> parameter_list) const;  // new-style
    virtual QuantumGate_SingleParameter* copy()
        const override = 0;  // both-style
    virtual void update_quantum_state(
        QuantumStateBase* state) override = 0;  // old-style
    virtual void update_quantum_state(QuantumStateBase* state,
        const std::vector<double>& parameter_list) = 0;  // new-style
};

/**
 * \~japanese-en
 * 1つの実数パラメータを持ち、1つの量子ビットに作用する量子ゲートのクラス
 */
class QuantumGate_SingleParameterOneQubitRotation
    : public QuantumGate_SingleParameter {
protected:
    typedef void(T_UPDATE_FUNC)(UINT, double, CTYPE*, ITYPE);
    typedef void(T_GPU_UPDATE_FUNC)(UINT, double, void*, ITYPE, void*, UINT);
    T_UPDATE_FUNC* _update_func = nullptr;
    T_UPDATE_FUNC* _update_func_dm = nullptr;
    T_GPU_UPDATE_FUNC* _update_func_gpu = nullptr;

    QuantumGate_SingleParameterOneQubitRotation(double angle);  // old-style
    QuantumGate_SingleParameterOneQubitRotation(const ParameterId& parameter_id,
        double parameter_coef = 1.);  // new-style
    void _update_quantum_state(
        QuantumStateBase* state, double angle);  // both-style

public:
    virtual void update_quantum_state(QuantumStateBase* state);  // old-style
    virtual void update_quantum_state(QuantumStateBase* state,
        const std::vector<double>& parameter_list);  // new-style
};

/**
 * \~japanese-en
 * ParametricRXゲートのクラス
 */
class ClsParametricRXGate : public QuantumGate_SingleParameterOneQubitRotation {
private:
    void _set_matrix(ComplexMatrix& matrix, double angle) const;  // both-style

public:
    ClsParametricRXGate(UINT target_qubit_index, double angle);  // old-style
    ClsParametricRXGate(UINT target_qubit_index,
        const ParameterId& parameter_id,
        double parameter_coef = 1.);  // new-style

    virtual void set_matrix(ComplexMatrix& matrix) const;  // old-style
    virtual void set_matrix(ComplexMatrix& matrix,
        const std::vector<double>& parameter_list) const;        // new-style
    virtual QuantumGate_SingleParameter* copy() const override;  // both-style
};

/**
 * \~japanese-en
 * ParametricRYゲートのクラス
 */
class ClsParametricRYGate : public QuantumGate_SingleParameterOneQubitRotation {
private:
    void _set_matrix(ComplexMatrix& matrix, double angle) const;  // both-style

public:
    ClsParametricRYGate(UINT target_qubit_index, double angle);  // old-style
    ClsParametricRYGate(UINT target_qubit_index,
        const ParameterId& parameter_id,
        double parameter_coef = 1.);  // new-style

    virtual void set_matrix(ComplexMatrix& matrix) const;  // old-style
    virtual void set_matrix(ComplexMatrix& matrix,
        const std::vector<double>& parameter_list) const;        // new-style
    virtual QuantumGate_SingleParameter* copy() const override;  // both-style
};

/**
 * \~japanese-en
 * ParametricRZゲートのクラス
 */
class ClsParametricRZGate : public QuantumGate_SingleParameterOneQubitRotation {
private:
    void _set_matrix(ComplexMatrix& matrix, double angle) const;  // both-style

public:
    ClsParametricRZGate(UINT target_qubit_index, double angle);  // old-style
    ClsParametricRZGate(UINT target_qubit_index,
        const ParameterId& parameter_id,
        double parameter_coef = 1.);  // new-style

    virtual void set_matrix(ComplexMatrix& matrix) const;  // old-style
    virtual void set_matrix(ComplexMatrix& matrix,
        const std::vector<double>& parameter_list) const;        // new-style
    virtual QuantumGate_SingleParameter* copy() const override;  // both-style
};

class ClsParametricPauliRotationGate : public QuantumGate_SingleParameter {
protected:
    PauliOperator* _pauli;
    void _update_quantum_state(
        QuantumStateBase* state, double angle);                   // both-style
    void _set_matrix(ComplexMatrix& matrix, double angle) const;  // both-style

public:
    ClsParametricPauliRotationGate(
        double angle, PauliOperator* pauli);  // old-style
    ClsParametricPauliRotationGate(PauliOperator* pauli,
        const ParameterId& parameter_id,
        double parameter_coef = 1.);                             // new-style
    virtual ~ClsParametricPauliRotationGate();                   // both-style
    virtual void update_quantum_state(QuantumStateBase* state);  // old-style
    virtual void update_quantum_state(QuantumStateBase* state,
        const std::vector<double>& parameter_list);              // new-style
    virtual QuantumGate_SingleParameter* copy() const override;  // both-style
    virtual void set_matrix(ComplexMatrix& matrix) const;        // old-style
    virtual void set_matrix(ComplexMatrix& matrix,
        const std::vector<double>& parameter_list) const;  // new-style
    virtual PauliOperator* get_pauli() const;              // both-style
};
