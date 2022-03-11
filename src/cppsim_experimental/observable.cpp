#include "observable.hpp"

#include <numeric>

#ifdef _USE_GPU
#include "../gpusim/stat_ops.h"
#endif

#include "pauli_operator.hpp"
#include "state.hpp"
#include "type.hpp"

CPPCTYPE Observable::calc_coef(
    const MultiQubitPauliOperator& a, const MultiQubitPauliOperator& b) const {
    auto x_a = a.get_x_bits();
    auto z_a = a.get_z_bits();
    auto x_b = b.get_x_bits();
    auto z_b = b.get_z_bits();
    const size_t max_size = std::max(x_a.size(), x_b.size());
    if (x_a.size() != x_b.size()) {
        x_a.resize(max_size);
        z_a.resize(max_size);
        x_b.resize(max_size);
        z_b.resize(max_size);
    }
    CPPCTYPE res = 1.0;
    CPPCTYPE I = 1.0i;
    ITYPE i;
#pragma omp parallel
    {
        // 各スレッドごとに独立な変数を用意する

        CPPCTYPE res_private = 1.0;
#pragma omp for nowait
        for (i = 0; i < x_a.size(); i++) {
            if (x_a[i] && !z_a[i]) {            // a = X
                if (!x_b[i] && z_b[i]) {        // b = Z
                    res_private *= -I;          // XZ = -iY
                } else if (x_b[i] && z_b[i]) {  // b = Y
                    res_private *= I;           // XY = iZ
                }
            } else if (!x_a[i] && z_a[i]) {     // a = Z
                if (x_b[i] && !z_b[i]) {        // b = X
                    res_private *= I;           // ZX = iY
                } else if (x_b[i] && z_b[i]) {  // b = Y
                    res_private *= -I;          // ZY = -iX
                }
            } else if (x_a[i] && z_a[i]) {       // a = Y
                if (x_b[i] && !z_b[i]) {         // b = X
                    res_private *= -I;           // YX = -iZ
                } else if (!x_b[i] && z_b[i]) {  // b = Z
                    res_private *= I;            // YZ = iX
                }
            }
        }
#pragma omp critical
        {
            // 各スレッドで計算した結果を結合する
            res *= res_private;
        }
    }
    return res;
}

std::pair<CPPCTYPE, MultiQubitPauliOperator> Observable::get_term(
    const UINT index) const {
    return std::make_pair(
        this->_coef_list.at(index), this->_pauli_terms.at(index));
}

std::unordered_map<std::string, ITYPE> Observable::get_dict() const {
    return _term_dict;
}

void Observable::add_term(const CPPCTYPE coef, MultiQubitPauliOperator op) {
    this->_coef_list.push_back(coef);
    this->_pauli_terms.push_back(op);
    this->_term_dict[op.to_string()] = _coef_list.size() - 1;
}

void Observable::add_term(const CPPCTYPE coef, std::string s) {
    MultiQubitPauliOperator op(s);
    add_term(coef, op);
}

void Observable::add_term(const std::vector<CPPCTYPE> coef_list,
    std::vector<MultiQubitPauliOperator> pauli_terms) {
    int base_size = this->_coef_list.size();
    int changed_size = base_size + coef_list.size();
    // vectorをindexでアクセス出来るようにする
    this->_coef_list.resize(changed_size);
    this->_pauli_terms.resize(changed_size);

#pragma omp parallel
    {
        std::unordered_map<std::string, ITYPE> term_dict_private;
#pragma omp for nowait
        for (ITYPE index = 0; index < coef_list.size(); index++) {
            int insert_pos = base_size + index;
            this->_coef_list[insert_pos] = coef_list[index];
            this->_pauli_terms[insert_pos] = pauli_terms[index];
            // dictには並列で同時に書き込めないので、一旦term_dict_privateに書き込む
            term_dict_private[pauli_terms[index].to_string()] = insert_pos;
        }

#pragma omp critical
        {
            _term_dict.insert(
                term_dict_private.begin(), term_dict_private.end());
        }
    }
}

void Observable::remove_term(UINT index) {
    this->_term_dict.erase(this->_pauli_terms.at(index).to_string());
    this->_coef_list.erase(this->_coef_list.begin() + index);
    this->_pauli_terms.erase(this->_pauli_terms.begin() + index);

    // index番目の項を削除したので、index番目以降の項のindexが1つずれる

    for (ITYPE i = 0; i < this->_coef_list.size() - index; i++) {
        this->_term_dict[this->_pauli_terms.at(index + i).to_string()] =
            index + i;
    }
}

CPPCTYPE Observable::get_expectation_value(
    const QuantumStateBase* state) const {
    CPPCTYPE sum = 0;
    for (ITYPE index = 0; index < this->_pauli_terms.size(); ++index) {
        sum += this->_coef_list.at(index) *
               this->_pauli_terms.at(index).get_expectation_value(state);
    }
    return sum;
}

CPPCTYPE Observable::get_transition_amplitude(const QuantumStateBase* state_bra,
    const QuantumStateBase* state_ket) const {
    CPPCTYPE sum = 0;
    for (ITYPE index = 0; index < this->_pauli_terms.size(); ++index) {
        sum += this->_coef_list.at(index) *
               this->_pauli_terms.at(index).get_transition_amplitude(
                   state_bra, state_ket);
    }
    return sum;
}

Observable* Observable::copy() const {
    Observable* res = new Observable();
    res->add_term(this->_coef_list, this->_pauli_terms);
    return res;
}

Observable Observable::operator+(const Observable& target) const {
    Observable res = *(this->copy());
    res += target;
    return res;
}

Observable& Observable::operator+=(const Observable& target) {
    auto u_map = target.get_dict();

    for (auto item : u_map) {
        if (_term_dict.find(item.first) != _term_dict.end()) {
            ITYPE id = _term_dict[item.first];
            this->_coef_list[id] += target.get_term(item.second).first;
        } else {
            auto term = target.get_term(item.second);
            this->add_term(term.first, term.second);
        }
    }
    return *this;
}

Observable Observable::operator-(const Observable& target) const {
    Observable res = *(this->copy());
    res -= target;
    return res;
}

Observable& Observable::operator-=(const Observable& target) {
    auto u_map = target.get_dict();

    for (auto item : u_map) {
        if (_term_dict.find(item.first) != _term_dict.end()) {
            ITYPE id = _term_dict[item.first];
            this->_coef_list[id] -= target.get_term(item.second).first;
        } else {
            auto term = target.get_term(item.second);
            this->add_term(-term.first, term.second);
        }
    }
    return *this;
}

Observable Observable::operator*(const Observable& target) const {
    Observable res;
    ITYPE i, j;
    for (i = 0; i < this->_pauli_terms.size(); i++) {
        for (j = 0; j < target.get_term_count(); j++) {
            Observable tmp;
            auto term = target.get_term(j);
            CPPCTYPE bits_coef = calc_coef(this->_pauli_terms[i], term.second);
            tmp.add_term(this->_coef_list[i] * term.first * bits_coef,
                this->_pauli_terms[i] * term.second);
            res += tmp;
        }
    }
    return res;
}

Observable Observable::operator*(const CPPCTYPE& target) const {
    Observable res = *(this->copy());
    res *= target;
    return res;
}

Observable& Observable::operator*=(const Observable& target) {
    Observable tmp = (*this) * target;
    this->_coef_list.clear();
    this->_pauli_terms.clear();
    this->_term_dict.clear();
    ITYPE i;
    for (i = 0; i < tmp.get_term_count(); i++) {
        auto term = tmp.get_term(i);
        this->add_term(term.first, term.second);
    }
    return *this;
}

Observable& Observable::operator*=(const CPPCTYPE& target) {
    ITYPE i;
#pragma omp parallel for
    for (i = 0; i < this->_coef_list.size(); i++) {
        this->_coef_list[i] *= target;
    }
    return *this;
}

std::string Observable::to_string() const {
    std::ostringstream ss;
    std::string res;
    ITYPE i;
    for (i = 0; i < get_term_count(); i++) {
        // (1.0-2.0j)
        ss << "(" << _coef_list[i].real();
        // 虚数部には符号をつける
        // +0j or -0j に対応させるためstd::showposを用いる
        ss << std::showpos << _coef_list[i].imag() << "j) ";

        // [X 0 Y 1 Z 2]
        ss << "[" << _pauli_terms[i].to_string() << "]";
        if (i != get_term_count() - 1) {
            ss << " +" << std::endl;
        }
        res += ss.str();
        ss.str("");
        ss.clear();
        // noshowposして、1項目以降の実数部に符号が付かないようにする
        ss << std::noshowpos;
    }
    return res;
}
