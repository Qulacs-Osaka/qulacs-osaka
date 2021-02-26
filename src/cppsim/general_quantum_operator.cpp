#include <cstring>
#include <fstream>

#include "type.hpp"
#include "utility.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif
#ifndef _MSC_VER
extern "C" {
#include <csim/stat_ops.h>
}
#else
#include <csim/stat_ops.h>
#endif
#include <Eigen/Dense>

#include "gate_factory.hpp"
#include "general_quantum_operator.hpp"
#include "pauli_operator.hpp"
#include "state.hpp"

GeneralQuantumOperator::GeneralQuantumOperator(const UINT qubit_count)
    : _qubit_count(qubit_count), _is_hermitian(true) {}

GeneralQuantumOperator::~GeneralQuantumOperator() {
    for (auto& term : this->_operator_list) {
        delete term;
    }
}

void GeneralQuantumOperator::add_operator(const PauliOperator* mpt) {
    PauliOperator* _mpt = mpt->copy();
    if (!check_Pauli_operator(this, _mpt)) {
        std::cerr << "Error: GeneralQuantumOperator::add_operator(const "
                     "PauliOperator*): pauli_operator applies target qubit of "
                     "which the index is larger than qubit_count"
                  << std::endl;
        return;
    }
    if (this->_is_hermitian && std::abs(_mpt->get_coef().imag()) > 0) {
        this->_is_hermitian = false;
    }
    this->_operator_list.push_back(_mpt);
}

void GeneralQuantumOperator::add_operator(
    CPPCTYPE coef, std::string pauli_string) {
    PauliOperator* _mpt = new PauliOperator(pauli_string, coef);
    if (!check_Pauli_operator(this, _mpt)) {
        std::cerr << "Error: "
                     "GeneralQuantumOperator::add_operator(double,std::string):"
                     " pauli_operator applies target qubit of which the index "
                     "is larger than qubit_count"
                  << std::endl;
        return;
    }
    if (this->_is_hermitian && std::abs(coef.imag()) > 0) {
        this->_is_hermitian = false;
    }
    this->add_operator(_mpt);
    delete _mpt;
}

CPPCTYPE GeneralQuantumOperator::get_expectation_value(
    const QuantumStateBase* state) const {
    if (this->_qubit_count != state->qubit_count) {
        std::cerr
            << "Error: GeneralQuantumOperator::get_expectation_value(const "
               "QuantumStateBase*): invalid qubit count"
            << std::endl;
        return 0.;
    }
    CPPCTYPE sum = 0;
    for (auto pauli : this->_operator_list) {
        sum += pauli->get_expectation_value(state);
    }
    return sum;
}

CPPCTYPE GeneralQuantumOperator::get_transition_amplitude(
    const QuantumStateBase* state_bra,
    const QuantumStateBase* state_ket) const {
    if (this->_qubit_count != state_bra->qubit_count ||
        this->_qubit_count != state_ket->qubit_count) {
        std::cerr
            << "Error: GeneralQuantumOperator::get_transition_amplitude(const "
               "QuantumStateBase*, const QuantumStateBase*): invalid qubit "
               "count"
            << std::endl;
        return 0.;
    }

    CPPCTYPE sum = 0;
    for (auto pauli : this->_operator_list) {
        sum += pauli->get_transition_amplitude(state_bra, state_ket);
    }
    return sum;
}

CPPCTYPE
GeneralQuantumOperator::solve_ground_state_eigenvalue_by_arnoldi_method(
    QuantumStateBase* state, const UINT iter_count) const {
    // Implemented based on
    // https://files.transtutors.com/cdn/uploadassignments/472339_1_-numerical-linear-aljebra.pdf
    const auto qubit_count = this->get_qubit_count();
    auto present_state = QuantumState(qubit_count);
    auto tmp_state = QuantumState(qubit_count);
    auto multiplied_state = QuantumState(qubit_count);

    // Vectors composing Krylov subspace.
    std::vector<QuantumStateBase*> state_list;
    state->normalize(state->get_squared_norm());
    state_list.push_back(state);

    ComplexMatrix hessenberg_matrix =
        ComplexMatrix::Zero(iter_count, iter_count);
    for (UINT i = 0; i < iter_count; i++) {
        this->apply_to_state(state_list[i], &multiplied_state);

        for (UINT j = 0; j < i + 1; j++) {
            const auto coef = state::inner_product(
                static_cast<QuantumState*>(state_list[j]), &multiplied_state);
            hessenberg_matrix(j, i) = coef;
            tmp_state.load(state_list[j]);
            tmp_state.multiply_coef(-coef);
            multiplied_state.add_state(&tmp_state);
        }

        const auto norm = multiplied_state.get_squared_norm();
        if (i != iter_count - 1) {
            hessenberg_matrix(i + 1, i) = std::sqrt(norm);
        }
        multiplied_state.normalize(norm);
        state_list.push_back(multiplied_state.copy());
    }

    Eigen::ComplexEigenSolver<ComplexMatrix> eigen_solver(hessenberg_matrix);
    const auto eigenvalues = eigen_solver.eigenvalues();
    const auto eigenvectors = eigen_solver.eigenvectors();

    // Find ground state vector.
    UINT minimum_eigenvalue_index = 0;
    auto minimum_eigenvalue = eigenvalues[0];
    for (UINT i = 0; i < eigenvalues.size(); i++) {
        if (eigenvalues[i].real() < minimum_eigenvalue.real()) {
            minimum_eigenvalue_index = i;
            minimum_eigenvalue = eigenvalues[i];
        }
    }

    // Compose ground state vector.
    present_state.multiply_coef(0.0);
    for (UINT i = 0; i < state_list.size(); i++) {
        tmp_state.load(state_list[i]);
        tmp_state.multiply_coef(eigenvectors(i, minimum_eigenvalue_index));
        present_state.add_state(&tmp_state);
    }
    state->load(&present_state);

    return minimum_eigenvalue;
}

CPPCTYPE GeneralQuantumOperator::solve_ground_state_eigenvalue_by_power_method(
    QuantumStateBase* state, const UINT iter_count, const CPPCTYPE mu) const {
    CPPCTYPE mu_;
    if (mu == 0.0) {
        // mu is not changed from default value.
        mu_ = this->calculate_default_mu();
    } else {
        mu_ = mu;
    }

    // Stores a result of A|a>
    auto multiplied_state = QuantumState(state->qubit_count);
    // Stores a result of -\mu|a>
    auto mu_timed_state = QuantumState(state->qubit_count);
    for (UINT i = 0; i < iter_count; i++) {
        mu_timed_state.load(state);
        mu_timed_state.multiply_coef(-mu_);

        multiplied_state.multiply_coef(0.0);
        this->apply_to_state(state, &multiplied_state);
        state->load(&multiplied_state);
        state->add_state(&mu_timed_state);
        state->normalize(state->get_squared_norm());
    }
    return this->get_expectation_value(state) + mu;
}

GeneralQuantumOperator* GeneralQuantumOperator::copy() const {
    auto quantum_operator = new GeneralQuantumOperator(_qubit_count);
    for (auto pauli : this->_operator_list) {
        quantum_operator->add_operator(pauli->copy());
    }
    return quantum_operator;
}

GeneralQuantumOperator GeneralQuantumOperator::operator+(
    const GeneralQuantumOperator& target) const {
    auto res = this->copy();
    *res += target;
    return *res;
}

GeneralQuantumOperator GeneralQuantumOperator::operator+(
    const PauliOperator& target) const {
    auto res = this->copy();
    *res += target;
    return *res;
}

GeneralQuantumOperator& GeneralQuantumOperator::operator+=(
    const GeneralQuantumOperator& target) {
#pragma omp parallel for
    for (UINT i = 0; i < _operator_list.size(); i++) {
        auto pauli_operator = _operator_list[i];
        for (UINT j = 0; j < target.get_terms().size(); j++) {
            auto target_operator = target.get_terms()[j];
            if (pauli_operator->get_x_bits() == target_operator->get_x_bits() &&
                pauli_operator->get_z_bits() == target_operator->get_z_bits()) {
                _operator_list[i]->change_coef(_operator_list[i]->get_coef() +
                                               target_operator->get_coef());
            }
        }
    }
    for (UINT i = 0; i < _operator_list.size(); i++) {
        auto pauli_operator = _operator_list[i];
        for (UINT j = 0; j < target.get_terms().size(); j++) {
            auto target_operator = target.get_terms()[j];
            if (pauli_operator->get_x_bits() == target_operator->get_x_bits() &&
                pauli_operator->get_z_bits() == target_operator->get_z_bits()) {
                continue;
            }
            this->add_operator(target_operator->copy());
        }
    }
    return *this;
}

GeneralQuantumOperator& GeneralQuantumOperator::operator+=(
    const PauliOperator& target) {
    bool flag = true;
#pragma omp parallel for
    for (UINT i = 0; i < _operator_list.size(); i++) {
        auto pauli_operator = _operator_list[i];
        if (pauli_operator->get_x_bits() == target.get_x_bits() &&
            pauli_operator->get_z_bits() == target.get_z_bits()) {
            _operator_list[i]->change_coef(
                _operator_list[i]->get_coef() + target.get_coef());
            flag = false;
        }
    }
    if (flag) {
        this->add_operator(target.copy());
    }
    return *this;
}

GeneralQuantumOperator GeneralQuantumOperator::operator-(
    const GeneralQuantumOperator& target) const {
    auto res = this->copy();
    *res -= target;
    return *res;
}

GeneralQuantumOperator GeneralQuantumOperator::operator-(
    const PauliOperator& target) const {
    auto res = this->copy();
    *res -= target;
    return *res;
}

GeneralQuantumOperator& GeneralQuantumOperator::operator-=(
    const GeneralQuantumOperator& target) {
#pragma omp parallel for
    for (UINT i = 0; i < _operator_list.size(); i++) {
        auto pauli_operator = _operator_list[i];
        for (UINT j = 0; j < target.get_terms().size(); j++) {
            auto target_operator = target.get_terms()[j];
            if (pauli_operator->get_x_bits() == target_operator->get_x_bits() &&
                pauli_operator->get_z_bits() == target_operator->get_z_bits()) {
                _operator_list[i]->change_coef(_operator_list[i]->get_coef() -
                                               target_operator->get_coef());
            }
        }
    }
    for (UINT i = 0; i < _operator_list.size(); i++) {
        auto pauli_operator = _operator_list[i];
        for (UINT j = 0; j < target.get_terms().size(); j++) {
            auto target_operator = target.get_terms()[j];
            if (pauli_operator->get_x_bits() == target_operator->get_x_bits() &&
                pauli_operator->get_z_bits() == target_operator->get_z_bits()) {
                continue;
            }
            auto copy = target_operator->copy();
            copy->change_coef(-copy->get_coef());
            this->add_operator(copy);
        }
    }
    return *this;
}

GeneralQuantumOperator& GeneralQuantumOperator::operator-=(
    const PauliOperator& target) {
    bool flag = true;
    for (UINT i = 0; i < _operator_list.size(); i++) {
        auto pauli_operator = _operator_list[i];
        if (pauli_operator->get_x_bits() == target.get_x_bits() &&
            pauli_operator->get_z_bits() == target.get_z_bits()) {
            _operator_list[i]->change_coef(
                _operator_list[i]->get_coef() - target.get_coef());
            flag = false;
        }
    }
    if (flag) {
        auto copy = target.copy();
        copy->change_coef(-copy->get_coef());
        this->add_operator(copy);
    }
    return *this;
}

GeneralQuantumOperator GeneralQuantumOperator::operator*(
    const GeneralQuantumOperator& target) const {
    auto res = this->copy();
    *res *= target;
    return *res;
}

GeneralQuantumOperator GeneralQuantumOperator::operator*(
    const PauliOperator& target) const {
    auto res = this->copy();
    *res *= target;
    return *res;
}
GeneralQuantumOperator GeneralQuantumOperator::operator*(
    CPPCTYPE target) const {
    auto res = this->copy();
    *res *= target;
    return *res;
}

GeneralQuantumOperator& GeneralQuantumOperator::operator*=(
    const GeneralQuantumOperator& target) {
    auto copy = this->copy();
    _operator_list.clear();
#pragma omp parallel for
    for (UINT i = 0; i < copy->get_terms().size(); i++) {
        auto pauli_operator = copy->get_terms()[i];
        for (UINT j = 0; j < target.get_terms().size(); j++) {
            auto target_operator = target.get_terms()[j];
            PauliOperator* product = new PauliOperator;
            *product = (*pauli_operator) * (*target_operator);
            *this += *product;
        }
    }
    // これ必要？スコープ抜けたら消える気もする
    delete copy;
    return *this;
}

GeneralQuantumOperator& GeneralQuantumOperator::operator*=(
    const PauliOperator& target) {
    auto copy = this->copy();
    _operator_list.clear();
#pragma omp parallel for
    for (UINT i = 0; i < copy->get_terms().size(); i++) {
        auto pauli_operator = copy->get_terms()[i];
        PauliOperator* product = new PauliOperator;
        *product = (*pauli_operator) * (target);
        *this += *product;
    }
    // これ必要？スコープ抜けたら消える気もする
    delete copy;
    return *this;
}

GeneralQuantumOperator& GeneralQuantumOperator::operator*=(
    CPPCTYPE target) {
#pragma omp parallel for
    for (UINT i = 0; i < _operator_list.size(); i++) {
        *_operator_list[i] *= target;
    }
    return *this;
}

void GeneralQuantumOperator::apply_to_state(
    QuantumStateBase* state_to_be_multiplied,
    QuantumStateBase* dst_state) const {
    dst_state->multiply_coef(0.0);
    auto work_state = QuantumState(state_to_be_multiplied->qubit_count);
    const auto term_count = this->get_term_count();
    for (UINT i = 0; i < term_count; i++) {
        work_state.load(state_to_be_multiplied);
        auto term = this->get_term(i);
        auto pauli_operator =
            gate::Pauli(term->get_index_list(), term->get_pauli_id_list());
        pauli_operator->update_quantum_state(&work_state);
        work_state.multiply_coef(term->get_coef());
        dst_state->add_state(&work_state);
    }
}

CPPCTYPE GeneralQuantumOperator::calculate_default_mu() const {
    CPPCTYPE mu = 0.0;
    const auto term_count = this->get_term_count();
    for (UINT i = 0; i < term_count; i++) {
        const auto term = this->get_term(i);
        mu += term->get_coef();
    }
    return mu;
}

namespace quantum_operator {
GeneralQuantumOperator* create_general_quantum_operator_from_openfermion_file(
    std::string file_path) {
    UINT qubit_count = 0;
    std::vector<CPPCTYPE> coefs;
    std::vector<std::string> ops;

    // loading lines and check qubit_count
    double coef_real, coef_imag;
    std::string str_buf;
    std::vector<std::string> index_list;

    std::ifstream ifs;
    std::string line;
    ifs.open(file_path);

    while (getline(ifs, line)) {
        std::tuple<double, double, std::string> parsed_items =
            parse_openfermion_line(line);
        coef_real = std::get<0>(parsed_items);
        coef_imag = std::get<1>(parsed_items);
        str_buf = std::get<2>(parsed_items);

        CPPCTYPE coef(coef_real, coef_imag);
        coefs.push_back(coef);
        ops.push_back(str_buf);
        index_list = split(str_buf, "IXYZ ");

        for (UINT i = 0; i < index_list.size(); ++i) {
            UINT n = std::stoi(index_list[i]) + 1;
            if (qubit_count < n) qubit_count = n;
        }
    }
    if (!ifs.eof()) {
        std::cerr << "ERROR: Invalid format" << std::endl;
        return (GeneralQuantumOperator*)NULL;
    }
    ifs.close();

    GeneralQuantumOperator* general_quantum_operator =
        new GeneralQuantumOperator(qubit_count);

    for (UINT i = 0; i < ops.size(); ++i) {
        general_quantum_operator->add_operator(
            new PauliOperator(ops[i].c_str(), coefs[i]));
    }

    return general_quantum_operator;
}

GeneralQuantumOperator* create_general_quantum_operator_from_openfermion_text(
    std::string text) {
    UINT qubit_count = 0;
    std::vector<CPPCTYPE> coefs;
    std::vector<std::string> ops;

    double coef_real, coef_imag;
    std::string str_buf;
    std::vector<std::string> index_list;

    std::vector<std::string> lines;
    lines = split(text, "\n");
    for (std::string line : lines) {
        std::tuple<double, double, std::string> parsed_items =
            parse_openfermion_line(line);
        coef_real = std::get<0>(parsed_items);
        coef_imag = std::get<1>(parsed_items);
        str_buf = std::get<2>(parsed_items);

        CPPCTYPE coef(coef_real, coef_imag);
        coefs.push_back(coef);
        ops.push_back(str_buf);
        index_list = split(str_buf, "IXYZ ");

        for (UINT i = 0; i < index_list.size(); ++i) {
            UINT n = std::stoi(index_list[i]) + 1;
            if (qubit_count < n) qubit_count = n;
        }
    }
    GeneralQuantumOperator* general_quantum_operator =
        new GeneralQuantumOperator(qubit_count);

    for (UINT i = 0; i < ops.size(); ++i) {
        general_quantum_operator->add_operator(
            new PauliOperator(ops[i].c_str(), coefs[i]));
    }

    return general_quantum_operator;
}

std::pair<GeneralQuantumOperator*, GeneralQuantumOperator*>
create_split_general_quantum_operator(std::string file_path) {
    UINT qubit_count = 0;
    std::vector<CPPCTYPE> coefs;
    std::vector<std::string> ops;

    std::ifstream ifs;
    ifs.open(file_path);

    if (!ifs) {
        std::cerr << "ERROR: Cannot open file" << std::endl;
        return std::make_pair(
            (GeneralQuantumOperator*)NULL, (GeneralQuantumOperator*)NULL);
    }

    // loading lines and check qubit_count
    double coef_real, coef_imag;
    std::string str_buf;
    std::vector<std::string> index_list;

    std::string line;
    while (getline(ifs, line)) {
        std::tuple<double, double, std::string> parsed_items =
            parse_openfermion_line(line);
        coef_real = std::get<0>(parsed_items);
        coef_imag = std::get<1>(parsed_items);
        str_buf = std::get<2>(parsed_items);
        if (str_buf == (std::string)NULL) {
            continue;
        }
        CPPCTYPE coef(coef_real, coef_imag);
        coefs.push_back(coef);
        ops.push_back(str_buf);
        index_list = split(str_buf, "IXYZ ");

        for (UINT i = 0; i < index_list.size(); ++i) {
            UINT n = std::stoi(index_list[i]) + 1;
            if (qubit_count < n) qubit_count = n;
        }
    }
    if (!ifs.eof()) {
        std::cerr << "ERROR: Invalid format" << std::endl;
        return std::make_pair(
            (GeneralQuantumOperator*)NULL, (GeneralQuantumOperator*)NULL);
    }
    ifs.close();

    GeneralQuantumOperator* general_quantum_operator_diag =
        new GeneralQuantumOperator(qubit_count);
    GeneralQuantumOperator* general_quantum_operator_non_diag =
        new GeneralQuantumOperator(qubit_count);

    for (UINT i = 0; i < ops.size(); ++i) {
        if (ops[i].find("X") != std::string::npos ||
            ops[i].find("Y") != std::string::npos) {
            general_quantum_operator_non_diag->add_operator(
                new PauliOperator(ops[i].c_str(), coefs[i]));
        } else {
            general_quantum_operator_diag->add_operator(
                new PauliOperator(ops[i].c_str(), coefs[i]));
        }
    }

    return std::make_pair(
        general_quantum_operator_diag, general_quantum_operator_non_diag);
}
}  // namespace quantum_operator

bool check_Pauli_operator(const GeneralQuantumOperator* quantum_operator,
    const PauliOperator* pauli_operator) {
    auto vec = pauli_operator->get_index_list();
    UINT val = 0;
    if (vec.size() > 0) {
        val = std::max(val, *std::max_element(vec.begin(), vec.end()));
    }
    return val < (quantum_operator->get_qubit_count());
}
