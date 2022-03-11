#include "observable.hpp"

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#undef NDEBUG
#include <Eigen/Dense>
#include <array>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>

#include "state.hpp"
#include "utility.hpp"

void HermitianQuantumOperator::add_operator(const PauliOperator* mpt) {
    if (std::abs(mpt->get_coef().imag()) > 0) {
        std::stringstream error_message_stream;
        error_message_stream
            << "Error: HermitianQuantumOperator::add_operator(const "
               "PauliOperator* mpt): PauliOperator must be Hermitian.";
        throw std::invalid_argument(error_message_stream.str());
    }
    GeneralQuantumOperator::add_operator(mpt);
}

void HermitianQuantumOperator::add_operator(
    CPPCTYPE coef, std::string pauli_string) {
    if (std::abs(coef.imag()) > 0) {
        std::stringstream error_message_stream;
        error_message_stream
            << "Error: HermitianQuantumOperator::add_operator(const "
               "PauliOperator* mpt): PauliOperator must be Hermitian.";
        throw std::invalid_argument(error_message_stream.str());
    }
    GeneralQuantumOperator::add_operator(coef, pauli_string);
}

CPPCTYPE HermitianQuantumOperator::get_expectation_value(
    const QuantumStateBase* state) const {
    return GeneralQuantumOperator::get_expectation_value(state).real();
}

CPPCTYPE
HermitianQuantumOperator::solve_ground_state_eigenvalue_by_lanczos_method(
    QuantumStateBase* init_state, const UINT iter_count,
    const CPPCTYPE mu) const {
    if (this->get_term_count() == 0) {
        std::stringstream error_message_stream;
        error_message_stream
            << "Error: "
               "HermitianQuantumOperator::solve_ground_state_eigenvalue_"
               "by_lanczos_method("
               "QuantumStateBase * state, const UINT iter_count, const "
               "CPPCTYPE mu): At least one PauliOperator is required.";
        throw std::invalid_argument(error_message_stream.str());
    }

    // Implemented based on
    // https://files.transtutors.com/cdn/uploadassignments/472339_1_-numerical-linear-aljebra.pdf
    CPPCTYPE mu_;
    if (mu == 0.0) {
        // mu is not changed from default value.
        mu_ = this->calculate_default_mu();
    } else {
        mu_ = mu;
    }

    const auto qubit_count = this->get_qubit_count();
    QuantumState tmp_state(qubit_count);
    QuantumState mu_timed_state(qubit_count);
    // work_states: [q_{i-1}, q_i, q_{i+1}]
    // q_0, q_1, q_2,... span Krylov subspace.
    init_state->normalize(init_state->get_squared_norm());
    std::array<QuantumState, 3> work_states = {QuantumState(qubit_count),
        QuantumState(qubit_count), QuantumState(qubit_count)};
    work_states.at(1).load(init_state);

    Eigen::VectorXd alpha_v(iter_count);
    Eigen::VectorXd beta_v(iter_count - 1);
    for (UINT i = 0; i < iter_count; i++) {
        // v = (A - μI) * q_i
        mu_timed_state.load(&work_states.at(1));
        mu_timed_state.multiply_coef(-mu_);
        this->apply_to_state(&tmp_state, work_states.at(1), &work_states.at(2));
        work_states.at(2).add_state(&mu_timed_state);

        // α_i = q_i^T * v
        alpha_v(i) =
            state::inner_product(&work_states.at(1), &work_states.at(2)).real();
        // In the last iteration, no need to calculate β.
        if (i == iter_count - 1) {
            break;
        }
        // v -= α_i * q_i
        tmp_state.load(&work_states.at(1));
        tmp_state.multiply_coef(-alpha_v(i));
        work_states.at(2).add_state(&tmp_state);
        if (i != 0) {
            // v -= β_{i-1} * q_{i-1}
            tmp_state.load(&work_states.at(0));
            tmp_state.multiply_coef(-beta_v(i - 1));
            work_states.at(2).add_state(&tmp_state);
        }

        // β_i = ||v||
        beta_v(i) = std::sqrt(work_states.at(2).get_squared_norm());
        // q_{i+1} = v / β_i
        work_states.at(2).multiply_coef(1 / beta_v(i));
        work_states.at(0).load(&work_states.at(1));
        work_states.at(1).load(&work_states.at(2));
    }

    // Compute eigenvalue of a symmetric matrix T whose diagonal elements are
    // `alpha_v` and subdiagonal elements are `beta_v`.
    Eigen::SelfAdjointEigenSolver<ComplexMatrix> solver;
    solver.computeFromTridiagonal(alpha_v, beta_v);
    const auto eigenvalues = solver.eigenvalues();
    // Find ground state eigenvalue.
    UINT minimum_eigenvalue_index = 0;
    auto minimum_eigenvalue = eigenvalues(0);
    for (UINT i = 0; i < eigenvalues.size(); i++) {
        if (eigenvalues(i) < minimum_eigenvalue) {
            minimum_eigenvalue_index = i;
            minimum_eigenvalue = eigenvalues(i);
        }
    }

    auto eigenvectors = solver.eigenvectors();
    auto eigenvector_in_krylov = eigenvectors.col(minimum_eigenvalue_index);
    // Store ground state eigenvector to `init_state`.
    // If λ is an eigenvalue of T and q is the eigenvector, Tq = λq.
    // And let V be a matrix whose column vectors span Krylov subspace.
    // Then, T = V^* AV where A is this observable.
    // Tq = λq, VV^* AVq = Vλq, A(Vq) = λ(Vq).
    // So, an eigenvector of A for λ is Vq.
    // q_0 = init_state
    work_states.at(1).load(init_state);
    init_state->set_zero_norm_state();
    assert(eigenvector_in_krylov.size() == iter_count);
    for (UINT i = 0; i < iter_count; i++) {
        // q += v_i * q_i, where q is eigenvector to compute
        tmp_state.load(&work_states.at(1));
        tmp_state.multiply_coef(eigenvector_in_krylov(i));
        init_state->add_state(&tmp_state);

        // v = (A - μI) * q_i
        mu_timed_state.load(&work_states.at(1));
        mu_timed_state.multiply_coef(-mu_);
        this->apply_to_state(&tmp_state, work_states.at(1), &work_states.at(2));
        work_states.at(2).add_state(&mu_timed_state);
        if (i == iter_count - 1) {
            break;
        }

        // v -= α_i * q_i
        tmp_state.load(&work_states.at(1));
        tmp_state.multiply_coef(-alpha_v(i));
        work_states.at(2).add_state(&tmp_state);
        if (i != 0) {
            // v -= β_{i-1} * q_{i-1}
            tmp_state.load(&work_states.at(0));
            tmp_state.multiply_coef(-beta_v(i - 1));
            work_states.at(2).add_state(&tmp_state);
        }

        // q_{i+1} = v / β_i
        work_states.at(2).multiply_coef(1 / beta_v[i]);
        work_states.at(0).load(&work_states.at(1));
        work_states.at(1).load(&work_states.at(2));
    }

    return minimum_eigenvalue + mu_;
}

std::string HermitianQuantumOperator::to_string() const {
    std::stringstream os;
    auto term_count = this->get_term_count();
    for (UINT index = 0; index < term_count; index++) {
        os << this->get_term(index)->get_coef().real() << " ";
        os << this->get_term(index)->get_pauli_string();
        if (index != term_count - 1) {
            os << " + ";
        }
    }
    return os.str();
}

namespace observable {
HermitianQuantumOperator* create_observable_from_openfermion_file(
    std::string file_path) {
    UINT qubit_count = 0;
    std::vector<CPPCTYPE> coefs;
    std::vector<std::string> ops;

    std::ifstream ifs;
    ifs.open(file_path);

    if (!ifs) {
        std::stringstream error_message_stream;
        error_message_stream << "ERROR: Cannot open file";
        throw std::runtime_error(error_message_stream.str());
    }

    // loading lines and check qubit_count
    std::string str_buf;
    std::vector<std::string> index_list;

    std::string line;
    while (getline(ifs, line)) {
        std::tuple<double, double, std::string> parsed_items =
            parse_openfermion_line(line);
        const auto coef_real = std::get<0>(parsed_items);
        const auto coef_imag = std::get<1>(parsed_items);
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
        std::stringstream error_message_stream;
        error_message_stream << "ERROR: Invalid format";
        throw std::runtime_error(error_message_stream.str());
    }
    ifs.close();

    HermitianQuantumOperator* observable =
        new HermitianQuantumOperator(qubit_count);

    for (UINT i = 0; i < ops.size(); ++i) {
        observable->add_operator(new PauliOperator(ops[i].c_str(), coefs[i]));
    }

    return observable;
}

HermitianQuantumOperator* create_observable_from_openfermion_text(
    const std::string& text) {
    UINT qubit_count = 0;
    std::vector<CPPCTYPE> coefs;
    std::vector<std::string> ops;

    std::vector<std::string> lines;
    std::string str_buf;
    std::vector<std::string> index_list;

    lines = split(text, "\n");
    for (std::string line : lines) {
        std::tuple<double, double, std::string> parsed_items =
            parse_openfermion_line(line);
        const auto coef_real = std::get<0>(parsed_items);
        const auto coef_imag = std::get<1>(parsed_items);
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
    HermitianQuantumOperator* hermitian_quantum_operator =
        new HermitianQuantumOperator(qubit_count);

    for (UINT i = 0; i < ops.size(); ++i) {
        hermitian_quantum_operator->add_operator(
            new PauliOperator(ops[i].c_str(), coefs[i]));
    }

    return hermitian_quantum_operator;
}

std::pair<HermitianQuantumOperator*, HermitianQuantumOperator*>
create_split_observable(std::string file_path) {
    UINT qubit_count = 0;
    std::vector<CPPCTYPE> coefs;
    std::vector<std::string> ops;

    std::ifstream ifs;
    ifs.open(file_path);

    if (!ifs) {
        std::stringstream error_message_stream;
        error_message_stream << "ERROR: Cannot open file";
        throw std::runtime_error(error_message_stream.str());
    }

    // loading lines and check qubit_count
    std::string str_buf;
    std::vector<std::string> index_list;

    std::string line;
    while (getline(ifs, line)) {
        std::tuple<double, double, std::string> parsed_items =
            parse_openfermion_line(line);
        const auto coef_real = std::get<0>(parsed_items);
        const auto coef_imag = std::get<1>(parsed_items);
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
        std::stringstream error_message_stream;
        error_message_stream << "ERROR: Invalid format";
        throw std::runtime_error(error_message_stream.str());
    }
    ifs.close();

    HermitianQuantumOperator* observable_diag =
        new HermitianQuantumOperator(qubit_count);
    HermitianQuantumOperator* observable_non_diag =
        new HermitianQuantumOperator(qubit_count);

    for (UINT i = 0; i < ops.size(); ++i) {
        if (ops[i].find("X") != std::string::npos ||
            ops[i].find("Y") != std::string::npos) {
            observable_non_diag->add_operator(
                new PauliOperator(ops[i].c_str(), coefs[i]));
        } else {
            observable_diag->add_operator(
                new PauliOperator(ops[i].c_str(), coefs[i]));
        }
    }

    return std::make_pair(observable_diag, observable_non_diag);
}
}  // namespace observable
