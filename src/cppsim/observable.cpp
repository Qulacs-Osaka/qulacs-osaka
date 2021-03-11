#include "observable.hpp"

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include <fstream>
#include <iostream>

#include "state.hpp"
#include "utility.hpp"

void HermitianQuantumOperator::add_operator(const PauliOperator* mpt) {
    if (std::abs(mpt->get_coef().imag()) > 0) {
        std::cerr << "Error: HermitianQuantumOperator::add_operator(const "
                     "PauliOperator* mpt): PauliOperator must be Hermitian."
                  << std::endl;
        return;
    }
    GeneralQuantumOperator::add_operator(mpt);
}

void HermitianQuantumOperator::add_operator(
    CPPCTYPE coef, std::string pauli_string) {
    if (std::abs(coef.imag()) > 0) {
        std::cerr << "Error: HermitianQuantumOperator::add_operator(const "
                     "PauliOperator* mpt): PauliOperator must be Hermitian."
                  << std::endl;
        return;
    }
    GeneralQuantumOperator::add_operator(coef, pauli_string);
}

CPPCTYPE HermitianQuantumOperator::get_expectation_value(
    const QuantumStateBase* state) const {
    return GeneralQuantumOperator::get_expectation_value(state).real();
}

CPPCTYPE
HermitianQuantumOperator::solve_ground_state_eigenvalue_by_arnoldi_method(
    QuantumStateBase* state, const UINT iter_count, const CPPCTYPE mu) const {
    return GeneralQuantumOperator::
        solve_ground_state_eigenvalue_by_arnoldi_method(state, iter_count, mu)
            .real();
}

CPPCTYPE
HermitianQuantumOperator::solve_ground_state_eigenvalue_by_power_method(
    QuantumStateBase* state, const UINT iter_count, const CPPCTYPE mu) const {
    return GeneralQuantumOperator::
        solve_ground_state_eigenvalue_by_power_method(state, iter_count, mu)
            .real();
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
        std::cerr << "ERROR: Cannot open file" << std::endl;
        return NULL;
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
        std::cerr << "ERROR: Invalid format" << std::endl;
        return NULL;
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
        std::cerr << "ERROR: Cannot open file" << std::endl;
        return std::make_pair(
            (HermitianQuantumOperator*)NULL, (HermitianQuantumOperator*)NULL);
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
        std::cerr << "ERROR: Invalid format" << std::endl;
        return std::make_pair(
            (HermitianQuantumOperator*)NULL, (HermitianQuantumOperator*)NULL);
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
