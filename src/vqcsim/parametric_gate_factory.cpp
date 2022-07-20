#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#define _USE_MATH_DEFINES
#include "parametric_gate_factory.hpp"

#include <cmath>
#include <cppsim/exception.hpp>
#include <cppsim/gate_factory.hpp>
#include <cppsim/utility.hpp>
#include <cstdlib>
#include <cstring>

#include "parametric_gate.hpp"

namespace gate {
QuantumGate_SingleParameter* ParametricRX(UINT target_qubit_index,
    const ParameterKey& parameter_id, double parameter_coef) {
    return new ClsParametricRXGate(
        target_qubit_index, parameter_id, parameter_coef);
}
QuantumGate_SingleParameter* ParametricRY(UINT target_qubit_index,
    const ParameterKey& parameter_id, double parameter_coef) {
    return new ClsParametricRYGate(
        target_qubit_index, parameter_id, parameter_coef);
}
QuantumGate_SingleParameter* ParametricRZ(UINT target_qubit_index,
    const ParameterKey& parameter_id, double parameter_coef) {
    return new ClsParametricRZGate(
        target_qubit_index, parameter_id, parameter_coef);
}
QuantumGate_SingleParameter* ParametricPauliRotation(std::vector<UINT> target,
    std::vector<UINT> pauli_id, const ParameterKey& parameter_id,
    double parameter_coef) {
    if (!check_is_unique_index_list(target)) {
        throw DuplicatedQubitIndexException(
            "Error: gate::ParametricPauliRotation(std::vector<UINT>, "
            "std::vector<UINT>, double): target qubit list contains "
            "duplicated values."
            "\nInfo: NULL used to be returned, "
            "but it changed to throw exception.");
    }
    auto pauli = new PauliOperator(target, pauli_id, 0.);
    return new ClsParametricPauliRotationGate(
        pauli, parameter_id, parameter_coef);
}

QuantumGateBase* create_parametric_quantum_gate_from_string(
    std::string gate_string, const ParameterKey& parameter_id) {
    auto non_parametric_gate =
        gate::create_quantum_gate_from_string(gate_string);
    if (non_parametric_gate != NULL) return non_parametric_gate;

    const char* gateString = gate_string.c_str();
    char* sbuf;
    // ITYPE elementCount;
    std::vector<CPPCTYPE> element;
    const char delim[] = " ";
    std::vector<UINT> targets;
    QuantumGateBase* gate = NULL;
    char* buf = (char*)calloc(strlen(gateString) + 1, sizeof(char));
    strcpy(buf, gateString);
    sbuf = strtok(buf, delim);

    if (strcasecmp(sbuf, "PRX") == 0) {
        unsigned int target = atoi(strtok(NULL, delim));
        gate = gate::ParametricRX(target, parameter_id);
    } else if (strcasecmp(sbuf, "PRY") == 0) {
        unsigned int target = atoi(strtok(NULL, delim));
        gate = gate::ParametricRY(target, parameter_id);
    } else if (strcasecmp(sbuf, "PRZ") == 0) {
        unsigned int target = atoi(strtok(NULL, delim));
        gate = gate::ParametricRZ(target, parameter_id);
    } else if (strcasecmp(sbuf, "PPR") == 0) {
        char* pauliStr = strtok(NULL, delim);
        unsigned int targetCount = (UINT)strlen(pauliStr);

        std::vector<UINT> pauli(targetCount, 0);
        for (unsigned int i = 0; i < targetCount; i++) {
            if (pauliStr[i] == 'x' || pauliStr[i] == 'X')
                pauli[i] = 1;
            else if (pauliStr[i] == 'y' || pauliStr[i] == 'Y')
                pauli[i] = 2;
            else if (pauliStr[i] == 'z' || pauliStr[i] == 'Z')
                pauli[i] = 3;
        }

        targets = std::vector<UINT>(targetCount, 0);
        for (unsigned int i = 0; i < targetCount; i++) {
            targets[i] = atoi(strtok(NULL, delim));
        }
        gate = gate::ParametricPauliRotation(targets, pauli, parameter_id);
    }
    free(buf);
    return gate;
}

}  // namespace gate
