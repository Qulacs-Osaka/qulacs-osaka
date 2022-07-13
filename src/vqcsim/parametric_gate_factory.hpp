#pragma once

#include <cppsim/type.hpp>
#include <string>
#include <vector>

#include "parametric_gate.hpp"

namespace gate {
DllExport QuantumGateBase* create_parametric_quantum_gate_from_string(
    std::string gate_string, const ParameterKey& parameter_id);
DllExport QuantumGate_SingleParameter* ParametricRX(
    UINT qubit_index, const ParameterKey& parameter_id);
DllExport QuantumGate_SingleParameter* ParametricRY(
    UINT qubit_index, const ParameterKey& parameter_id);
DllExport QuantumGate_SingleParameter* ParametricRZ(
    UINT qubit_index, const ParameterKey& parameter_id);
DllExport QuantumGate_SingleParameter* ParametricPauliRotation(
    std::vector<UINT> target, std::vector<UINT> pauli_id,
    const ParameterKey& parameter_id);
}  // namespace gate
