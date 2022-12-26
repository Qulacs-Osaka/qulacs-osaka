#pragma once

#include <cppsim/type.hpp>
#include <string>
#include <vector>

#include "parametric_gate.hpp"

namespace gate {
DllExport QuantumGateBase* create_parametric_quantum_gate_from_string(
    std::string gate_string);  // old-style
DllExport QuantumGate_SingleParameter* ParametricRX(
    UINT qubit_index, double initial_angle = 0.);  // old-style
DllExport QuantumGate_SingleParameter* ParametricRX_existing_parameter(
    UINT qubit_index, const ParameterId& parameter_id,
    double parameter_coef = 1.);  // new-style
DllExport QuantumGate_SingleParameter* ParametricRY(
    UINT qubit_index, double initial_angle = 0.);  // old-style
DllExport QuantumGate_SingleParameter* ParametricRY_existing_parameter(
    UINT qubit_index, const ParameterId& parameter_id,
    double parameter_coef = 1.);  // new-style
DllExport QuantumGate_SingleParameter* ParametricRZ(
    UINT qubit_index, double initial_angle = 0.);  // old-style
DllExport QuantumGate_SingleParameter* ParametricRZ_existing_parameter(
    UINT qubit_index, const ParameterId& parameter_id,
    double parameter_coef = 1.);  // new-style
DllExport QuantumGate_SingleParameter* ParametricPauliRotation(
    std::vector<UINT> target, std::vector<UINT> pauli_id,
    double initial_angle = 0.);  // old-style
DllExport QuantumGate_SingleParameter*
ParametricPauliRotation_existing_parameter(std::vector<UINT> target,
    std::vector<UINT> pauli_id, const ParameterId& parameter_id,
    double parameter_coef = 1.);  // new-style
}  // namespace gate
