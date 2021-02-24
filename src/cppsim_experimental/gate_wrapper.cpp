
#include "gate_wrapper.hpp"

#include "gate_basic.hpp"

namespace gate {
DllExport QuantumGateWrapped* DepolarizingNoise(UINT index, double prob) {
    auto ptr = QuantumGateWrapped::ProbabilisticGate(

        {1 - prob, prob / 3, prob / 3, prob / 3},
        {gate::Identity(index), gate::X(index), gate::Y(index), gate::Z(index)},
        true);
    return ptr;
}
DllExport QuantumGateWrapped* IndependentXZNoise(UINT index, double prob) {
    auto ptr = QuantumGateWrapped::ProbabilisticGate(
        {(1 - prob) * (1 - prob), prob * (1 - prob), (1 - prob) * prob,
            prob * prob},
        {gate::Identity(index), gate::X(index), gate::Z(index), gate::Y(index)},
        true);
    return ptr;
}

DllExport QuantumGateWrapped* TwoQubitDepolarizingNoise(
    UINT index1, UINT index2, double prob) {
    std::vector<QuantumGateBase*> gates;
    std::vector<double> probs;
    probs.push_back(1 - prob);
    gates.push_back(gate::Identity(index1));
    for (UINT i = 1; i < 16; ++i) {
        auto gate = QuantumGateBasic::PauliMatrixGate(
            {index1, index2}, {i % 4, i / 4}, 0.);
        gates.push_back(gate);
        probs.push_back(prob / 15);
    }
    auto ptr = QuantumGateWrapped::ProbabilisticGate(probs, gates, true);
    return ptr;
}
DllExport QuantumGateWrapped* BitFlipNoise(UINT index, double prob) {
    auto ptr = QuantumGateWrapped::ProbabilisticGate(
        {1 - prob, prob}, {gate::Identity(index), gate::X(index)}, true);
    return ptr;
}
DllExport QuantumGateWrapped* DephasingNoise(UINT index, double prob) {
    auto ptr = QuantumGateWrapped::ProbabilisticGate(
        {1 - prob, prob}, {gate::Identity(index), gate::Z(index)}, true);
    return ptr;
}
DllExport QuantumGateWrapped* AmplitudeDampingNoise(
    UINT target_index, double prob) {
    ComplexMatrix damping_matrix_0(2, 2), damping_matrix_1(2, 2);
    damping_matrix_0 << 1, 0, 0, sqrt(1 - prob);
    damping_matrix_1 << 0, sqrt(prob), 0, 0;
    auto gate0 =
        QuantumGateBasic::DenseMatrixGate({target_index}, damping_matrix_0);
    auto gate1 =
        QuantumGateBasic::DenseMatrixGate({target_index}, damping_matrix_1);
    auto new_gate = QuantumGateWrapped::CPTP({gate0, gate1});
    delete gate0;
    delete gate1;
    return new_gate;
}
DllExport QuantumGateWrapped* CPTP(std::vector<QuantumGateBase*> gate_list) {
    return QuantumGateWrapped::CPTP(gate_list);
}
DllExport QuantumGateWrapped* Probabilistic(
    std::vector<double> distribution, std::vector<QuantumGateBase*> gate_list) {
    auto ptr =
        QuantumGateWrapped::ProbabilisticGate(distribution, gate_list, false);
    return ptr;
}
}  // namespace gate
