#include <gtest/gtest.h>

#include <cppsim/exception.hpp>
#include <cppsim/gate_factory.hpp>
#include <cppsim/state_dm.hpp>
#include <vqcsim/GradCalculator.hpp>
#include <vqcsim/causalcone_simulator.hpp>
#include <vqcsim/parametric_circuit_builder.hpp>
#include <vqcsim/problem.hpp>
#include <vqcsim/solver.hpp>

#include "../util/util.hpp"

class ClsParametricNullUpdateGate
    : public QuantumGate_SingleParameterOneQubitRotation {
public:
    ClsParametricNullUpdateGate(UINT target_qubit_index,
        SingleParameter* parameter, decltype(_angle_func) angle_func)
        : QuantumGate_SingleParameterOneQubitRotation(parameter, angle_func) {
        this->_name = "ParametricNullUpdate";
        this->_target_qubit_list.push_back(TargetQubitInfo(target_qubit_index));
    }
    virtual void set_matrix(ComplexMatrix& matrix) const override {}
    virtual QuantumGate_SingleParameter* copy() const override {
        return new ClsParametricNullUpdateGate(*this);
    };
};

TEST(ParametricGate, NullUpdateFunc) {
    SingleParameter parameter(0., 0);
    ClsParametricNullUpdateGate gate(0, &parameter, identity_map);
    QuantumState state(1);
    ASSERT_THROW(
        gate.update_quantum_state(&state), UndefinedUpdateFuncException);
}

TEST(ParametricCircuit, GateApply) {
    const UINT n = 3;
    const UINT depth = 10;
    ParametricQuantumCircuit* circuit = new ParametricQuantumCircuit(n);
    Random random;
    for (UINT d = 0; d < depth; ++d) {
        for (UINT i = 0; i < n; ++i) {
            circuit->add_parametric_RX_gate(i, random.uniform());
            circuit->add_parametric_RY_gate(i, random.uniform());
            circuit->add_parametric_RZ_gate(i, random.uniform());
        }
        for (UINT i = d % 2; i + 1 < n; i += 2) {
            circuit->add_parametric_multi_Pauli_rotation_gate(
                {i, i + 1}, {3, 3}, random.uniform());
        }
    }

    UINT param_count = circuit->get_parameter_count();
    for (UINT p = 0; p < param_count; ++p) {
        double current_angle = circuit->get_parameter(p);
        circuit->set_parameter(p, current_angle + random.uniform());
    }

    QuantumState state(n);
    circuit->update_quantum_state(&state);
    // std::cout << state << std::endl;
    // std::cout << circuit << std::endl;
    delete circuit;
}

TEST(ParametricCircuit, GateApplyDM) {
    const UINT n = 3;
    const UINT depth = 10;
    ParametricQuantumCircuit* circuit = new ParametricQuantumCircuit(n);
    Random random;
    for (UINT d = 0; d < depth; ++d) {
        for (UINT i = 0; i < n; ++i) {
            circuit->add_parametric_RX_gate(i, random.uniform());
            circuit->add_parametric_RY_gate(i, random.uniform());
            circuit->add_parametric_RZ_gate(i, random.uniform());
        }
        for (UINT i = d % 2; i + 1 < n; i += 2) {
            circuit->add_parametric_multi_Pauli_rotation_gate(
                {i, i + 1}, {3, 3}, random.uniform());
        }
    }

    UINT param_count = circuit->get_parameter_count();
    for (UINT p = 0; p < param_count; ++p) {
        double current_angle = circuit->get_parameter(p);
        circuit->set_parameter(p, current_angle + random.uniform());
    }

    DensityMatrix state(n);
    circuit->update_quantum_state(&state);
    // std::cout << state << std::endl;
    // std::cout << circuit << std::endl;
    delete circuit;
}

TEST(ParametricCircuit, ParametricGatePosition) {
    auto circuit = ParametricQuantumCircuit(3);
    circuit.add_parametric_RX_gate(0, 0.);
    circuit.add_H_gate(0);
    circuit.add_parametric_RZ_gate(0, 2.);
    circuit.add_gate_copy(gate::CNOT(0, 1));
    circuit.add_parametric_RY_gate(1, 1U);
    circuit.add_gate_copy(gate::X(0), 2);
    circuit.add_parametric_RZ_gate(1, 0U);
    circuit.remove_gate(4);
    circuit.remove_gate(5);
    circuit.add_parametric_multi_Pauli_rotation_gate({1}, {0}, 3.);
    circuit.set_parameter(0, 1.);

    // RX(0), H, X, RZ(1), RZ(0), multi(2)
    ASSERT_EQ(circuit.get_parameter_count(), 3);
    ASSERT_NEAR(circuit.get_parameter(0), 1., 1e-12);
    ASSERT_NEAR(circuit.get_parameter(1), 2., 1e-12);
    ASSERT_NEAR(circuit.get_parameter(2), 3., 1e-12);
    ASSERT_EQ(circuit.get_parametric_gate_position(0), 0);
    ASSERT_EQ(circuit.get_parametric_gate_position(1), 3);
    ASSERT_EQ(circuit.get_parametric_gate_position(2), 4);
    ASSERT_EQ(circuit.get_parametric_gate_position(3), 5);
}

class MyRandomCircuit : public ParametricCircuitBuilder {
    ParametricQuantumCircuit* create_circuit(
        UINT output_dim, UINT param_count) const override {
        ParametricQuantumCircuit* circuit =
            new ParametricQuantumCircuit(output_dim);
        UINT depth = param_count / output_dim;
        if (param_count % output_dim > 0) depth++;
        UINT param_index = 0;
        for (UINT d = 0; d < depth; ++d) {
            for (UINT i = 0; i < output_dim; ++i) {
                if (param_index < param_count) {
                    circuit->add_parametric_RX_gate(1, 0.);
                    param_index++;
                } else {
                    circuit->add_gate(gate::RX(i, 0.0));
                }
            }
            for (UINT i = depth % 2; i + 1 < output_dim; ++i) {
                circuit->add_gate(gate::CNOT(0, 1));
            }
        }
        return circuit;
    }
};

TEST(EnergyMinimization, SingleQubitClassical) {
    const UINT n = 1;

    // define quantum circuit as prediction model
    std::function<ParametricQuantumCircuit*(UINT, UINT)> func =
        [](unsigned int qubit_count,
            unsigned int param_count) -> ParametricQuantumCircuit* {
        ParametricQuantumCircuit* circuit =
            new ParametricQuantumCircuit(qubit_count);
        for (unsigned int i = 0; i < qubit_count; ++i) {
            circuit->add_parametric_RX_gate(i, 0.);
        }
        return circuit;
    };

    Observable* observable = new Observable(n);
    observable->add_operator(1.0, "Z 0");

    EnergyMinimizationProblem* emp = new EnergyMinimizationProblem(observable);

    QuantumCircuitEnergyMinimizationSolver qcems(&func, 0);
    qcems.solve(emp, 1000, "GD");
    double qc_loss = qcems.get_loss();

    DiagonalizationEnergyMinimizationSolver dems;
    dems.solve(emp);
    double diag_loss = dems.get_loss();

    EXPECT_NEAR(qc_loss, diag_loss, 1e-2);
}

TEST(EnergyMinimization, SingleQubitComplex) {
    const UINT n = 1;

    // define quantum circuit as prediction model
    std::function<ParametricQuantumCircuit*(UINT, UINT)> func =
        [](unsigned int qubit_count,
            unsigned int param_count) -> ParametricQuantumCircuit* {
        ParametricQuantumCircuit* circuit =
            new ParametricQuantumCircuit(qubit_count);
        for (unsigned int i = 0; i < qubit_count; ++i) {
            circuit->add_parametric_RX_gate(i, 0.);
            circuit->add_parametric_RY_gate(i, 0.);
            circuit->add_parametric_RX_gate(i, 0.);
        }
        return circuit;
    };

    Observable* observable = new Observable(n);
    observable->add_operator(1.0, "Z 0");
    observable->add_operator(1.0, "X 0");
    observable->add_operator(1.0, "Y 0");

    EnergyMinimizationProblem* emp = new EnergyMinimizationProblem(observable);

    QuantumCircuitEnergyMinimizationSolver qcems(&func, 0);
    qcems.solve(emp, 1000, "GD");
    double qc_loss = qcems.get_loss();

    DiagonalizationEnergyMinimizationSolver dems;
    dems.solve(emp);
    double diag_loss = dems.get_loss();

    EXPECT_NEAR(qc_loss, diag_loss, 1e-2);
}

TEST(EnergyMinimization, MultiQubit) {
    const UINT n = 2;

    // define quantum circuit as prediction model
    std::function<ParametricQuantumCircuit*(UINT, UINT)> func =
        [](unsigned int qubit_count,
            unsigned int param_count) -> ParametricQuantumCircuit* {
        ParametricQuantumCircuit* circuit =
            new ParametricQuantumCircuit(qubit_count);
        for (unsigned int i = 0; i < qubit_count; ++i) {
            circuit->add_parametric_RX_gate(i, 0.);
            circuit->add_parametric_RY_gate(i, 0.);
            circuit->add_parametric_RX_gate(i, 0.);
        }
        for (unsigned int i = 0; i + 1 < qubit_count; i += 2) {
            circuit->add_CNOT_gate(i, i + 1);
        }
        for (unsigned int i = 0; i < qubit_count; ++i) {
            circuit->add_parametric_RX_gate(i, 0.);
            circuit->add_parametric_RY_gate(i, 0.);
            circuit->add_parametric_RX_gate(i, 0.);
        }
        return circuit;
    };

    Observable* observable = new Observable(n);
    observable->add_operator(1.0, "Z 0 X 1");
    observable->add_operator(-1.0, "Z 0 Y 1");
    observable->add_operator(0.2, "Y 0 Y 1");

    EnergyMinimizationProblem* emp = new EnergyMinimizationProblem(observable);

    QuantumCircuitEnergyMinimizationSolver qcems(&func, 0);
    qcems.solve(emp, 1000, "GD");
    double qc_loss = qcems.get_loss();

    DiagonalizationEnergyMinimizationSolver dems;
    dems.solve(emp);
    double diag_loss = dems.get_loss();
    // std::cout << qc_loss << " " << diag_loss << std::endl;
    ASSERT_GT(qc_loss, diag_loss);
    EXPECT_NEAR(qc_loss, diag_loss, 1e-1);
}

TEST(GradCalculator, BasicCheck) {
    Random rnd;
    unsigned int n = 5;
    Observable observable(n);
    std::string Pauli_string = "";
    for (int i = 0; i < n; ++i) {
        double coef = rnd.uniform();
        std::string Pauli_string = "Z ";
        Pauli_string += std::to_string(i);
        observable.add_operator(coef, Pauli_string.c_str());
    }

    ParametricQuantumCircuit circuit(n);
    int cnter_parametric_gate = 0;
    for (int depth = 0; depth < 2; ++depth) {
        for (int i = 0; i < n; ++i) {
            circuit.add_parametric_RX_gate(i, 0.);
            circuit.add_parametric_RZ_gate(i, 0.);
            cnter_parametric_gate += 2;
        }

        for (int i = 0; i + 1 < n; i += 2) {
            circuit.add_CNOT_gate(i, i + 1);
        }

        for (int i = 1; i + 1 < n; i += 2) {
            circuit.add_CNOT_gate(i, i + 1);
        }
    }

    // Calculate using GradCalculator
    GradCalculator hoge;
    std::vector<double> theta;
    for (int i = 0; i < cnter_parametric_gate; ++i) {
        theta.push_back(rnd.uniform() * 5.0);
    }
    auto GradCalculator_ans_theta_specified =
        hoge.calculate_grad(circuit, observable, theta);

    for (UINT i = 0; i < cnter_parametric_gate; ++i) {
        ASSERT_EQ(circuit.get_parameter(i), 0);
        circuit.set_parameter(i, theta[i]);
    }
    auto GradCalculator_ans_without_theta =
        hoge.calculate_grad(circuit, observable);

    // Calculate using normal Greedy.
    std::vector<std::complex<double>> Greedy_ans;
    {
        for (int i = 0; i < circuit.get_parameter_count(); ++i) {
            std::complex<double> y, z;
            {
                for (int q = 0; q < circuit.get_parameter_count(); ++q) {
                    float diff = 0;
                    if (i == q) {
                        diff = 0.001;
                    }
                    circuit.set_parameter(q, theta[q] + diff);
                }
                CausalConeSimulator cone(circuit, observable);
                y = cone.get_expectation_value();
            }
            {
                for (int q = 0; q < circuit.get_parameter_count(); ++q) {
                    float diff = 0;
                    if (i == q) {
                        diff = 0.001;
                    }
                    circuit.set_parameter(q, theta[q] - diff);
                }
                CausalConeSimulator cone(circuit, observable);
                z = cone.get_expectation_value();
            }
            Greedy_ans.push_back((y - z) / 0.002);
        }
    }
    for (int i = 0; i < GradCalculator_ans_without_theta.size(); ++i) {
        ASSERT_LT(abs(GradCalculator_ans_theta_specified[i] - Greedy_ans[i]),
            1e-6);  // Difference should be lower than 1e-7
        ASSERT_LT(abs(GradCalculator_ans_without_theta[i] - Greedy_ans[i]),
            1e-6);  // Difference should be lower than 1e-7
    }
}
