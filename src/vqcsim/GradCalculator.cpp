#define _USE_MATH_DEFINES
#include "GradCalculator.hpp"

#include <math.h>

#include "causalcone_simulator.hpp"

std::vector<std::complex<double>> GradCalculator::calculate_grad(
    ParametricQuantumCircuit& x, Observable& obs,
    const ParameterSet& parameter_set) {
    ParametricQuantumCircuit* x_copy = x.copy();
    std::vector<ParameterKey> parameter_id_list = x.get_parameter_id_list();

    std::vector<std::complex<double>> grad;
    for (ParameterKey& i : parameter_id_list) {
        std::complex<double> y, z;
        {
            for (ParameterKey& q : parameter_id_list) {
                double diff = 0;
                if (i == q) {
                    diff = M_PI_2;
                }
                auto parameter_it = parameter_set.find(q);
                if (parameter_it == parameter_set.end()) {
                    throw ParameterIdNotFoundException(
                        "Error: "
                        "GradCalculator::calculate_grad("
                        "ParametricQuantumCircuit&, Observable&, ParameterSet, "
                        "const ParameterSet&): given parameter_set does not "
                        "contain a parameter: " +
                        q);
                }
                x_copy->set_parameter(q, parameter_it->second + diff);
            }
            CausalConeSimulator hoge(*x_copy, obs);
            y = hoge.get_expectation_value();
        }
        {
            for (ParameterKey& q : parameter_id_list) {
                double diff = 0;
                if (i == q) {
                    diff = M_PI_2;
                }
                auto parameter_it = parameter_set.find(q);
                if (parameter_it == parameter_set.end()) {
                    throw ParameterIdNotFoundException(
                        "Error: "
                        "GradCalculator::calculate_grad("
                        "ParametricQuantumCircuit&, Observable&, ParameterSet, "
                        "const ParameterSet&): given parameter_set does not "
                        "contain a parameter: " +
                        q);
                }
                x_copy->set_parameter(q, parameter_it->second - diff);
            }
            CausalConeSimulator hoge(*x_copy, obs);
            z = hoge.get_expectation_value();
        }
        grad.push_back((y - z) / 2.0);
    }
    delete x_copy;
    return grad;
};

std::vector<std::complex<double>> GradCalculator::calculate_grad(
    ParametricQuantumCircuit& x, Observable& obs, std::vector<double> theta) {
    ParameterSet parameter_set;
    for (UINT i = 0; i < theta.size(); ++i) {
        parameter_set[x.generate_parameter_id_from_index(i)] = theta[i];
    }
    calculate_grad(x, obs, parameter_set);
};

std::vector<std::complex<double>> GradCalculator::calculate_grad(
    ParametricQuantumCircuit& x, Observable& obs) {
    return calculate_grad(x, obs, x.get_parameter_set());
};
