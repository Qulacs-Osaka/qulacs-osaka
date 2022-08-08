#pragma once

#include <algorithm>
#include <cassert>
#include <cppsim/exception.hpp>
#include <csim/type.hpp>
#include <map>
#include <string>
#include <vector>

using ParameterType = std::string;
using ParameterId = std::string;
using ParameterKey = std::string;
using ParameterSet = std::map<ParameterId, double>;

namespace parameter {
ParameterType get_parameter_type(const ParameterKey& parameter_key);
ParameterId get_parameter_id(const ParameterKey& parameter_key);
ParameterKey create_parameter(double parameter);
extern std::vector<double> parameter_list;
void set_parameter_value(const ParameterKey& parameter_key, double parameter,
    ParameterSet& parameter_set);
void set_parameter_value(const ParameterKey& parameter_key, double parameter);
double get_parameter_value(
    const ParameterKey& parameter_key, const ParameterSet& parameter_set);

}  // namespace parameter
