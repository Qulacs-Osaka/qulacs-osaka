#include "parameter.hpp"

namespace parameter {
std::vector<double> parameter_list;
ParameterType get_parameter_type(const ParameterKey& parameter_key) {
    auto it = std::find(parameter_key.begin(), parameter_key.end(), ':');
    assert(it != parameter_key.end());
    return std::string(parameter_key.begin(), it);
}
ParameterId get_parameter_id(const ParameterKey& parameter_key) {
    auto it = std::find(parameter_key.begin(), parameter_key.end(), ':');
    assert(it != parameter_key.end());
    return std::string(it + 1, parameter_key.end());
}

ParameterKey create_parameter(double parameter) {
    UINT index = parameter_list.size();
    parameter_list.push_back(parameter);
    return "global:" + std::to_string(index);
}
void set_parameter_value(const ParameterKey& parameter_key, double parameter,
    ParameterSet& parameter_set) {
    ParameterType ptype = get_parameter_type(parameter_key);
    ParameterId pid = get_parameter_id(parameter_key);
    if (ptype == "user") {
        auto it = parameter_set.find(pid);
        if (it != parameter_set.end()) {
            parameter_set[pid] = parameter;
            return;
        }
    }
    if (ptype == "global") {
        UINT pid_int = std::stoi(pid);
        if (pid_int < parameter::parameter_list.size()) {
            parameter::parameter_list[pid_int] = parameter;
            return;
        }
    }
    throw ParameterKeyNotFoundException(
        "Error: parameter_internal::set_parameter_value("
        "const ParameterKey&, double, parameter_set&): parameter_key" +
        parameter_key + "is not found in parameter_set or global");
}
void set_parameter_value(const ParameterKey& parameter_key, double parameter) {
    ParameterSet dummy;
    set_parameter_value(parameter_key, parameter, dummy);
}
double get_parameter_value(
    const ParameterKey& parameter_key, const ParameterSet& parameter_set = {}) {
    ParameterType ptype = get_parameter_type(parameter_key);
    ParameterId pid = get_parameter_id(parameter_key);
    if (ptype == "user") {
        auto it = parameter_set.find(pid);
        if (it != parameter_set.end()) return it->second;
    }
    if (ptype == "global") {
        UINT pid_int = std::stoi(pid);
        if (pid_int < parameter::parameter_list.size()) {
            return parameter::parameter_list[pid_int];
        }
    }
    throw ParameterKeyNotFoundException(
        "Error: parameter_internal::get_parameter_value(const ParameterKey&, "
        "const ParameterSet&): parameter_key" +
        parameter_key + "is not found in parameter_set or global");
}

}  // namespace parameter