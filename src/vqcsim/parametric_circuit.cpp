#define _USE_MATH_DEFINES
#include "parametric_circuit.hpp"

#include <cppsim/exception.hpp>
#include <cppsim/gate_factory.hpp>
#include <cppsim/gate_matrix.hpp>
#include <cppsim/gate_merge.hpp>
#include <cppsim/state.hpp>
#include <cppsim/type.hpp>
#include <iostream>

#include "parametric_gate.hpp"
#include "parametric_gate_factory.hpp"

ParametricQuantumCircuit::ParametricQuantumCircuit(UINT qubit_count_)
    : QuantumCircuit(qubit_count_){};
ParametricQuantumCircuit::ParametricQuantumCircuit(
    UINT qubit_count_, const ParameterSet& parameter_set)
    : QuantumCircuit(qubit_count_), _parameter_set(parameter_set){};

ParametricQuantumCircuit* ParametricQuantumCircuit::copy() const {
    ParametricQuantumCircuit* new_circuit =
        new ParametricQuantumCircuit(this->qubit_count, this->_parameter_set);
    new_circuit->_next_parameter_index = this->_next_parameter_index;
    new_circuit->_parametric_gate_position =
        std::vector<UINT>(this->_parametric_gate_position.begin(),
            this->_parametric_gate_position.end());
    for (UINT gate_pos = 0; gate_pos < this->gate_list.size(); gate_pos++) {
        new_circuit->add_gate(this->gate_list[gate_pos]->copy());
    }
    std::transform(this->_parametric_gate_position.begin(),
        this->_parametric_gate_position.end(),
        std::back_inserter(new_circuit->_parametric_gate_list),
        [&](UINT index) {
            return dynamic_cast<QuantumGate_SingleParameter*>(
                this->_gate_list[index]);
        });
    return new_circuit;
}

void ParametricQuantumCircuit::add_parametric_gate(
    QuantumGate_SingleParameter* gate) {
    ParameterKey parameter_id = gate->get_parameter_id();
    if (parameter_id.substr(0, 5) == "gate_") {
        UINT parameter_index = std::stoi(parameter_id.substr(5));
        if (parameter_index < gate::internal::initial_angle_list.size()) {
            _parameter_set[parameter_id] =
                gate::internal::initial_angle_list[parameter_index];
        }
    }
    _parametric_gate_position.push_back((UINT)gate_list.size());
    this->add_gate(gate);
    _parametric_gate_list.push_back(gate);
};
void ParametricQuantumCircuit::add_parametric_gate(
    QuantumGate_SingleParameter* gate, UINT index) {
    ParameterKey parameter_id = gate->get_parameter_id();
    if (parameter_id.substr(0, 5) == "gate_") {
        UINT parameter_index = std::stoi(parameter_id.substr(5));
        if (parameter_index < gate::internal::initial_angle_list.size()) {
            _parameter_set[parameter_id] =
                gate::internal::initial_angle_list[parameter_index];
        }
    }
    _parametric_gate_position.push_back(index);
    this->add_gate(gate, index);
    _parametric_gate_list.push_back(gate);
}
void ParametricQuantumCircuit::add_parametric_gate_copy(
    QuantumGate_SingleParameter* gate) {
    ParameterKey parameter_id = gate->get_parameter_id();
    if (parameter_id.substr(0, 5) == "gate_") {
        UINT parameter_index = std::stoi(parameter_id.substr(5));
        if (parameter_index < gate::internal::initial_angle_list.size()) {
            _parameter_set[parameter_id] =
                gate::internal::initial_angle_list[parameter_index];
        }
    }
    _parametric_gate_position.push_back((UINT)gate_list.size());
    QuantumGate_SingleParameter* copied_gate = gate->copy();
    QuantumCircuit::add_gate(copied_gate);
    _parametric_gate_list.push_back(copied_gate);
};
void ParametricQuantumCircuit::add_parametric_gate_copy(
    QuantumGate_SingleParameter* gate, UINT index) {
    ParameterKey parameter_id = gate->get_parameter_id();
    if (parameter_id.substr(0, 5) == "gate_") {
        UINT parameter_index = std::stoi(parameter_id.substr(5));
        if (parameter_index < gate::internal::initial_angle_list.size()) {
            _parameter_set[parameter_id] =
                gate::internal::initial_angle_list[parameter_index];
        }
    }
    for (auto& val : _parametric_gate_position)
        if (val >= index) val++;
    _parametric_gate_position.push_back(index);
    QuantumGate_SingleParameter* copied_gate = gate->copy();
    QuantumCircuit::add_gate(copied_gate, index);
    _parametric_gate_list.push_back(copied_gate);
}
UINT ParametricQuantumCircuit::get_parameter_id_count() const {
    return (UINT)_parameter_set.size();
}
UINT ParametricQuantumCircuit::get_parametric_gate_count() const {
    return (UINT)_parametric_gate_list.size();
}
std::vector<ParameterKey> ParametricQuantumCircuit::get_parameter_id_list()
    const {
    std::vector<ParameterKey> keys;
    keys.reserve(_parameter_set.size());
    std::transform(_parameter_set.begin(), _parameter_set.end(),
        std::back_inserter(keys), [](auto& pair) { return pair.first; });
    return keys;
}
ParameterKey ParametricQuantumCircuit::generate_parameter_id_from_index(
    UINT index) const {
    return "index_" + std::to_string(index);
}
void ParametricQuantumCircuit::create_parameter(
    const ParameterKey& parameter_id, double initial_parameter) {
    if (_parameter_set.count(parameter_id) > 0) {
        throw ParameterIdDuplicatedException(
            "Error: ParametricQuantumCircuit::create_parameter(ParameterKey&): "
            "parameter_id \"" +
            parameter_id + "\" already exists");
    }
    _parameter_set[parameter_id] = initial_parameter;
}
void ParametricQuantumCircuit::remove_parameter(
    const ParameterKey& parameter_id) {
    if (_parameter_set.count(parameter_id) == 0) {
        throw ParameterIdNotFoundException(
            "Error: ParametricQuantumCircuit::create_parameter(ParameterKey&): "
            "parameter_id \"" +
            parameter_id + "\" is not found");
    }
    _parameter_set.erase(parameter_id);
}
bool ParametricQuantumCircuit::contains_parameter(
    const ParameterKey& parameter_id) const {
    return _parameter_set.count(parameter_id) > 0;
}
double ParametricQuantumCircuit::get_parameter(
    const ParameterKey& parameter_id) const {
    auto it = _parameter_set.find(parameter_id);
    if (it == _parameter_set.end()) {
        throw ParameterIdNotFoundException(
            "Error: ParametricQuantumCircuit::create_parameter(ParameterKey&): "
            "parameter_id \"" +
            parameter_id + "\" is not found");
    }
    return it->second;
}
void ParametricQuantumCircuit::set_parameter(
    const ParameterKey& parameter_id, double value) {
    if (_parameter_set.count(parameter_id) == 0) {
        throw ParameterIdNotFoundException(
            "Error: ParametricQuantumCircuit::create_parameter(ParameterKey&): "
            "parameter_id \"" +
            parameter_id + "\" is not found");
        return;
    }
    _parameter_set[parameter_id] = value;
}
ParameterSet ParametricQuantumCircuit::get_parameter_set() const {
    return _parameter_set;
}
void ParametricQuantumCircuit::set_parameter_set(
    const ParameterSet& parameter_set) {
    _parameter_set = parameter_set;
}

std::string ParametricQuantumCircuit::to_string() const {
    std::stringstream os;
    os << QuantumCircuit::to_string();
    os << "*** Parameter Info ***" << std::endl;
    os << "# of parameter: " << this->get_parameter_id_count() << std::endl;
    return os.str();
}

std::ostream& operator<<(
    std::ostream& stream, const ParametricQuantumCircuit& circuit) {
    stream << circuit.to_string();
    return stream;
}
std::ostream& operator<<(
    std::ostream& stream, const ParametricQuantumCircuit* gate) {
    stream << *gate;
    return stream;
}

UINT ParametricQuantumCircuit::get_parametric_gate_position(UINT index) const {
    if (index >= this->_parametric_gate_list.size()) {
        throw GateIndexOutOfRangeException(
            "Error: "
            "ParametricQuantumCircuit::get_parametric_gate_position(UINT): "
            "index is out of range");
    }

    return _parametric_gate_position[index];
}
void ParametricQuantumCircuit::add_gate(QuantumGateBase* gate) {
    QuantumCircuit::add_gate(gate);
}
void ParametricQuantumCircuit::add_gate(QuantumGateBase* gate, UINT index) {
    QuantumCircuit::add_gate(gate, index);
    for (auto& val : _parametric_gate_position)
        if (val >= index) val++;
}
void ParametricQuantumCircuit::add_gate_copy(const QuantumGateBase* gate) {
    QuantumCircuit::add_gate(gate->copy());
}
void ParametricQuantumCircuit::add_gate_copy(
    const QuantumGateBase* gate, UINT index) {
    QuantumCircuit::add_gate(gate->copy(), index);
    for (auto& val : _parametric_gate_position)
        if (val >= index) val++;
}

void ParametricQuantumCircuit::remove_gate(UINT index) {
    auto ite = std::find(_parametric_gate_position.begin(),
        _parametric_gate_position.end(), (unsigned int)index);
    if (ite != _parametric_gate_position.end()) {
        UINT dist = (UINT)std::distance(_parametric_gate_position.begin(), ite);
        _parametric_gate_position.erase(
            _parametric_gate_position.begin() + dist);
        _parametric_gate_list.erase(_parametric_gate_list.begin() + dist);
    }
    QuantumCircuit::remove_gate(index);
    for (auto& val : _parametric_gate_position)
        if (val >= index) val--;
}

void ParametricQuantumCircuit::add_parametric_RX_gate_new_parameter(
    UINT target_index, const ParameterKey& parameter_id,
    double initial_parameter, double parameter_coef) {
    this->create_parameter(parameter_id, initial_parameter);
    this->add_parametric_gate(
        gate::ParametricRX(target_index, parameter_id, parameter_coef));
}
void ParametricQuantumCircuit::add_parametric_RX_gate_share_parameter(
    UINT target_index, const ParameterKey& parameter_id,
    double parameter_coef) {
    this->add_parametric_gate(
        gate::ParametricRX(target_index, parameter_id, parameter_coef));
}
void ParametricQuantumCircuit::add_parametric_RY_gate_new_parameter(
    UINT target_index, const ParameterKey& parameter_id,
    double initial_parameter, double parameter_coef) {
    this->create_parameter(parameter_id, initial_parameter);
    this->add_parametric_gate(
        gate::ParametricRY(target_index, parameter_id, parameter_coef));
}
void ParametricQuantumCircuit::add_parametric_RY_gate_share_parameter(
    UINT target_index, const ParameterKey& parameter_id,
    double parameter_coef) {
    this->add_parametric_gate(
        gate::ParametricRY(target_index, parameter_id, parameter_coef));
}
void ParametricQuantumCircuit::add_parametric_RZ_gate_new_parameter(
    UINT target_index, const ParameterKey& parameter_id,
    double initial_parameter, double parameter_coef) {
    this->create_parameter(parameter_id, initial_parameter);
    this->add_parametric_gate(
        gate::ParametricRZ(target_index, parameter_id, parameter_coef));
}
void ParametricQuantumCircuit::add_parametric_RZ_gate_share_parameter(
    UINT target_index, const ParameterKey& parameter_id,
    double parameter_coef) {
    this->add_parametric_gate(
        gate::ParametricRZ(target_index, parameter_id, parameter_coef));
}

void ParametricQuantumCircuit::
    add_parametric_multi_Pauli_rotation_gate_new_parameter(
        std::vector<UINT> target, std::vector<UINT> pauli_id,
        const ParameterKey& parameter_id, double initial_angle,
        double parameter_coef) {
    this->create_parameter(parameter_id, initial_angle);
    this->add_parametric_gate(gate::ParametricPauliRotation(
        target, pauli_id, parameter_id, parameter_coef));
}
void ParametricQuantumCircuit::
    add_parametric_multi_Pauli_rotation_gate_share_parameter(
        std::vector<UINT> target, std::vector<UINT> pauli_id,
        const ParameterKey& parameter_id, double parameter_coef) {
    this->add_parametric_gate(gate::ParametricPauliRotation(
        target, pauli_id, parameter_id, parameter_coef));
}

std::vector<double> ParametricQuantumCircuit::backprop_inner_product(
    QuantumState* bistate) {
    // circuitを実行した状態とbistateの、inner_productを取った結果を「値」として、それを逆誤差伝搬します
    // bistateはノルムが1のやつでなくてもよい
    int n = this->qubit_count;
    QuantumState* state = new QuantumState(n);
    //これは、ゲートを前から適用したときの状態を示す
    state->set_zero_state();
    this->update_quantum_state(state);  //一度最後までする

    int num_gates = this->gate_list.size();
    std::vector<int> inverse_parametric_gate_position(num_gates, -1);
    for (UINT i = 0; i < this->get_parametric_gate_count(); i++) {
        inverse_parametric_gate_position[this->_parametric_gate_position[i]] =
            i;
    }
    std::vector<double> ans(this->get_parametric_gate_count());

    /*
    現在、2番のゲートを見ているとする
    ゲート 0 1 2 3 4 5
         state | bistate
    前から2番までのゲートを適用した状態がstate
    最後の微分値から逆算して3番までgateの逆行列を掛けたのがbistate

    1番まで掛けて、YのΘ微分した行列を掛けたやつと、bistateの内積の実数部分をとれば答えが出ることが知られている(知られてないかも)

    ParametricR? の微分値を計算した行列は、Θに180°を足した行列/2 と等しい

    だから、2番まで掛けて、 R?(π) を掛けたやつと、bistateの内積を取る

    さらに、見るゲートは逆順である。
    だから、最初にstateを最後までやって、ゲートを進めるたびにstateに逆行列を掛けている
    さらに、bistateが複素共役になっていることを忘れると、bistateに転置行列を掛ける必要がある。
    しかしこのプログラムではbistateはずっと複素共役なので、転置して共役な行列を掛ける必要がある。
    ユニタリ性より、転置して共役な行列 = 逆行列
    なので、両社にadjoint_gateを掛けている
    */
    QuantumState* Astate = new QuantumState(n);  //一時的なやつ
    for (int i = num_gates - 1; i >= 0; i--) {
        QuantumGateBase* gate_now = this->gate_list[i];  // sono gate
        if (inverse_parametric_gate_position[i] != -1) {
            Astate->load(state);
            QuantumGateBase* RcPI;
            if (gate_now->get_name() == "ParametricRX") {
                RcPI = gate::RX(gate_now->get_target_index_list()[0], M_PI);
            } else if (gate_now->get_name() == "ParametricRY") {
                RcPI = gate::RY(gate_now->get_target_index_list()[0], M_PI);
            } else if (gate_now->get_name() == "ParametricRZ") {
                RcPI = gate::RZ(gate_now->get_target_index_list()[0],
                    M_PI);  // 本当はここで2で割りたいけど、行列を割るのは実装が面倒
            } else if (gate_now->get_name() == "ParametricPauliRotation") {
                ClsParametricPauliRotationGate* pauli_gate_now =
                    (ClsParametricPauliRotationGate*)gate_now;
                RcPI =
                    gate::PauliRotation(pauli_gate_now->get_target_index_list(),
                        pauli_gate_now->get_pauli()->get_pauli_id_list(), M_PI);
            } else {
                std::stringstream error_message_stream;
                error_message_stream
                    << "Error: " << gate_now->get_name()
                    << " does not support backprop in parametric";
                throw NotImplementedException(error_message_stream.str());
            }
            RcPI->update_quantum_state(Astate);
            ans[inverse_parametric_gate_position[i]] =
                state::inner_product(bistate, Astate).real() /
                2.0;  //だからここで2で割る
            delete RcPI;
        }
        auto Agate = gate::get_adjoint_gate(gate_now);
        Agate->update_quantum_state(bistate);
        Agate->update_quantum_state(state);
        delete Agate;
    }
    delete Astate;
    delete state;
    return ans;
}  // CPP

std::vector<double> ParametricQuantumCircuit::backprop(
    GeneralQuantumOperator* obs) {
    //オブザーバブルから、最終段階での微分値を求めて、backprop_from_stateに流す関数
    //上側から来た変動量 * 下側の対応する微分値 =
    //最終的な変動量になるようにする。

    int n = this->qubit_count;
    QuantumState* state = new QuantumState(n);
    state->set_zero_state();
    this->update_quantum_state(state);  //一度最後までする
    QuantumState* bistate = new QuantumState(n);
    QuantumState* Astate = new QuantumState(n);  //一時的なやつ

    obs->apply_to_state(Astate, *state, bistate);
    bistate->multiply_coef(2);
    /*一度stateを最後まで求めてから、さらにapply_to_state している。
    なぜなら、量子のオブザーバブルは普通の機械学習と違って、二乗した値の絶対値が観測値になる。
    二乗の絶対値を微分したやつと、値の複素共役*2は等しい

    オブザーバブルよくわからないけど、テストしたらできてた
    */

    //ニューラルネットワークのbackpropにおける、後ろからの微分値的な役目を果たす
    auto ans = backprop_inner_product(bistate);
    delete bistate;
    delete state;
    delete Astate;
    return ans;

}  // CPP

// 旧メソッド
double ParametricQuantumCircuit::get_parameter(UINT index) const {
    return this->get_parameter(this->generate_parameter_id_from_index(index));
}
void ParametricQuantumCircuit::set_parameter(UINT index, double value) {
    return this->set_parameter(
        this->generate_parameter_id_from_index(index), value);
}
void ParametricQuantumCircuit::add_parametric_RX_gate(
    UINT target_index, double initial_angle) {
    std::string parameter_id =
        this->generate_parameter_id_from_index(_next_parameter_index++);
    this->add_parametric_RX_gate_new_parameter(
        target_index, parameter_id, initial_angle);
}
void ParametricQuantumCircuit::add_parametric_RY_gate(
    UINT target_index, double initial_angle) {
    std::string parameter_id =
        this->generate_parameter_id_from_index(_next_parameter_index++);
    this->add_parametric_RY_gate_new_parameter(
        target_index, parameter_id, initial_angle);
}
void ParametricQuantumCircuit::add_parametric_RZ_gate(
    UINT target_index, double initial_angle) {
    std::string parameter_id =
        this->generate_parameter_id_from_index(_next_parameter_index++);
    this->add_parametric_RZ_gate_new_parameter(
        target_index, parameter_id, initial_angle);
}
void ParametricQuantumCircuit::add_parametric_multi_Pauli_rotation_gate(
    std::vector<UINT> target, std::vector<UINT> pauli_id,
    double initial_angle) {
    std::string parameter_id =
        this->generate_parameter_id_from_index(_next_parameter_index++);
    this->add_parametric_multi_Pauli_rotation_gate_new_parameter(
        target, pauli_id, parameter_id, initial_angle);
}
UINT ParametricQuantumCircuit::get_parameter_count() const {
    throw NotImplementedException(
        "Error: ParametricQuantumCircuit::get_parameter_count(): This function "
        "is no longer supported. Use get_parameter_id_count() or "
        "get_parametric_gate_count() instead.");
}