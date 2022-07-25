このドキュメントはチュートリアル作成のための情報をまとめるために置いてます。プルリクを送る前に削除する予定です。
## 一覧
```
o : qulacs qulacs-osaka両方
+ : qulacs-osakaのみ
- : qulacsのみ

0 : 解説無し
1 : Python教材
2 : 上級編
```
### PauliOperator
- o0 init(coef)
- o2 init(pauli_string, coef)
- +0 init(x_bits, z_bits, coef)
- o2 get_index_list
- o2 get_pauli_id_list
- o2 get_coef
- o2 add_single_Pauli
- o2 get_expectation_value
- +0 get_expectation_value_single_thread
- o2 get_transition_amplitude
- o2 copy
- +0 get_pauli_string
- +0 change_coef
- +0 get_x_bits
- +0 get_z_bits
- +0 self * self
- +0 \_\_mul__
- +0 update_quantum_state
### GeneralQuantumOperator
- o2 init
- o2 add_operator(pauli_operator)
- o2 add_operator(coef, pauli_string)
- o0 is_hermitian
- o2 get_qubit_count
- o0 get_state_dim
- o2 get_term_count
- +0 apply_to_state
- o2 get_term
- o0 get_expectation_value
- +0 get_expectation_value_single_thread
- o2 get_transition_amplitude
- +0 \_\_str__
- +0 copy
- +0 self + self
- +0 \_\_add__
- +0 self += self
- +0 \_\_IADD__
- +0 self - self
- +0 \_\_sub__
- +0 self -= self
- +0 \_\_ISUB__
- +0 self * self
- +0 \_\_mul__(GeneralQuantumOperator, PauliOperator)
- +0 \_\_mul__(PauliOperator, std::complex<double>)
- +0 self * self
- +0 \_\_IMUL__(GeneralQuantumOperator, PauliOperator)
- +0 \_\_IMUL__(PauliOperator, std::complex<double>)
#### quantum_operator
- o2 create_quantum_operator_from_openfermion_file
- o2 create_quantum_operator_from_openfermion_text
- o0 create_split_quantum_operator
### Observable
- o0 init
- o0 add_operator(pauli_operator)
- o0 add_operator(coef, string)
- o0 get_qubit_count
- o0 get_state_dim
- o0 get_term_count
- o0 get_term
- o0 get_expectation_value
- +0 get_expectation_value_single_thread
- o0 get_transition_amplitude
- +0 add_random_operator(operator_count)
- +0 add_random_operator(operator_count, seed)
- +0 solve_ground_state_eigenvalue_by_arnoldi_method
- +0 solve_ground_state_eigenvalue_by_power_method
- +0 solve_ground_state_eigenvalue_by_lanczos_method
- +0 apply_to_state
- +0 \_\_str__
#### observable
- +2 create_observable_from_openfermion_file
- +0 create_observable_from_openfermion_text
- +2 create_split_observable
### QuantumStateBase
### QuantumState
- o1 init
- o1 set_zero_state
- o1 set_computational_basis
- o1 set_Haar_random_state()
- o1 set_Haar_random_state(seed)
- o1 get_zero_probability
- o2 get_marginal_probability
- o2 get_entropy
- o2 get_squared_norm
- o2 normalize
- o2 allocate_buffer
- o0 copy
- o2 load(state : QuantumStateBase)
- o1 load(state : std::vector)
- o2 get_device_name
- o2 add_state
- o2 multiply_coef
- o0 multiply_elementwise_function
- o2 get_classical_value
- o2 set_classical_value
- o0 to_string
- o1 sampling(count)
- o0 sampling(count, seed)
- o1 get_vector
- +0 get_amplitude
- o2 get_qubit_count
- o1 \_\_repr__
- o0 StateVector
### DensityMatrix
- o0 init
- o0 set_zero_state
- o0 set_computational_basis
- o0 set_Haar_random_state()
- o0 set_Haar_random_state(seed)
- o0 get_zero_probability
- o0 get_merginal_probability
- o0 get_entropy
- o0 get_squared_norm
- o0 normalize
- o0 allocate_buffer
- o0 copy
- o0 load(state : QuantumStateBase)
- o0 load(state : std::vector)
- o0 load(state : ComplexMatrix)
- o0 get_device_name
- o0 add_state
- o0 multiply_coef
- o0 get_classical_value
- o0 set_classical_value
- o0 to_string
- o0 sampling(count)
- o0 sampling(count, seed)
- o0 get_matrix
- o0 \_\_repr__
#### state
- o1 inner_product
- o0 tensor_product(QuantumState)
- o0 permutate_qubit(DensityMatrix)
- o0 permutate_qubit(QuantumState)
- o0 permutate_qubit(DensityMatrix)
- o0 drop_qubit
- o0 partial_trace(QuantumState)
- o0 partial_trace(DensityMatrix)
### QuantumGateBase
- o1 update_quantum_state
- o0 copy
- o0 to_string
- o1 get_matrix
- o1 \_\_repr__
- o0 get_target_index_list
- o0 get_control_index_list
- +0 get_control_value_list
- +0 get_control_index_value_list
- o0 get_name
- o0 is_commute
- o0 is_Pauli
- o0 is_Clifford
- o0 is_Gaussian
- o0 is_parametric
- o0 is_diagonal
- +0 get_gate_list
- +0 optimize_ProbabilicGate
- +0 get_distribution
- +0 get?cumulative_distribution
### QuantumGateMatrix
- -0 update_quantum_state
- o1 add_control_qubit
- o2 multiply_scalar
- -0 copy
- -0 to_string
- -0 get_matrix
- -0 \_\_repr__
#### gate
- o2 Identity
- o1 X,Y,Z,H,T
- o2 S,Sdag,Tdag
- o2 sqrtX,sqrtXdag,sqrtY,sqrtYdag
- o2 P0,P1
- o2 U1,U2,U3
- o1 RX,RY,RZ
- o1 CNOT,CZ,SWAP
- o0 TOFFOLI,FREDKIN
- o2 Pauli
- o2 PauliRotation
- o1 DenseMatrix(target_qubit_index, matrix)
- o1 DenseMatrix(target_qubit_index_list, matrix)
- o2 SparseMatrix
- o0 DiagonalMatrix
- o2 RandomUnitary(target_qubit_index_list)
- +0 RandomUnitary(target_qubit_index_list, seed)
- o2 ReversibleBoolean
- o2 StateReflection
- o2 BitFlipNoise
- o2 DephasingNoise
- o2 IndependentXZNoise
- o2 DepolarizingNoise
- o2 TwoQubitDepolarizingNoise
- o2 AmplitudeDampingNoise
- o2 Measurement
- o2 merge(gate1, gate2)
- o0 merge(gate_list)
- o0 add
- o1 to_matrix_gate
- o2 Probabilistic
- o0 ProbabilisticInstrument
- o2 CPTP
- o0 CP
- o2 Instrument
- o2 Adaptive(gate,condition)
- +0 Adaptive(gate,condition,id)
- +0 NoisyEvolution
### QuantumGate_SingleParameter
- o0 get_parameter_value
- o0 set_parameter_value
- o0 copy
### gate
  - o0 ParametricRX,ParametricRY,ParametricRZ
  - o0 ParametricPauliRotation
### QuantumCircuit
- o1 init
- o0 copy
- o1 add_gate(gate)
- o0 add_gate(gate,position)
- +0 add_noise_gate
- o0 remove_gate
- o0 get_gate
- +0 merge_circuit
- o0 get_gate_count
- o0 get_qubit_count
- o1 update_quantum_state(state)
- o0 update_quantum_state(state,start,end)
- o2 calculate_depth
- o0 to_string
- o1 add_X_gate,add_Y_gate,add_Z_gate
- o1 add_H_gate,add_S_gate,add_Sdag_gate,add_T_gate,add_Tdag_gate
- o0 add_sqrtX_gate,add_sqrtXdag_gate,add_sqrtY_gate,add_sqrtYdag_gate
- o0 add_P0_gate,add_P1_gate
- o0 add_CNOT_gate,add_CZ_gate,add_SWAP_gate
- o1 add_RX_gate,add_RY_gate,add_RZ_gate
- o0 add_U1,add_U2,add_U3_gate
- o0 add_multi_Pauli_gate(index_list,pauli_ids)
- o0 add_multi_Pauli_gate(pauli)
- o0 add_multi_Pauli_rotation_gate(index_list,pauli_ids,angle)
- o0 add_multi_Pauli_rotation_gate(pauli)
- o0 add_dense_matrix_gate(index,matrix)
- o0 add_dense_matrix_gate(index_list,matrix)
- o0 add_random_unitary_gate(index_list)
- +0 add_random_unitary_gate(index_list,seed)
- o0 add_diagonal_observable_gate
- o0 add_observable_rotation_gate
- o1 \_\_repr__
### ParametricQuantumCircuit
- o2 init
- o0 copy
- o0 add_parametric_gate(gate)
- o0 add_parametric_gate(gate,position)
- o0 add_gate(gate)
- o0 add_gate(gate,position)
- o2 get_parameter_count
- o2 get_parameter
- o2 set_parameter
- o0 get_parametric_gate_position
- o0 remove_gate
- o2 add_parametric_RX_gate,add_parametric_RY_gate,add_parametric_RZ_gate
- o2 add_parametric_multi_Pauli_rotation_gate
- +0 backprop
- o0 backprop_inner_product
- \_\_repr__
### GradCalculator
- +0 init
- +0 calculate_grad(parametric_circuit,observable)
- +0 calculate_grad(parametric_circuit,observable,angles of gates)
#### circuit
### QuantumCircuitOptimizer
- o2 init
- o2 optimize
- o0 optimize_light
- o0 merge_all
### QuantumCircuitSimulator
- o0 init
- o0 initialize_random_state()
- +0 initialize_random_state(seed)
- o0 simulate
- o0 simulate_range
- o0 get_expectation_value
- o0 get_gate_count
- o0 copy_state_to_buffer
- o0 copy_state_from_buffer
- o0 swap_state_and_buffer
### CausalConeSimulator
- +0 init
- +0 build
- +0 get_expectation_value
- +0 get_circuit_list
- +0 get_pauli_operator_list
- +0 get_coef_list
### NoiseSimulator
- +0 init
- +0 execute

``` plantuml
@startuml
class PauliOperator
class GeneralQuantumOperator
class Observable
GeneralQuantumOperator <|-- Observable
class QuantumStateBase
class QuantumState
QuantumStateBase <|-- QuantumState
class DensityMatrix
QuantumStateBase <|-- DensityMatrix
class QuantumGateBase
class QuantumGateMatrix
QuantumGateBase <|-- QuantumGateMatrix
class QuantumCircuit
class ParametricQuantumCircuit
QuantumCircuit <|-- ParametricQuantumCircuit
class GradCalculator
class QuantumCircuitOptimizer
class QuantumCircuitSimulator
class CausalConeSimulator
class NoiseSimulator
@enduml
```
