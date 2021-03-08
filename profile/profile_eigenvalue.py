from memory_profiler import profile, memory_usage
from qulacs import Observable, QuantumState
import argparse

def profile_computing_eigenvalue(method, qubit_count, operator_count, iter_count):
    state = QuantumState(qubit_count)
    state.set_Haar_random_state()
    observable = Observable(qubit_count)
    observable.add_random_operator(operator_count)
    if method == "power":
        observable.solve_ground_state_eigenvalue_by_power_method(state, iter_count)
    elif method == "arnoldi":
        observable.solve_ground_state_eigenvalue_by_arnoldi_method(state, iter_count)
    else:
        print("Invalid method; method should be power or arnoldi.")
        exit(1)

def profile_over_qubit_count(method, qubit_count_start, qubit_count_end, operator_count, iter_count):
    # Holds (qubit_count, memory_usage)
    memory_usage_list = []
    for qubit_count in range(qubit_count_start, qubit_count_end):
        usage_history = memory_usage((profile_computing_eigenvalue, (method, qubit_count, operator_count, iter_count)))
        usage = max(usage_history)
        memory_usage_list.append((qubit_count, usage))
        print("{}MB in {} method with {} qubit".format(usage, method, qubit_count))

    with open("memory_usage_{}_method.csv".format(method), "w") as f:
        f.write("qubit_count,memory_usage\n")
        for usage in memory_usage_list:
            f.write("{},{}\n".format(usage[0], usage[1]))
    return memory_usage_list

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--qubit_count_start", type=int, help="Lower bound of qubit count of observable to profile")
    parser.add_argument("-e", "--qubit_count_end", type=int, help="Upper bound of qubit count of observable to profile")
    parser.add_argument("-o", "--operator_count", type=int, help="Operator count of observable to profile")
    parser.add_argument("-i", "--iter_count", type=int, help="Iteration count for each profile")
    args = parser.parse_args()
    qubit_count_start = args.qubit_count_start
    qubit_count_end = args.qubit_count_end
    operator_count = args.operator_count
    iter_count = args.iter_count

    method = "power"
    profile_over_qubit_count(method, qubit_count_start, qubit_count_end, operator_count, iter_count)

    method = "arnoldi"
    profile_over_qubit_count(method, qubit_count_start, qubit_count_end, operator_count, iter_count)

if __name__ == "__main__":
    main()
