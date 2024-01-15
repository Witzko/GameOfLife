import numpy as np
from matplotlib import pyplot as plt
import os

pwd = os.path.dirname(os.path.realpath(__file__))

path_to_benchmark_1 = pwd + "/../benchmark/parallel_benchmark_weakscaling.txt"
path_to_benchmark_2 = pwd + "/../benchmark/sequential_benchmark.txt"

processes_rows = []
processes_cols = []
times_parallel = []
times_sequential = []

with open(path_to_benchmark_1, "r") as input:
    grid_sizes = input.readline().strip().split()
    for line in input.readlines():
        values = line.strip().split()
        processes_rows.append(int(values[0]))
        processes_cols.append(int(values[1]))
        times_parallel.append(float(values[2]))

with open(path_to_benchmark_2, "r") as input:
    grid_sizes = input.readline().strip().split()
    for line in input.readlines():
        values = line.strip().split()
        times_sequential.append(float(values[0]))

tot_processes = np.array(processes_cols) * np.array(processes_rows)
n = ["({},{})".format(val[0], val[1]) for _, val in enumerate(zip(processes_rows, processes_cols))]

# Runtime plot
plt.title("Weak scaling runtime, gridsize=[{},{}]".format(grid_sizes[0], grid_sizes[1]))
plt.plot(np.array(processes_cols) * np.array(processes_rows), times_parallel, "-o")
for i, txt in enumerate(n):
    plt.annotate(txt, (tot_processes[i], times_parallel[i]))
plt.xlabel("Processes")
plt.ylabel("Avg time per iteration")
plt.savefig("plots/weak_scaling/runtime_{}.svg".format(grid_sizes[0]))
plt.close()
