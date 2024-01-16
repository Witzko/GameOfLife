import numpy as np
from matplotlib import pyplot as plt
import os

pwd = os.path.dirname(os.path.realpath(__file__))

path_to_benchmark_parallel_weakscaling = pwd + "/../benchmark/parallel_benchmark_weakscaling_alltoall.txt"

processes_rows = []
processes_cols = []
times_parallel = []

with open(path_to_benchmark_parallel_weakscaling, "r") as input:
    grid_sizes = input.readline().strip().split()
    for line in input.readlines():
        values = line.strip().split()
        processes_rows.append(int(values[0]))
        processes_cols.append(int(values[1]))
        times_parallel.append(float(values[2]))

tot_processes = np.array(processes_cols) * np.array(processes_rows)
n = ["({},{})".format(val[0], val[1]) for _, val in enumerate(zip(processes_rows, processes_cols))]

# Runtime plot
plt.title("Weak scaling runtime alltoall, gridsize=[{},{}]".format(grid_sizes[0], grid_sizes[1]))
plt.plot(tot_processes, times_parallel, "-o")
for i, txt in enumerate(n):
    plt.annotate(txt, (tot_processes[i], times_parallel[i]))
plt.xlabel("Processes")
plt.ylabel("Avg time per iteration [seconds]")
plt.savefig("plots/alltoall_weak_scaling/runtime_{}.svg".format(grid_sizes[0]))
plt.close()
