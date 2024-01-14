import numpy as np
from matplotlib import pyplot as plt
import os
import sys

pwd = os.path.dirname(os.path.realpath(__file__))

path_to_benchmark_1 = pwd + "/parallel_benchmark_strongscaling.txt"
path_to_benchmark_2 = pwd + "/sequential_benchmark_strongscaling.txt"

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

# Runtime plot
# fig, ax = plt.subplots()
# ax.scatter(z, y)
#
# for i, txt in enumerate(n):
#     ax.annotate(txt, (z[i], y[i]))

plt.title("Runtime, gridsize=[{},{}]".format(grid_sizes[0], grid_sizes[1]))
plt.plot(np.array(processes_cols) * np.array(processes_rows), times_parallel, label="ble")
plt.xlabel("Processes")
plt.ylabel("Avg time per iteration")
plt.legend()
plt.savefig("runtime.svg")
plt.close()

# Speedup plot
plt.title("Speedup, gridsize=[{},{}]".format(grid_sizes[0], grid_sizes[1]))
plt.plot(np.array(processes_cols) * np.array(processes_rows), np.array(times_sequential) / np.array(times_parallel), label="ble")
plt.xlabel("Processes")
plt.ylabel("Speedup [seconds]")
plt.legend()
plt.savefig("speedup.svg")
plt.close()

# Parallel efficiency plot
plt.title("Speedup, gridsize=[{},{}]".format(grid_sizes[0], grid_sizes[1]))
plt.plot(np.array(processes_cols) * np.array(processes_rows), (np.array(times_sequential) / np.array(times_parallel)) / (np.array(processes_cols) * np.array(processes_rows)), label="ble")
# plt.plot(np.array(processes_cols) * np.array(processes_rows), np.array(processes_cols) * np.array(processes_rows), "--", label="bleblelbelble")
plt.xlabel("Numbers")
plt.ylabel("Operations per second")
plt.legend()
plt.savefig("parallel_efficiency.svg")
plt.close()
