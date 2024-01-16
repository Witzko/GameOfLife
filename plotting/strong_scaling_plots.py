import numpy as np
from matplotlib import pyplot as plt
import os

pwd = os.path.dirname(os.path.realpath(__file__))

path_to_benchmark_parallel_strongscaling = pwd + "/../benchmark/parallel_benchmark_strongscaling_alltoall.txt"
path_to_benchmark_sequential = pwd + "/../benchmark/sequential_benchmark.txt"

processes_rows = []
processes_cols = []
times_parallel = []
time_sequential = .0

with open(path_to_benchmark_parallel_strongscaling, "r") as input:
    grid_sizes = input.readline().strip().split()
    for line in input.readlines():
        values = line.strip().split()
        processes_rows.append(int(values[0]))
        processes_cols.append(int(values[1]))
        times_parallel.append(float(values[2]))

with open(path_to_benchmark_sequential, "r") as input:
    grid_sizes = input.readline().strip().split()
    time_sequential = float(input.readline().strip().split()[0])

tot_processes = np.array(processes_cols) * np.array(processes_rows)
n = ["({},{})".format(val[0], val[1]) for _, val in enumerate(zip(processes_rows, processes_cols))]

# Runtime plot
plt.title("Strong scaling runtime alltoall, gridsize=[{},{}]".format(grid_sizes[0], grid_sizes[1]))
plt.plot(tot_processes, times_parallel, "-o")
for i, txt in enumerate(n):
    plt.annotate(txt, (tot_processes[i], times_parallel[i]))
plt.xlabel("Processes")
plt.ylabel("Avg time per iteration [seconds]")
plt.savefig("plots/alltoall_strong_scaling/runtime_{}.svg".format(grid_sizes[0]))
plt.close()

# Speedup plot
speedup = time_sequential / np.array(times_parallel)
plt.title("Speedup alltoall, gridsize=[{},{}]".format(grid_sizes[0], grid_sizes[1]))
plt.plot(tot_processes, speedup, "-o", label="Speedup")
plt.plot(tot_processes, tot_processes, "-o", label="Linear speedup")
for i, txt in enumerate(n):
    plt.annotate(txt, (tot_processes[i], speedup[i]))
    plt.annotate(txt, (tot_processes[i], tot_processes[i]))
plt.xlabel("Processes")
plt.ylabel("Speedup [$T_1 / T_p$]")
plt.legend()
plt.savefig("plots/alltoall_strong_scaling/speedup_{}.svg".format(grid_sizes[0]))
plt.close()

# Parallel efficiency plot
parallel_efficiency = speedup / tot_processes
plt.title("Parallel efficiency alltoall, gridsize=[{},{}]".format(grid_sizes[0], grid_sizes[1]))
plt.plot(tot_processes, parallel_efficiency, "-o", label="Parallel efficency")
plt.plot(tot_processes, tot_processes / tot_processes, "-o", label="Optimal parallel efficiency")
for i, txt in enumerate(n):
    plt.annotate(txt, (tot_processes[i], parallel_efficiency[i]))
    plt.annotate(txt, (tot_processes[i], tot_processes[i]))
plt.xlabel("Processes")
plt.ylabel("Parallel efficiency [speedup / num_proc]")
plt.legend()
plt.savefig("plots/alltoall_strong_scaling/parallel_efficiency_{}.svg".format(grid_sizes[0]))
plt.close()
