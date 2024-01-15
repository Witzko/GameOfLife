#!/bin/bash

# Run parallel app in debug mode with 2 processes (note --args produces a warning but that's the only
# way GBD seems to find the file symbol inside build/parallel).
# make clean
# make debug parallel
# mpirun -n 2 xterm -hold -e gdb -ex --args ./build/parallel

# In the GBD terminal, for example put:
# $ break <line_number>
# $ run <grid_size> <grid_size> <prob_of_life> <num_iterations> <weak_scaling>
# Example -> run 64 64 0.4 10 false

# Comment everything below if you just want to run a trivial parallel program for testing:
make clean
make debug parallel
mpirun -n 2 ./build/parallel 4 4 0.5 1 false 1 2

# Run below with valgrind
# make clean
# make debug parallel
# mpirun -n 2 valgrind ./build/parallel 4 4 0.5 1 false
