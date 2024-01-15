#!/bin/bash

# Run parallel app in debug mode with 2 processes (note --args produces a warning but that's the only
# way GBD seems to find the file symbol inside build/parallel).
# make clean
# make debug parallel
# mpirun -n 4 xterm -hold -e gdb -ex --args ./build/parallel

# In the GBD terminal, for example put:
# $ break <line_number>
# $ run <grid_size> <grid_size> <prob_of_life> <num_iterations> <weak_scaling>
# Example -> run 64 64 0.4 10 false

# Comment everything below if you just want to run a trivial parallel program for testing:
make clean
make debug parallel
mpirun -n 9 ./build/parallel 6 6 0.5 1 false 3 3

# Run below with valgrind
# make clean
# make debug parallel
# mpirun -n 2 valgrind ./build/parallel 4 4 0.5 1 false


TRÄFF STUFF SOURCES0
2 , left
1 , right
3 , down
6 , up
5 , left_down
8 , left_up
7 , right_up
4 , right_down
TRÄFF STUFF DESTINATIONS0
1 , right
2 , left
6 , up
3 , down
7 , right_up
4 , right_down
5 , left_down
8 , left_up

