# HPC example project: Game of Life Stencil

A parallized game of life implementation using MPI, for educational purposes. Using both p2p communication and neighbourhood collectives.

## 1. Project Structure

The project folders are:

    - /include: hpp files
    - /src: cpp files
    - /build: .o files, .exe files
    - /debug: .csv debug files
    - /pics: .png elements for documentation use
    - /plotting: benchmark plotting python scripts etc
    - /benchmark: benchmark raw result data etc
    - /data: .odp file for documentation use, other files used during implementation

Also in the root folder, we have:

    - Makefile
    - run.sh
    - Readme.pdf

## 2. How to build and run it

The project can be build with the use of a Makefile, both in the optimized version and the debug version:

----------------------------------------------------------------------------------------------------------------------
    make clean

    make sequential 
    make parallel

    make debug sequential
    make debug parallel
----------------------------------------------------------------------------------------------------------------------

The object files and executables are then located in the /build folder.

To run the program, we included a shell file **run.sh** for automation. The executable can also be called with the following commands and CL arguments for the sequential version:

----------------------------------------------------------------------------------------------------------------------
    mpirun -n 1 ./build/sequential <matrix_size_row> <matrix_size_col> <prob_of_life> <number_of_repetitions>

----------------------------------------------------------------------------------------------------------------------


And the following arguments for the parallel version:

----------------------------------------------------------------------------------------------------------------------

    mpirun -n <num_of_processes> ./build/parallel <matrix_size_row> <matrix_size_col> <prob_of_life> <number_of_repetitions> <weak_scaling_flag> <num_procs_by_col> <num_procs_by_row>
----------------------------------------------------------------------------------------------------------------------


with the command line arguments:

    - matrix_size_row: int
    - matrix_size_col: int
    - prob_of_life: float in range [0, 1.0]
    - number_of_repetitions: int
    - weak_scaling_flag: true | false
    - num_procs_by_col: int
    - num_procs_by_row: int

e.g. parallel execution: *mpirun -n 16 ./build/parallel 16 16 0.6 100 false 4 4*

This runs the program with 16 processes in a 4x4 grid for a matrix size of 16x16 with probability of life of 60%, 100 iterations and weak-scaling off.

**Note**: We did not include a command line argument to decide if we run it with collective communication or with P2P. The resepective function must be commented/uncommented in line 218 in /src/parallel/main.cpp.

