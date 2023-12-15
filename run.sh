#!/bin/bash

# Check if an argument is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <executable_type>"
    exit 1
fi

executable_type="$1"
executable=""

# Set the executable based on the provided argument
if [ "$executable_type" == "sequential" ]; then
    executable="./build/sequential"
    make clean
    make sequential
elif [ "$executable_type" == "parallel" ]; then
    executable="./build/parallel"
    make clean
    make parallel
else
    echo "Invalid exe type. Supported types are: sequential, parallel"
    exit 1
fi

# Define the parameter values
num_rows=(64 1024)
num_cols=(64 1024)
prob_lives=(0.4 0.6)
repetitions=(1 10)
processes=(2 4)

# Loop over parameter combinations
for ((i=0; i<${#num_rows[@]}; i++)); do
    N="${num_rows[i]}"
    M="${num_cols[i]}"
    for prob_of_life in "${prob_lives[@]}"; do
        for num_of_repetitions in "${repetitions[@]}"; do
            if [ "$executable" == "./build/sequential" ]; then
                echo "Running sequentially with num_rows=$N, num_cols=$M, prob_of_life=$prob_of_life, repetitions=$num_of_repetitions"
                mpirun -np 1 $executable $N $M $prob_of_life $num_of_repetitions
            else
                for num_of_processes in "${processes[@]}"; do
                    echo "Running in parallel with num_rows=$N, num_cols=$M, prob_of_life=$prob_of_life, repetitions=$num_of_repetitions, processes=$num_of_processes"
                    mpirun -np $num_of_processes $executable $N $M $prob_of_life $num_of_repetitions false
                done
            fi
        done
    done
done
