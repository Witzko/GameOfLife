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
    executable="/build/sequential"
    make sequential
elif [ "$executable_type" == "parallel" ]; then
    executable="/build/parallel"
    make parallel
else
    echo "Invalid exe type. Supported types are: sequential, parallel"
    exit 1
fi

# Define the parameter values
Ns=(1024 10240)
prob_lives=(0.4 0.6)
repetitions=(10 100)

# Loop over parameter combinations
for N in "${Ns[@]}"; do
    for prob_of_life in "${prob_lives[@]}"; do
        for num_of_repetitions in "${repetitions[@]}"; do
            echo "Running with N=$N, prob_of_life=$prob_of_life, repetitions=$num_of_repetitions"
            mpirun -np 1 $executable $N $prob_of_life $num_of_repetitions
        done
    done
done