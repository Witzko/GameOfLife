#! /bin/bash

#SBATCH -p q_student
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --cpu-freq=High
#SBATCH --time=1:00
#SBATCH --output=benchmark/unformatted_outputs/test_benchmark_ugly.out

executable="./build/parallel"

num_rows=(4)
num_cols=(4)
prob_lives=(0.6)
repetitions=(10)
processes=(4)

for ((i=0; i<${#num_rows[@]}; i++)); do
    N="${num_rows[i]}"
    M="${num_cols[i]}"
    for prob_of_life in "${prob_lives[@]}"; do
        for num_of_repetitions in "${repetitions[@]}"; do
            for num_of_processes in "${processes[@]}"; do
                echo "Running in parallel with num_rows=$N, num_cols=$M, prob_of_life=$prob_of_life, repetitions=$num_of_repetitions, num_processes=$num_of_processes"
                mpirun -np $num_of_processes $executable $N $M $prob_of_life $num_of_repetitions false 1 $num_of_processes
            done
        done
    done
done
