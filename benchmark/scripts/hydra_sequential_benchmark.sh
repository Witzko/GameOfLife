#! /bin/bash

#SBATCH -p q_student
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpu-freq=High
#SBATCH --time=10:00
#SBATCH --output=benchmark/unformatted_outputs/sequential_benchmark_ugly.out

executable="./build/sequential"

num_rows=(1024 10240)
num_cols=(1024 10240)
prob_lives=(0.4 0.6)
repetitions=(10)

for ((i=0; i<${#num_rows[@]}; i++)); do
    N="${num_rows[i]}"
    M="${num_cols[i]}"
    for prob_of_life in "${prob_lives[@]}"; do
        for num_of_repetitions in "${repetitions[@]}"; do
            echo "Running sequentially with num_rows=$N, num_cols=$M, prob_of_life=$prob_of_life, repetitions=$num_of_repetitions"
            mpirun -np 1 $executable $N $M $prob_of_life $num_of_repetitions
        done
    done
done
