#! /bin/bash

#SBATCH -p q_student
#SBATCH -N 32
#SBATCH --ntasks-per-node=32
#SBATCH --cpu-freq=High
#SBATCH --time=10:00
#SBATCH --output=benchmark/unformatted_outputs/parallel_benchmark_strongscaling_alltoall_ugly.out

executable="./build/parallel"

num_rows=(1024 10240)
num_cols=(1024 10240)
prob_lives=(0.4 0.6)
repetitions=(10)
processes=(1 8 16 32)

for ((i=0; i<${#num_rows[@]}; i++)); do
    N="${num_rows[i]}"
    M="${num_cols[i]}"
    for prob_of_life in "${prob_lives[@]}"; do
        for num_of_repetitions in "${repetitions[@]}"; do
            for num_of_processes in "${processes[@]}"; do
                echo "Running in parallel with num_rows=$N, num_cols=$M, prob_of_life=$prob_of_life, repetitions=$num_of_repetitions, num_proc_x=32, num_proc_y=$num_of_processes"
                mpirun -np $(($num_of_processes*32)) $executable $N $M $prob_of_life $num_of_repetitions false $num_of_processes 32
            done
        done
    done
done
