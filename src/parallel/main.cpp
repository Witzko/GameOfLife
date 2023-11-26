#include <iostream>
#include <vector>
#include <mpi.h>
#include <cassert>
#include "../../include/Matrix.hpp"
#include "../../include/Generation.hpp"
#include "../../include/functions.hpp"

int main(int argc, char **argv)
{

    if (argc != 5)
    {
        std::cerr << "Usage: " << argv[0] << " <matrix_dim_N> <prob_of_life> <num_of_repetitions> <num_of_processes>\n";
        return 1;
    }

    int N = std::atoi(argv[1]);
    float prob_of_life = std::atof(argv[2]);
    int num_of_repetitions = std::atoi(argv[3]);
    int number_of_processes = std::atoi(argv[4]);

    /*
        Initialization:
    */
    Generation current_gen{Matrix(N, prob_of_life)};
    
#ifdef DEBUG
    current_gen.printGeneration("first_gen");
#endif

    Generation next_gen{};

    /*
        MPI Section Start:
    */
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    assert(0 <= rank && rank < size);

    double start_time, end_time;
    std::vector<double> times;

    for (int i = 0; i < num_of_repetitions; i++)
    {
        start_time = MPI_Wtime();

        if (rank == 0)
        {
            next_gen = calculateNextGenParallel(current_gen, number_of_processes);
        }

        end_time = MPI_Wtime();
        times.push_back(end_time - start_time);

#ifdef DEBUG
        Generation next_gen_sequential = calculateNextGenSequentially(current_gen);
        if (!areGenerationsEqual(next_gen, next_gen_sequential))
        {
            next_gen.printGeneration("parallel_gen");
            next_gen_sequential.printGeneration("sequential_gen");
        }
#endif

        current_gen = next_gen;
    }

    MPI_Finalize();
    /*
        MPI Section End

        Post Processing Start
    */

#ifdef DEBUG
    next_gen.printGeneration("last_gen");
#endif

    int alive_cells{0};
    int dead_cells{0};
    countAliveAndDeadCells(current_gen, alive_cells, dead_cells);
    std::cout << "Last generation \n"
                 "Alive Cells: "
              << alive_cells << " Dead Cells: " << dead_cells << std::endl;

    std::cout << "Average calculation time per generation: " << averageVectorElements(times) << std::endl;
    return 0;
}