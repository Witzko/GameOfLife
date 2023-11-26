#include <iostream>
#include <vector>
#include <mpi.h>
#include <cassert>
#include "../../include/Matrix.hpp"
#include "../../include/Generation.hpp"
#include "../../include/functions.hpp"

int main(int argc, char **argv)
{

    /*
        Initialization
    */
    if (argc != 4)
    {
        std::cerr << "Usage: " << argv[0] << " <matrix_dim_N> <prob_of_life> <num_of_repetitions>\n";
        return 1;
    }

    int N = std::atoi(argv[1]);
    float prob_of_life = std::atof(argv[2]);
    int num_of_repetitions = std::atoi(argv[3]);

    Generation current_gen{Matrix(N, prob_of_life)};

#ifdef DEBUG
    current_gen.printGeneration("first_gen");
#endif

    Generation next_gen = current_gen; //tmp same values as first gen for first calculation

    /*
        MPI Section Start:
    */
    MPI_Init(&argc, &argv);

    int rank, size;
    const int ndim = 2;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    assert(0 <= rank && rank < size);

    int dims[ndim] = {0, 0};
    MPI_Dims_create(size, 2, dims);

    MPI_Comm cart_comm;
    int periods[ndim] = {1, 1};
    int reorder = 1;
    MPI_Cart_create(MPI_COMM_WORLD, ndim, dims, periods, reorder, &cart_comm);
    assert(cart_comm != MPI_COMM_NULL);

    int coords[ndim];
    MPI_Cart_coords(cart_comm, rank, ndim, coords);

#ifdef DEBUG
    for (int i = 0; i < size; ++i)
    {
        MPI_Barrier(cart_comm);
        if (rank == i)
        {
            std::cout << "Rank " << rank << " has coordinates (" << coords[0] << ", " << coords[1] << ")" << std::endl;
        }
        MPI_Barrier(cart_comm);
    }
#endif

        double start_time, end_time;
        std::vector<double> times;

        for (int i = 0; i < num_of_repetitions; i++)
        {
            start_time = MPI_Wtime();

            calculateNextGenParallel(current_gen, next_gen, cart_comm);

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

    MPI_Comm_free(&cart_comm);
    MPI_Finalize();

    /*
        MPI Section End

        POST PROCESSING Section Start
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

    std::cout << "Average calculation time in micro seconds per generation: " << averageVectorElements(times) << std::endl;
    /*
        POST PROCESSING Section End
    */
    return 0;
}