#include <iostream>
#include <vector>
#include <mpi.h>
#include <cassert>
#include <string>
#include <sstream>
#include "../../include/Generation.hpp"
#include "../../include/functions.hpp"

int main(int argc, char **argv)
{

    /*
        Initialization
    */
    if (argc != 6)
    {
        std::cerr << "Usage: " << argv[0] << " <matrix_size_row> <matrix_size_col> <prob_of_life> <number_of_repetitions> <weak_scaling_flag>\n";
        return 1;
    }

    int row_size = std::atoi(argv[1]);
    int col_size = std::atoi(argv[2]);
    float prob_of_life = std::atof(argv[3]);
    int number_of_repetitions = std::atoi(argv[4]);
    std::string weak_scaling_string = argv[5];
    std::istringstream iss(weak_scaling_string);
    bool weak_scaling_flag;
    if (!(iss >> std::boolalpha >> weak_scaling_flag))
    {
        std::cerr << "Invalid boolean argument: " << weak_scaling_flag << std::endl;
        return 1;
    }
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

    /*
        Initialization of Generation object: Depends on weak/strong scaling flag
    */
    Generation current_gen;

    if (weak_scaling_flag == true)
    {
        current_gen = Generation{row_size * size, col_size * size, prob_of_life}; // if weak scaling, the matrix size of one process is given by num_rows, num_cols and the full matrix size is therefore
                                                                                  // row_size * size (= num of processes), col_size * size
    }
    else
    {
        current_gen = Generation{row_size, col_size, prob_of_life}; // if strong scaling, the matrix size of one process is dynamic. It gets evaluated in the function calculateNextGenParallel
    }

#ifdef DEBUG
    current_gen.printGeneration("first_gen");
#endif

    Generation next_gen = current_gen; // tmp same values as first gen for first calculation

    double start_time, end_time;
    std::vector<double> times;

    MPI_Barrier(cart_comm);
    for (int i = 0; i < number_of_repetitions; i++)
    {
        start_time = MPI_Wtime();

        calculateNextGenParallel(current_gen, next_gen, cart_comm, weak_scaling_flag);

        MPI_Barrier(cart_comm); // Träff said this is kind of okay in the lecture, but he would prefer a solution without two MPI Barriers but just one before the for loop.
                                // Stuff to think about when implementing Exercise 3 + 4 ...
        end_time = MPI_Wtime();
        times.push_back(end_time - start_time);

#ifdef DEBUG
        MPI_Barrier(cart_comm);
        if (rank == 0)
        {
            Generation next_gen_sequential = calculateNextGenSequentially(current_gen);
            if (!areGenerationsEqual(next_gen, next_gen_sequential))
            {
                std::cout << "The sequential and parallel solution are not the same! Check the /debug folder." << std::endl;
                next_gen.printGeneration("parallel_gen");
                next_gen_sequential.printGeneration("sequential_gen");
            }
        }
        MPI_Barrier(cart_comm);
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