#include <iostream>
#include <vector>
#include <mpi.h>
#include <numeric>
#include <cassert>
#include <string>
#include <sstream>
#include <unistd.h>
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
    Generation current_gen, next_gen;

    /*
        If weak scaling, then the dimensions provided should be divided by the amount of processors in
        each x-y direction based on the cartesian communicator domain.
    */
    if (!weak_scaling_flag)
    {
        row_size = row_size / dims[0];
        col_size = col_size / dims[1];
    }

    current_gen = Generation(row_size, col_size, prob_of_life);

    double start_time, end_time;
    std::vector<double> times;

    /*
        Define the MPI_CELL datatype which represents the Cell class for the MPI communication.
    */
    MPI_Datatype MPI_CELL;
    MPI_Datatype types[1] = {MPI_CHAR};
    int blocklength[1] = {1};
    MPI_Aint displacement[1] = {0};
    MPI_Type_create_struct(1, blocklength, displacement, types, &MPI_CELL);
    MPI_Type_commit(&MPI_CELL);

    /*
        Define the MPI vector padding datatype for the MPI communication of the left and right borders (for indexing).
    */
    const int halo_layer_size = 1; // 2024-01-03 -> For now only size 1 works anyway.
    int col_size_whalo = col_size + 2 * halo_layer_size; // As the column is added by one layer on each side

    MPI_Datatype MPI_COL_PADDING_WHALO;
    MPI_Type_vector(row_size, 1, col_size_whalo, MPI_CELL, &MPI_COL_PADDING_WHALO);
    MPI_Type_commit(&MPI_COL_PADDING_WHALO);

    /*
        Create the custom datatype used for the MPI gather of all sub-grids into one global grid.
    */
    std::vector<Cell> global_grid;
    MPI_Datatype recvGridBlock, recvGlobalGridBlock;

    int global_row_size = row_size * dims[0];
    int global_col_size = col_size * dims[1];
    int global_sizes_direction[2] = {global_row_size, global_col_size};
    int local_sizes_direction[2] = {row_size, col_size};
    int local_start_indexes[2] = {0,0};

    /*
        Create MPI datatype for the master receiving all sub-grids and combining them into the global grid.
    */
    if (rank == 0)
    {
        global_grid = std::vector<Cell>(global_row_size * global_col_size, Cell('d'));

        /*
            Create a datatype of the receive block for all the sub-grids
        */
        MPI_Type_create_subarray(ndim, global_sizes_direction, local_sizes_direction, local_start_indexes,
                                 MPI_ORDER_C, MPI_CELL, &recvGridBlock);

        /*
            Resize the recvGridBlock so that when col_size elements have been recevied/placed (1 block) by the master, start
            to place the next rank/senders elements into the next block.

            For example a 6x6 grid of 9 processors
            |block1|block2|block3|
            |--------------------|
            | A  A | B  B | C  C |
            | A  A | B  B | C  C |
            |--------------------|
            | D  D | E  E | F  F |
            | D  D | E  E | F  F |
            |--------------------|
            | G  G | H  H | I  I |
            | G  G | H  H | I  I |
            |--------------------|
        */
        int block_size = col_size * sizeof(Cell);
        MPI_Type_create_resized(recvGridBlock, 0, block_size, &recvGlobalGridBlock);
        MPI_Type_commit(&recvGlobalGridBlock);
    }

    /*
        Each process handles one single block.
    */
    int num_blocks_per_proc[size];
    for (int i = 0; i < size; i++) {
        num_blocks_per_proc[i] = 1;
    }

    /*
        col_padding_block contains the padding (num of col_size) of each sub-grid/block in the global grid.

        A 6x6 grid of 9 processes contains the following col_padding_block:
        col_padding_block[9] = {0, 1, 2, 6, 7, 8, 12, 13, 14}

       [0]    [1]    [2]
        | A  A | B  B | C  C |
        | A  A | B  B | C  C |
       [3]    [4]    [5]
        | D  D | E  E | F  F |
        | D  D | E  E | F  F |
       [6]    [7]    [8]
        | G  G | H  H | I  I |
        | G  G | H  H | I  I |
    */
    int col_padding_block[size];
    int k = 0;
    for (int i = 0; i < dims[0]; i++)
    {
        for (int j = 0; j < dims[1]; j++)
        {
            col_padding_block[k++] = j  +  i * (row_size * dims[1]);
        }
    }

    /*
        Collect the global grid using Gatherv
    */
    MPI_Gatherv(&current_gen.getCell(0,0), row_size * col_size, MPI_CELL,
                &global_grid[0], num_blocks_per_proc, col_padding_block, recvGlobalGridBlock, 0,
                cart_comm);

#ifdef DEBUG
    MPI_Barrier(cart_comm);
    if (rank == 0){
        printf("\n--ORIGINAL GLOBAL GRID--\n");
        printGrid(global_grid, global_row_size, global_col_size);
    }
    MPI_Barrier(cart_comm);
#endif

    /*
        Start the algorithm
     */
    std::vector<Cell> next_gen_cells, curr_gen_cells;
    MPI_Barrier(cart_comm);
    for (int i = 0; i < number_of_repetitions; i++)
    {
        start_time = MPI_Wtime();

        /*
            For Peter boi. You can just comment/uncomment for now.
         */
        next_gen = calculateNextGenParallel(std::move(current_gen), cart_comm, MPI_CELL, MPI_COL_PADDING_WHALO, halo_layer_size);
//        next_gen = calculateNextGenParallelWCollNeighbourComm(current_gen, cart_comm);

        end_time = MPI_Wtime();
        times.push_back(end_time - start_time);

#ifdef DEBUG
        MPI_Barrier(cart_comm);

        if (rank == 0) {
            curr_gen_cells = global_grid; // As the global grid contains the global current_gen entries
        }

        MPI_Gatherv(&next_gen.getCell(0,0), row_size * col_size, MPI_CELL,
                    &global_grid[0], num_blocks_per_proc, col_padding_block, recvGlobalGridBlock, 0,
                    cart_comm);
        /*
            Print global grid and compare it with sequential version
         */
        if (rank == 0){
            printf("\n--GLOBAL GRID--\n");
            printf("Parallel Iteration: %d\n", i);
            printGrid(global_grid, global_row_size, global_col_size);

            next_gen_cells = global_grid;

            Generation global_next_gen(std::move(next_gen_cells), global_row_size, global_col_size);
            Generation global_curr_gen(std::move(curr_gen_cells), global_row_size, global_col_size);

            Generation next_gen_sequential = calculateNextGenSequentially(global_curr_gen);

            if (!areGenerationsEqual(global_next_gen, next_gen_sequential))
            {
                std::cout << "The sequential and parallel solution are not the same! Check the /debug folder." << std::endl;
                global_next_gen.printGeneration("parallel_gen");
                next_gen_sequential.printGeneration("sequential_gen");
            }
        }
        MPI_Barrier(cart_comm);
#endif

        current_gen = next_gen;
    }

    /*
        MPI Section End

        POST PROCESSING Section Start
    */
    int alive_cells{0};
    int dead_cells{0};
    countAliveAndDeadCells(current_gen, alive_cells, dead_cells);

    /*
        Find the total number of alive and dead cells
    */
    int global_alive_cells, global_dead_cells;
    MPI_Reduce(&alive_cells, &global_alive_cells, 1, MPI_INT, MPI_SUM, 0, cart_comm);
    MPI_Reduce(&dead_cells, &global_dead_cells, 1, MPI_INT, MPI_SUM, 0, cart_comm);

    /*
        Find the time of slowest processor as well as the average time per iteration.
    */
    double global_average_time, global_total_time;
    double average_time = averageVectorElements(times);
    double total_time = std::reduce(times.begin(), times.end());

    // Find the average time per repetition
    MPI_Reduce(&average_time, &global_average_time, 1, MPI_DOUBLE, MPI_SUM, 0, cart_comm);
    global_average_time = global_average_time / size;

    // Find the maximum time (i.e. the slowest processor)
    MPI_Reduce(&total_time, &global_total_time, 1, MPI_DOUBLE, MPI_MAX, 0, cart_comm);

    if (rank == 0) {
        std::cout << "Total alive cells: " << global_alive_cells
                  << ", Total dead cells: " << global_dead_cells << std::endl;

        std::cout << "Average time per iteration/repetition: " << averageVectorElements(times)
                  << ", Total time taken (in parallel app the finish time for slowest proc): " << global_total_time << std::endl;
    }

    MPI_Type_free(&MPI_CELL);
    MPI_Type_free(&MPI_COL_PADDING_WHALO);
    if (rank == 0) {
        MPI_Type_free(&recvGlobalGridBlock);
        MPI_Type_free(&recvGridBlock);
    }

    MPI_Comm_free(&cart_comm);
    MPI_Finalize();

    /*
        POST PROCESSING Section End
    */
    return 0;
}
