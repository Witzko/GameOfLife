#include <iostream>
#include <vector>
#include <mpi.h>
#include "../../include/Matrix.hpp"
#include "../../include/Generation.hpp"
#include "../../include/functions.hpp"

int main(int argc, char **argv)
{

    /*
        INITIALIZATION:
        initialize first Generation with size NxN, liveness probability for each cell

    */
    Generation first_gen{Matrix(16, 0.5)};
    first_gen.printGeneration("first_gen");

    int alive_cells{0};
    int dead_cells{0};
    countAliveAndDeadCells(first_gen, alive_cells, dead_cells);
    std::cout << "First generation \n"
                 "Alive Cells: "
              << alive_cells << " Dead Cells: " << dead_cells << std::endl;


    /* 
        MPI Section Start
        calculate next generation in a sequential MPI program
    
    */
    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0)
    {
    Generation next_gen = calculateNextGenSequentially(first_gen);
    next_gen.printGeneration("next_gen");

    alive_cells = 0;
    dead_cells = 0;
    countAliveAndDeadCells(next_gen, alive_cells, dead_cells);
    std::cout << "Second generation \n"
                 "Alive Cells: "
              << alive_cells << " Dead Cells: " << dead_cells << std::endl;
    }

    MPI_Finalize();
    return 0;
}