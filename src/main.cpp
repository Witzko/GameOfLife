#include <iostream>
#include <vector>
#include "../include/Matrix.hpp"
#include "../include/Generation.hpp"
#include "../include/functions.hpp"

int main()
{
    Generation first_gen{Matrix(16, 0.8)};
    int alive_cells{0}, dead_cells{0};
    first_gen.getGeneration().countAliveAndDeadCells(alive_cells, dead_cells);
    std::cout << "First generation \n" "Alive Cells: " << alive_cells << " Dead Cells: " << dead_cells << std::endl;

    Generation next_gen = calculateNextGen(first_gen);
    first_gen.printGeneration("first_gen");
    next_gen.printGeneration("next_gen");

    alive_cells = 0;
    dead_cells = 0;
    next_gen.getGeneration().countAliveAndDeadCells(alive_cells, dead_cells);
    std::cout << "Second generation \n" "Alive Cells: " << alive_cells << " Dead Cells: " << dead_cells << std::endl;
    return 0;
}