#include <iostream>
#include <vector>
#include "../include/Matrix.hpp"
#include "../include/Generation.hpp"

int main()
{
    Generation first_gen{std::make_unique<Matrix>(16, 0.8)};
    int alive_cells{0}, dead_cells{0};
    first_gen.getCurrentGen().countAliveAndDeadCells(alive_cells, dead_cells);
    // std::cout << "First generation \n" "Alive Cells: " << alive_cells << " Dead Cells: " << dead_cells << std::endl;

    first_gen.calculateNextGen();
    first_gen.printGenerations();
    alive_cells = 0;
    dead_cells = 0;
    first_gen.getNextGen().countAliveAndDeadCells(alive_cells, dead_cells);
    // std::cout << "Second generation \n" "Alive Cells: " << alive_cells << " Dead Cells: " << dead_cells << std::endl;
    return 0;
}