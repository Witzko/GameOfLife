#include "../include/Generation.hpp"
#include <fstream>
#include <iostream>

/**
 Default Constructor
 Constructs the Generation object from a Matrix object
 @param _generation Matrix object
 @return Generation object
*/
Generation::Generation(Matrix _generation) : generation(_generation)
{
}

/**
 Gets the Matrix  of the Generation object

 @return Matrix
*/
const Matrix &Generation::getGeneration() const
{
    return generation;
}

/**
 Gets the underlying Matrix member of the Generation object

 @return Matrix member
*/
const std::vector<std::vector<Cell>> &Generation::getGenerationCells() const
{
    return generation.getMatrix();
}

/**
 Counts Alive Neighbours of cell. Each cell has 8 neighbours in the Grid

 @param left_col_idx Left col index
 @param right_col_idx Right col index
 @param lower_row_idx Lower row index
 @param upper_row_idx Upper col index
 @param col_idx col index
 @param row_idx row index
 @return integer alive count
*/
int Generation::countAliveNeighbours(const int &left_col_idx, const int &right_col_idx, const int &lower_row_idx, const int &upper_row_idx, const int &col_idx, const int &row_idx) const
{

    int alive_neighbours_count{0};
    const std::vector<std::vector<Cell>> &current_gen_cells = generation.getMatrix();

    // clockwise from bottom left
    current_gen_cells[lower_row_idx][left_col_idx].isAlive() && alive_neighbours_count++;
    current_gen_cells[lower_row_idx][col_idx].isAlive() && alive_neighbours_count++;
    current_gen_cells[lower_row_idx][right_col_idx].isAlive() && alive_neighbours_count++;
    current_gen_cells[row_idx][right_col_idx].isAlive() && alive_neighbours_count++;
    current_gen_cells[upper_row_idx][right_col_idx].isAlive() && alive_neighbours_count++;
    current_gen_cells[upper_row_idx][col_idx].isAlive() && alive_neighbours_count++;
    current_gen_cells[upper_row_idx][left_col_idx].isAlive() && alive_neighbours_count++;
    current_gen_cells[row_idx][left_col_idx].isAlive() && alive_neighbours_count++;

    return alive_neighbours_count;
}

/**
 Prints the Generation Object

 @param filename takes a filename without filetype ending
 @return void
*/
void Generation::printGeneration(const std::string &filename) const
{
    std::ofstream file("./debug/" + filename + ".csv");

    if (!file.is_open())
    {
        std::cout << "Error opening file." << std::endl;
        return;
    }

    const auto &matrix_current_gen = generation.getMatrix();
    file << "Generation  \n";

    for (const auto &row : matrix_current_gen)
    {
        for (const auto &cell : row)
        {
            if (cell.isAlive())
            {
                file << "1 ";
            }
            else
            {
                file << "0 ";
            }
        }
        file << "\n";
    }

    file.close();
}
