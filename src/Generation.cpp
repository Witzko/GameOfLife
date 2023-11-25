#include "../include/Generation.hpp"
#include <fstream>
#include <iostream>

Generation::Generation(Matrix _generation) : generation(_generation)
{
}

const Matrix &Generation::getGeneration() const
{
    return generation;
}

const std::vector<std::vector<Cell>>& Generation::getGenerationCells() const{
    return generation.getMatrix();
}

int Generation::countAliveNeighbours(const int &left_col_idx, const int &right_col_idx, const int &lower_row_idx, const int &upper_row_idx,  const int&col_idx, const int &row_idx) const
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

void Generation::printGeneration(const std::string& filename) const
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
