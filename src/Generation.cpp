#include "../include/Generation.hpp"
#include <fstream>
#include <iostream>
#include <random>
#include <memory>

Generation::Generation(std::vector<std::vector<Cell>> _generation) : generation(_generation)
{
}

const std::vector<std::vector<Cell>> &Generation::getGeneration() const
{
    return generation;
}

std::vector<std::vector<Cell>> &Generation::getGeneration()
{
    return generation;
}

Generation::Generation(int row_size, int col_size, float prob_of_life)
{

    std::vector<std::vector<Cell>> _generation;
    std::default_random_engine random_engine;
    std::uniform_real_distribution<> uniform_zero_to_one(0.0, 1.0);

    for (int i = 0; i < row_size; i++)
    {
        std::vector<Cell> row;

        for (int j = 0; j < col_size; j++)
        {
            if (uniform_zero_to_one(random_engine) > prob_of_life)
            {
                row.push_back(Cell{'d'});
            }
            else
            {
                row.push_back(Cell{'a'});
            }
        }

        _generation.push_back(std::move(row));
    }

    this->generation = _generation;
}

int Generation::getRowSize() const
{
    return this->generation.size();
}

int Generation::getColSize() const
{
    return this->generation[0].size();
}

int Generation::countAliveNeighbours(const int &left_col_idx, const int &right_col_idx, const int &lower_row_idx, const int &upper_row_idx, const int &col_idx, const int &row_idx) const
{

    int alive_neighbours_count = 0;
    const std::vector<std::vector<Cell>> &current_gen_cells = generation;

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

void Generation::printGeneration(const std::string &filename) const
{
    std::ofstream file("./debug/" + filename + ".csv");

    if (!file.is_open())
    {
        std::cout << "Error opening file." << std::endl;
        return;
    }

    const auto &matrix_current_gen = generation;
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
