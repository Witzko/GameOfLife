#include "../include/Generation.hpp"
#include <fstream>
#include <iostream>

Generation::Generation(std::unique_ptr<Matrix> _current_gen) : current_gen(std::move(_current_gen))
{
}

Matrix &Generation::getCurrentGen() const
{
    return *current_gen;
}

Matrix &Generation::getNextGen() const
{
    return *next_gen;
}

void Generation::calculateNextGen(/* int num_of_processes */)
{   
    /*
    int num_of_row_processes = num_of_processes; // num of processes in one row of the matrix
    int num_of_col_processes = num_of_processes; // num of processes in one column of the matrix
    */

    int N = current_gen->getSize();
    std::vector<std::vector<Cell>> next_gen_cells;
    

    for (int i = 0; i < N; ++i)
    {
        std::vector<Cell> next_gen_rows;
        int upper_row_idx = (i - 1 + N) % N;
        int lower_row_idx = (i + 1) % N;

        for (int j = 0; j < N; ++j)
        {
            int left_col_idx = (j - 1 + N) % N;
            int right_col_idx = (j + 1) % N;
            int alive_neighbours_count = getAliveNeighboursCount(left_col_idx, right_col_idx, lower_row_idx, upper_row_idx, j, i);
            bool isAlive = current_gen->getMatrix()[i][j].isAlive();

            if (!isAlive && alive_neighbours_count == 3)
            {
                next_gen_rows.push_back(Cell{'a'});
            }
            else if (isAlive && (alive_neighbours_count == 3 || alive_neighbours_count == 2))
            {
                next_gen_rows.push_back(Cell('a'));
            }
            else
            {
                next_gen_rows.push_back(Cell{'d'});
            }
        };
        next_gen_cells.push_back(next_gen_rows);
    }

    Matrix next_gen_matrix{std::move(next_gen_cells)};
    next_gen = std::make_unique<Matrix>(std::move(next_gen_matrix));
}

int Generation::getAliveNeighboursCount(int &left_col_idx, int &right_col_idx, int &lower_row_idx, int &upper_row_idx,  int&col_idx, int &row_idx)
{

    int alive_neighbours_count{0};
    const std::vector<std::vector<Cell>> &current_gen_cells = current_gen->getMatrix();

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

void Generation::printGenerations()
{
    std::ofstream file("./debug/generation_debug.csv");

    if (!file.is_open())
    {
        std::cout << "Error opening file." << std::endl;
        return;
    }

    const auto &matrix_current_gen = current_gen->getMatrix();
    file << "Current gen  \n";

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

    const auto &matrix_next_gen = next_gen->getMatrix();
    file << "Next gen  \n";

    for (const auto &row : matrix_next_gen)
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
