#include "../include/functions.hpp"

Generation calculateNextGen(const Generation& current_gen)
{   
    /*
    int num_of_row_processes = num_of_processes; // num of processes in one row of the matrix
    int num_of_col_processes = num_of_processes; // num of processes in one column of the matrix
    */

    int N = current_gen.getGeneration().getSize();
    std::vector<std::vector<Cell>> next_gen_cells;

    for (int i = 0; i < N; ++i)
    {
        int upper_row_idx = (i - 1 + N) % N;
        int lower_row_idx = (i + 1) % N;
        std::vector<Cell> next_gen_rows;

        for (int j = 0; j < N; ++j)
        {
            const int left_col_idx = (j - 1 + N) % N;
            int right_col_idx = (j + 1) % N;
            int alive_neighbours_count = current_gen.countAliveNeighbours(left_col_idx, right_col_idx, lower_row_idx, upper_row_idx, j, i);
            bool isAlive = current_gen.getGeneration().getMatrix()[i][j].isAlive();

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
        next_gen_cells.push_back(std::move(next_gen_rows));
    }

    return Generation(Matrix{std::move(next_gen_cells)});
}

void countAliveAndDeadCells(const Generation& gen, int &alive_count, int &dead_count)
{
    for (const auto &row : gen.getGeneration().getMatrix())
    {
        for (const auto &cell : row)
        {
            if (cell.isAlive())
            {
                alive_count++;
            }
            else
            {
                dead_count++;
            }
        }
    }
}