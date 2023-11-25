#include "../include/functions.hpp"

/**
    Calculates the next Generation, based on the current Generation in a single process

    @param current_gen Generation reference object
    @return Generation Object
*/
Generation calculateNextGenSequentially(const Generation &current_gen)
{
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
            bool isAlive = current_gen.getGenerationCells()[i][j].isAlive();

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

/**
    Calculates the next Generation, based on the current Generation in multiple processes using MPI

    @param current_gen Generation reference object
    @param num_of_processes number of processes
    @return Generation Object
*/
Generation calculateNextGenParallel(const Generation &current_gen, int num_of_processes)
{

    int num_of_row_processes = num_of_processes; // num of processes in one row of the matrix
    int num_of_col_processes = num_of_processes; // num of processes in one column of the matrix

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
            bool isAlive = current_gen.getGenerationCells()[i][j].isAlive();

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

/**
    Iterates over the grid of a generation and increments the counters correspondingly

    @param gen Generation object reference
    @param alive_count alive count reference
    @param dead_count dead count reference
    @return void
*/
void countAliveAndDeadCells(const Generation &gen, int &alive_count, int &dead_count)
{
    for (const auto &row : gen.getGenerationCells())
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

/**
    Takes two generations and compares them entry by entry. generations must be of equal size

    @param gen_one Generation object reference
    @param gen_two Generation object reference
    @return boolean
*/
bool areGenerationsEqual(const Generation &gen_one, const Generation &gen_two)
{

    const std::vector<std::vector<Cell>> &gen_one_cells = gen_one.getGenerationCells();
    const std::vector<std::vector<Cell>> &gen_two_cells = gen_two.getGenerationCells();

    if (gen_one_cells.size() != gen_two_cells.size())
    {
        return false;
    }

    for (size_t i = 0; i < gen_one_cells.size(); i++)
    {
        if (gen_one_cells[i].size() != gen_two_cells[i].size())
        {
            return false;
        }

        if (gen_one_cells[i] != gen_two_cells[i]) // std::vector supports "==" operator, as well as our Cell class
        {
            return false;
        }
    }

    return true;
}

/**
    Retrieves subMatrix of specified size from an input matrix

    @param matrix Matrix object reference
    @param start_row Start row of main matrix
    @param start_col Start col of main matrix
    @param num_rows Number of rows for resulting submatrix
    @param num_cols Number of cols for resulting submatrix
    @return Matrix Object
*/
Matrix getSubMatrix(const Matrix &matrix, int start_row, int start_col, int num_rows, int num_cols)
{

    std::vector<std::vector<Cell>> sub_matrix_cells;

    for (int i = start_row; i < start_row + num_rows; i++)
    {
        std::vector<Cell>::const_iterator first_row_entry = matrix.getMatrix()[i].begin() + start_col;
        std::vector<Cell>::const_iterator last_row_entry = matrix.getMatrix()[i].begin() + (start_col + num_cols);
        sub_matrix_cells.push_back(std::vector<Cell>(first_row_entry, last_row_entry));
    }

    return Matrix{sub_matrix_cells};
}