#include "../include/functions.hpp"
#include <numeric>
#include <iostream>
#include <mpi.h>

Generation calculateNextGenSequentially(const Generation &current_gen)
{
    int row_size = current_gen.getRowSize();
    int col_size = current_gen.getColSize();
    std::vector<std::vector<Cell>> next_gen_cells;

    for (int i = 0; i < row_size; ++i)
    {
        int upper_row_idx = (i - 1 + row_size) % row_size;
        int lower_row_idx = (i + 1) % row_size;
        std::vector<Cell> next_gen_rows;

        for (int j = 0; j < col_size; ++j)
        {
            const int left_col_idx = (j - 1 + col_size) % col_size;
            int right_col_idx = (j + 1) % col_size;
            int alive_neighbours_count = current_gen.countAliveNeighbours(left_col_idx, right_col_idx, lower_row_idx, upper_row_idx, j, i);
            bool isAlive = current_gen.getGeneration()[i][j].isAlive();

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

    return Generation(std::move(next_gen_cells));
}

void calculateNextGenParallel(const Generation &current_gen, Generation &next_gen, MPI_Comm &cart_comm, bool weak_scaling_flag)
{

    int row_size = current_gen.getRowSize();
    int col_size = current_gen.getColSize();

    int dims[2];
    int rank, size;
    MPI_Comm_rank(cart_comm, &rank);
    MPI_Comm_size(cart_comm, &size);

    const int ndim = 2;
    int coords[ndim];
    int periods[ndim] = {1, 1};

    MPI_Cart_get(cart_comm, ndim, dims, periods, coords);

    int num_local_rows, num_local_cols = 0;

    if (weak_scaling_flag)
    {
        num_local_rows = row_size / size;
        num_local_cols = col_size / size;
    }
    else
    {
        num_local_rows = row_size / dims[0];
        num_local_cols = col_size / dims[1];
    }

    int start_row = coords[0] * num_local_rows;
    int start_col = coords[1] * num_local_cols;

    int end_row = (coords[0] + 1) * num_local_rows - 1;
    int end_col = (coords[1] + 1) * num_local_cols - 1;

#ifdef DEBUG
    MPI_Barrier(cart_comm);
    std::cout << "Rank " << rank << " has coordinates (" << coords[0] << ", " << coords[1] << ")" << std::endl;
    std::cout << "Start row: " << start_row << ", End row: " << end_row << ", Start col: " << start_col << ", End col: " << end_col << std::endl;
    MPI_Barrier(cart_comm);
#endif

    /*  !!! Delete this dummy value. Only there cause next_gen can not stay unused when compiling !!!*/
    char dummy = next_gen.getGeneration()[0][0].getState();
    std::cout << dummy << std::endl;

    /*
            ### Exercise 3 + 4  ###

            My boyz take over and can start implementing the logic here

            Note that you can do similar stuff like in the sequential version:

            just set the state inside the for loops with next_gen.getGeneration()[i][j].setStateToAlive();

    */
}

void countAliveAndDeadCells(const Generation &gen, int &alive_count, int &dead_count)
{
    for (const auto &row : gen.getGeneration())
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

bool areGenerationsEqual(const Generation &gen_one, const Generation &gen_two)
{

    const std::vector<std::vector<Cell>> &gen_one_cells = gen_one.getGeneration();
    const std::vector<std::vector<Cell>> &gen_two_cells = gen_two.getGeneration();

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

std::vector<std::vector<Cell>> getSubMatrix(const std::vector<std::vector<Cell>> &matrix, int start_row, int start_col, int num_rows, int num_cols)
{

    std::vector<std::vector<Cell>> sub_matrix_cells;

    for (int i = start_row; i < start_row + num_rows; i++)
    {
        std::vector<Cell>::const_iterator first_row_entry = matrix[i].begin() + start_col;
        std::vector<Cell>::const_iterator last_row_entry = matrix[i].begin() + (start_col + num_cols);
        sub_matrix_cells.push_back(std::vector<Cell>(first_row_entry, last_row_entry));
    }

    return sub_matrix_cells;
}

double averageVectorElements(const std::vector<double> &_vector)
{
    auto const _size = static_cast<double>(_vector.size());
    return std::reduce(_vector.begin(), _vector.end()) / _size;
}