#include "../include/functions.hpp"
#include <numeric>
#include <iostream>
#include <mpi.h>

Generation calculateNextGenSequentially(const Generation &current_gen)
{
    int row_size = current_gen.getRowSize();
    int col_size = current_gen.getColSize();
    std::vector<Cell> next_gen_cells(row_size * col_size, Cell{'d'});
    Generation next_gen(std::move(next_gen_cells), row_size, col_size);

    for (int i = 0; i < row_size; ++i)
    {
        int upper_row_idx = (i - 1 + row_size) % row_size;
        int lower_row_idx = (i + 1) % row_size;

        for (int j = 0; j < col_size; ++j)
        {
            int left_col_idx = (j - 1 + col_size) % col_size;
            int right_col_idx = (j + 1) % col_size;
            int alive_neighbours_count = current_gen.countAliveNeighbours(left_col_idx, right_col_idx, lower_row_idx, upper_row_idx, j, i);
            bool isAlive = current_gen.getCell(i, j).isAlive();

            if (!isAlive && alive_neighbours_count == 3)
            {
                next_gen.getCell(i, j).setState('a');
            }
            else if (isAlive && (alive_neighbours_count == 3 || alive_neighbours_count == 2))
            {
                next_gen.getCell(i, j).setState('a');
            }
        }
    }

    return next_gen;
}

void printGrid(std::vector<Cell> &vector, int rows, int columns)
{
    int i, j;
    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < columns; j++)
        {
            printf("%c ", vector[i * columns + j].getState());
        }
        printf("\n");
    }
    fflush(stdout);
}

Generation calculateNextGenParallelWCollNeighbourComm(Generation &&current_gen, MPI_Comm &cart_comm, MPI_Datatype &MPI_CELL, MPI_Datatype &MPI_COL_PADDING_WHALO, int halo_layer_size)
{
    int row_size = current_gen.getRowSize();
    int col_size = current_gen.getColSize();
    int row_size_whalo = row_size + 2 * halo_layer_size;
    int col_size_whalo = col_size + 2 * halo_layer_size;

    int dims[2];
    int rank, size;
    MPI_Comm_rank(cart_comm, &rank);
    MPI_Comm_size(cart_comm, &size);

    const int ndim = 2;
    int coords[ndim];
    int periods[ndim] = {1, 1};

    MPI_Cart_get(cart_comm, ndim, dims, periods, coords);

    /*
        Create new generation with halo layer and move the current generation to it:
        | halo   halo   halo |
        | halo curr_gen halo |
        | halo   halo   halo |
     */

    std::vector<Cell> current_gen_cells_whalo(row_size_whalo * col_size_whalo, Cell('d'));

    for (int i = 0; i < row_size; i++)
    {
        int padding_1 = i * col_size;
        int padding_2 = (i + 1) * col_size;
        int padding_3 = halo_layer_size + (i + 1) * col_size_whalo;
        std::copy(current_gen.getGeneration().begin() + padding_1, // Left col starting index
                  current_gen.getGeneration().begin() + padding_2, // Right col end index
                  current_gen_cells_whalo.begin() + padding_3);    // Starting index with halo
    }

    Generation current_gen_whalo(std::move(current_gen_cells_whalo), row_size_whalo, col_size_whalo);

    // Find the ranks of the neighbours
    int upper_left_rank, upper_rank, upper_right_rank, left_rank, right_rank, lower_left_rank, lower_rank, lower_right_rank;

    MPI_Cart_shift(cart_comm, 0, 1, &upper_rank, &lower_rank);
    MPI_Cart_shift(cart_comm, 1, 1, &left_rank, &right_rank);

    // Find ranks of the corners
    int upper_left_coords[2] = {coords[0] - 1, coords[1] - 1};
    int upper_right_coords[2] = {coords[0] - 1, coords[1] + 1};
    int lower_left_coords[2] = {coords[0] + 1, coords[1] - 1};
    int lower_right_coords[2] = {coords[0] + 1, coords[1] + 1};

    MPI_Cart_rank(cart_comm, upper_left_coords, &upper_left_rank);
    MPI_Cart_rank(cart_comm, upper_right_coords, &upper_right_rank);
    MPI_Cart_rank(cart_comm, lower_left_coords, &lower_left_rank);
    MPI_Cart_rank(cart_comm, lower_right_coords, &lower_right_rank);

    // Convert to Distributed Graph Comm
    MPI_Comm dist_graph_comm;

    int sources[] = {
        upper_rank,
        lower_rank,
        left_rank,
        right_rank,
        upper_left_rank,
        upper_right_rank,
        lower_left_rank,
        lower_right_rank,
    };

    // int source_weights[] = {0};
    // int rcv_weights[] = {0};

    MPI_Dist_graph_create_adjacent(cart_comm, 8, sources, MPI_UNWEIGHTED, 8,
                                   sources, MPI_UNWEIGHTED, MPI_INFO_NULL, false, &dist_graph_comm);

#ifdef DEBUG
    int test_sources[8];
    int test_destinations[8];
    int test_src_weights[8];
    int test_dest_weights[8];

    MPI_Dist_graph_neighbors(dist_graph_comm, 8, test_sources, test_src_weights, 8, test_destinations, test_dest_weights);

    MPI_Barrier(cart_comm);
    if (rank == 0)
    {
        std::cout << "col_size " << col_size << "\n";
        std::cout << "row_size " << row_size << "\n";
        std::cout << "col_size_whalo " << col_size_whalo << "\n";
        std::cout << "row_size_whalo " << row_size_whalo << "\n";
    }
    MPI_Barrier(cart_comm);
#endif

    int sendcounts[] = {col_size,
                        col_size,
                        1,
                        1,
                        1,
                        1,
                        1,
                        1};

    int sdispls[] = {
        static_cast<int>((&current_gen_whalo.getCell(1, 1) - &current_gen_whalo.getCell(0, 0))),
        static_cast<int>((&current_gen_whalo.getCell(row_size_whalo - 2, 1) - &current_gen_whalo.getCell(0, 0))),
        static_cast<int>((&current_gen_whalo.getCell(1, 1) - &current_gen_whalo.getCell(0, 0))),
        static_cast<int>((&current_gen_whalo.getCell(1, col_size_whalo - 2) - &current_gen_whalo.getCell(0, 0))),
        static_cast<int>((&current_gen_whalo.getCell(1, 1) - &current_gen_whalo.getCell(0, 0))),
        static_cast<int>((&current_gen_whalo.getCell(1, col_size_whalo - 2) - &current_gen_whalo.getCell(0, 0))),
        static_cast<int>((&current_gen_whalo.getCell(row_size_whalo - 2, 1) - &current_gen_whalo.getCell(0, 0))),
        static_cast<int>((&current_gen_whalo.getCell(row_size_whalo - 2, col_size_whalo - 2) - &current_gen_whalo.getCell(0, 0))),
    };

    int rdispls[] = {
        static_cast<int>((&current_gen_whalo.getCell(0, 1) - &current_gen_whalo.getCell(0, 0))),
        static_cast<int>((&current_gen_whalo.getCell(row_size_whalo - 1, 1) - &current_gen_whalo.getCell(0, 0))),
        static_cast<int>((&current_gen_whalo.getCell(1, 0) - &current_gen_whalo.getCell(0, 0))),
        static_cast<int>((&current_gen_whalo.getCell(1, col_size_whalo - 1) - &current_gen_whalo.getCell(0, 0))),
        static_cast<int>((&current_gen_whalo.getCell(0, 0) - &current_gen_whalo.getCell(0, 0))),
        static_cast<int>((&current_gen_whalo.getCell(0, col_size_whalo - 1) - &current_gen_whalo.getCell(0, 0))),
        static_cast<int>((&current_gen_whalo.getCell(row_size_whalo - 1, 0) - &current_gen_whalo.getCell(0, 0))),
        static_cast<int>((&current_gen_whalo.getCell(row_size_whalo - 1, col_size_whalo - 1) - &current_gen_whalo.getCell(0, 0))),
    };

    // int sdispls[] = {
    //     (1 * col_size_whalo + 1) * (int)sizeof(MPI_CELL),
    //     ((row_size_whalo - 2) * col_size_whalo + 1),
    //     (1 * col_size_whalo + 1),
    //     (2 * col_size_whalo - 1),
    //     (1 * col_size + 1),
    //     (1 * col_size + (col_size_whalo - 2)),
    //     ((row_size_whalo - 2) * col_size + 1),
    //     ((row_size_whalo - 2) * col_size + (col_size_whalo - 2)),
    // };

    // int rdispls[] = {
    //     (row_size_whalo * (col_size_whalo - 1) + 1) * (int)sizeof(MPI_CELL),
    //     (0 * col_size_whalo + 1),
    //     (2 * col_size_whalo + 0),
    //     (1 * col_size_whalo + 1),
    //     (0 * col_size + 0),
    //     (0 * col_size + (col_size_whalo - 1)),
    //     ((row_size_whalo - 1) * col_size + 0),
    //     ((row_size_whalo - 1) * col_size + (col_size_whalo - 1)),
    // };

    MPI_Datatype sendtypes[] = {
        MPI_CELL,
        MPI_CELL,
        MPI_COL_PADDING_WHALO,
        MPI_COL_PADDING_WHALO,
        MPI_CELL,
        MPI_CELL,
        MPI_CELL,
        MPI_CELL,
    };

#ifdef DEBUG
    MPI_Barrier(cart_comm);
    if (rank == 0)
    {
        printf("Rank 0 printed original\n");
        current_gen_whalo.printGeneration("original_0");
        for (size_t i = 0; i < 8; i++)
        {
            std::cout << sdispls[i] << "\n";
        }
    }
    else
    {
        printf("Rank 1 printed original\n");
        current_gen_whalo.printGeneration("original_1");
    }

    MPI_Barrier(cart_comm);
#endif

    MPI_Alltoallw(&current_gen_whalo.getCell(0, 0), sendcounts, sdispls, sendtypes, &current_gen_whalo.getCell(0, 0), sendcounts, rdispls, sendtypes, dist_graph_comm);

#ifdef DEBUG
    MPI_Barrier(cart_comm);
    if (rank == 0)
    {
        printf("Rank 0 printed after\n");
        current_gen_whalo.printGeneration("alltoall_0");
    }
    else
    {
        printf("Rank 1 printed after\n");
        current_gen_whalo.printGeneration("alltoall_1");
    }
    MPI_Barrier(cart_comm);
#endif

    // Reuse the Generation object current_gen which is already defined and store the next generation in it.
    std::vector<Cell> next_gen_cells(row_size * col_size, Cell('d'));
    current_gen.setGenerationAndProperties(std::move(next_gen_cells), row_size, col_size);

    // Iterate over the inner grid of the halo layer generation and count neighbours. Store the iterated cells into the next_gen
    for (int i = 1; i < row_size_whalo - 1; ++i)
    {
        int upper_row_idx = i - 1;
        int lower_row_idx = i + 1;

        for (int j = 1; j < col_size_whalo - 1; ++j)
        {
            int left_col_idx = j - 1;
            int right_col_idx = j + 1;
            int alive_neighbours_count = current_gen_whalo.countAliveNeighbours(left_col_idx, right_col_idx, lower_row_idx, upper_row_idx, j, i);
            bool isAlive = current_gen_whalo.getCell(i, j).isAlive();

            if (!isAlive && alive_neighbours_count == 3)
            {
                current_gen.getCell(i - halo_layer_size, j - halo_layer_size).setStateToAlive();
            }
            else if (isAlive && (alive_neighbours_count == 3 || alive_neighbours_count == 2))
            {
                current_gen.getCell(i - halo_layer_size, j - halo_layer_size).setStateToAlive();
            }
        }
    }

    return current_gen;
}

Generation calculateNextGenParallel(Generation &&current_gen, MPI_Comm &cart_comm,
                                    MPI_Datatype &MPI_CELL, MPI_Datatype &MPI_COL_PADDING_WHALO, int halo_layer_size)
{
    int row_size = current_gen.getRowSize();
    int col_size = current_gen.getColSize();
    int row_size_whalo = row_size + 2 * halo_layer_size;
    int col_size_whalo = col_size + 2 * halo_layer_size;

    int dims[2];
    int rank, size;
    MPI_Comm_rank(cart_comm, &rank);
    MPI_Comm_size(cart_comm, &size);

    const int ndim = 2;
    int coords[ndim];
    int periods[ndim] = {1, 1};

    MPI_Cart_get(cart_comm, ndim, dims, periods, coords);

    /*
        Create new generation with halo layer and move the current generation to it:
        | halo   halo   halo |
        | halo curr_gen halo |
        | halo   halo   halo |
     */

    std::vector<Cell> current_gen_cells_whalo(row_size_whalo * col_size_whalo, Cell('d'));

    for (int i = 0; i < row_size; i++)
    {
        int padding_1 = i * col_size;
        int padding_2 = (i + 1) * col_size;
        int padding_3 = halo_layer_size + (i + 1) * col_size_whalo;
        std::copy(current_gen.getGeneration().begin() + padding_1, // Left col starting index
                  current_gen.getGeneration().begin() + padding_2, // Right col end index
                  current_gen_cells_whalo.begin() + padding_3);    // Starting index with halo
    }

    Generation current_gen_whalo(std::move(current_gen_cells_whalo), row_size_whalo, col_size_whalo);

    // Find the ranks of the neighbours
    int upper_left_rank, upper_rank, upper_right_rank, left_rank, right_rank, lower_left_rank, lower_rank, lower_right_rank;

    MPI_Cart_shift(cart_comm, 0, 1, &upper_rank, &lower_rank);
    MPI_Cart_shift(cart_comm, 1, 1, &left_rank, &right_rank);

    // Find ranks of the corners
    int upper_left_coords[2] = {coords[0] - 1, coords[1] - 1};
    int upper_right_coords[2] = {coords[0] - 1, coords[1] + 1};
    int lower_left_coords[2] = {coords[0] + 1, coords[1] - 1};
    int lower_right_coords[2] = {coords[0] + 1, coords[1] + 1};

    MPI_Cart_rank(cart_comm, upper_left_coords, &upper_left_rank);
    MPI_Cart_rank(cart_comm, upper_right_coords, &upper_right_rank);
    MPI_Cart_rank(cart_comm, lower_left_coords, &lower_left_rank);
    MPI_Cart_rank(cart_comm, lower_right_coords, &lower_right_rank);

    MPI_Request send_request[8], recv_request[8];
    int tag[8] = {0, 1, 2, 3, 4, 5, 6, 7};

    // Send (only the "inner" matrix of current_gen_whalo)
    MPI_Isend(&current_gen_whalo.getCell(1, 1), col_size, MPI_CELL, upper_rank, tag[0], cart_comm, &send_request[0]);                        // Upper border
    MPI_Isend(&current_gen_whalo.getCell(row_size_whalo - 2, 1), col_size, MPI_CELL, lower_rank, tag[1], cart_comm, &send_request[1]);       // Lower border
    MPI_Isend(&current_gen_whalo.getCell(1, 1), 1, MPI_COL_PADDING_WHALO, left_rank, tag[2], cart_comm, &send_request[2]);                   // Left border
    MPI_Isend(&current_gen_whalo.getCell(1, col_size_whalo - 2), 1, MPI_COL_PADDING_WHALO, right_rank, tag[3], cart_comm, &send_request[3]); // Right border

    MPI_Isend(&current_gen_whalo.getCell(1, 1), 1, MPI_CELL, upper_left_rank, tag[4], cart_comm, &send_request[4]);                                    // Upper-left corner
    MPI_Isend(&current_gen_whalo.getCell(1, col_size_whalo - 2), 1, MPI_CELL, upper_right_rank, tag[5], cart_comm, &send_request[5]);                  // Upper-right corner
    MPI_Isend(&current_gen_whalo.getCell(row_size_whalo - 2, 1), 1, MPI_CELL, lower_left_rank, tag[6], cart_comm, &send_request[6]);                   // Lower-left corner
    MPI_Isend(&current_gen_whalo.getCell(row_size_whalo - 2, col_size_whalo - 2), 1, MPI_CELL, lower_right_rank, tag[7], cart_comm, &send_request[7]); // Lower-right corner

    // Receive (address of the halo layer indexes)
    MPI_Irecv(&current_gen_whalo.getCell(0, 1), col_size, MPI_CELL, upper_rank, tag[1], cart_comm, &recv_request[0]);                        // Upper border
    MPI_Irecv(&current_gen_whalo.getCell(row_size_whalo - 1, 1), col_size, MPI_CELL, lower_rank, tag[0], cart_comm, &recv_request[1]);       // Lower border
    MPI_Irecv(&current_gen_whalo.getCell(1, 0), 1, MPI_COL_PADDING_WHALO, left_rank, tag[3], cart_comm, &recv_request[2]);                   // Left border
    MPI_Irecv(&current_gen_whalo.getCell(1, col_size_whalo - 1), 1, MPI_COL_PADDING_WHALO, right_rank, tag[2], cart_comm, &recv_request[3]); // Right border

    MPI_Irecv(&current_gen_whalo.getCell(0, 0), 1, MPI_CELL, upper_left_rank, tag[7], cart_comm, &recv_request[4]);                                    // Upper-left corner
    MPI_Irecv(&current_gen_whalo.getCell(0, col_size_whalo - 1), 1, MPI_CELL, upper_right_rank, tag[6], cart_comm, &recv_request[5]);                  // Upper-right corner
    MPI_Irecv(&current_gen_whalo.getCell(row_size_whalo - 1, 0), 1, MPI_CELL, lower_left_rank, tag[5], cart_comm, &recv_request[6]);                   // Lower-left corner
    MPI_Irecv(&current_gen_whalo.getCell(row_size_whalo - 1, col_size_whalo - 1), 1, MPI_CELL, lower_right_rank, tag[4], cart_comm, &recv_request[7]); // Lower-right corner

    MPI_Waitall(8, recv_request, MPI_STATUS_IGNORE);

    // Reuse the Generation object current_gen which is already defined and store the next generation in it.
    std::vector<Cell> next_gen_cells(row_size * col_size, Cell('d'));
    current_gen.setGenerationAndProperties(std::move(next_gen_cells), row_size, col_size);

    // Iterate over the inner grid of the halo layer generation and count neighbours. Store the iterated cells into the next_gen
    for (int i = 1; i < row_size_whalo - 1; ++i)
    {
        int upper_row_idx = i - 1;
        int lower_row_idx = i + 1;

        for (int j = 1; j < col_size_whalo - 1; ++j)
        {
            int left_col_idx = j - 1;
            int right_col_idx = j + 1;
            int alive_neighbours_count = current_gen_whalo.countAliveNeighbours(left_col_idx, right_col_idx, lower_row_idx, upper_row_idx, j, i);
            bool isAlive = current_gen_whalo.getCell(i, j).isAlive();

            if (!isAlive && alive_neighbours_count == 3)
            {
                current_gen.getCell(i - halo_layer_size, j - halo_layer_size).setStateToAlive();
            }
            else if (isAlive && (alive_neighbours_count == 3 || alive_neighbours_count == 2))
            {
                current_gen.getCell(i - halo_layer_size, j - halo_layer_size).setStateToAlive();
            }
        }
    }

    return current_gen;
}

void countAliveAndDeadCells(const Generation &gen, int &alive_count, int &dead_count)
{
    for (int i = 0; i < gen.getRowSize(); i++)
    {
        for (int j = 0; j < gen.getColSize(); j++)
        {
            if (gen.getCell(i, j).isAlive())
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
    const std::vector<Cell> &gen_one_cells = gen_one.getGeneration();
    const std::vector<Cell> &gen_two_cells = gen_two.getGeneration();

    if (gen_one_cells.size() != gen_two_cells.size())
    {
        return false;
    }

    for (size_t i = 0; i < gen_one_cells.size(); i++)
    {
        if (gen_one_cells[i] == gen_two_cells[i])
            continue; // std::vector supports "==" operator, as well as our Cell class
        return false;
    }

    return true;
}

/**
    I'm to lazy to fix it because it's not really used atm. Plez fix if needed.
    // Adam
*/
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
