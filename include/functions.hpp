#pragma once
#include "Generation.hpp"
#include <mpi.h>

/**
    Calculates the next Generation, based on the current Generation in a single process

    @param current_gen Generation reference object
    @return Generation Object
*/
Generation calculateNextGenSequentially(const Generation &current_gen);

/**
    Calculates the next Generation, based on the current Generation in multiple processes using MPI

    @param current_gen Generation reference object
    @param num_of_processes number of processes
    @return Generation Object
*/
// Generation calculateNextGenParallel(const Generation &current_gen, const MPI_Comm &cart_comm);

/**
    Iterates over the grid of a generation and increments the counters correspondingly

    @param gen Generation object reference
    @param alive_count alive count reference
    @param dead_count dead count reference
    @return void
*/
void countAliveAndDeadCells(const Generation &gen, int &alive_count, int &dead_count);

/**
    Takes two generations and compares them entry by entry. generations must be of equal size

    @param gen_one Generation object reference
    @param gen_two Generation object reference
    @return boolean
*/
bool areGenerationsEqual(const Generation &gen_one, const Generation &gen_two);

/**
    Retrieves subMatrix of specified size from an input matrix

    @param matrix Matrix object reference
    @param start_row Start row of main matrix
    @param start_col Start col of main matrix
    @param num_rows Number of rows for resulting submatrix
    @param num_cols Number of cols for resulting submatrix
    @return Matrix Object
*/
Matrix getSubMatrix(const Matrix &matrix, int start_row, int start_col, int num_rows, int num_cols);

/**
    Averages the values in a std::vector

    @param vector std::vector const reference
    @return average of all values
*/
double averageVectorElements(const std::vector<double> &_vector);