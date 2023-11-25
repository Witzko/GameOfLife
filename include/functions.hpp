#pragma once
#include "Generation.hpp"
  
Generation calculateNextGenSequentially(const Generation& current_gen);
Generation calculateNextGenParallel(const Generation& current_gen, int num_of_processes);
void countAliveAndDeadCells(const Generation& gen, int &alive_count, int &dead_count);
Matrix getSubMatrix(const Matrix& matrix, int start_row, int start_col, int num_rows, int num_cols);
bool areGenerationsEqual(const Generation& gen_one, const Generation& gen_two);