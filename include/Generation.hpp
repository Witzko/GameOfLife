#pragma once
#include "./Matrix.hpp"
#include <memory>

class Generation
{

  std::unique_ptr<Matrix> current_gen;
  std::unique_ptr<Matrix> next_gen;

public:
  Generation(std::unique_ptr<Matrix> _current_gen);
  Matrix& getCurrentGen() const;
  Matrix& getNextGen() const;
  void calculateNextGen(/* int num_of_processes */);
  int getAliveNeighboursCount(int &left_col_idx, int &right_col_idx, int &lower_row_idx, int &upper_row_idx, int &row_idx, int &column_idx);
  void printGenerations();
};