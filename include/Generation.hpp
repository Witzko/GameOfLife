#pragma once
#include "./Matrix.hpp"
#include <memory>
#include <string>

class Generation
{

  Matrix generation;

public:

  explicit Generation(Matrix _generation);
  const Matrix& getGeneration() const;
  const std::vector<std::vector<Cell>>& getGenerationCells() const;
  int countAliveNeighbours(const int &left_col_idx, const int &right_col_idx, const int &lower_row_idx, const int &upper_row_idx, const int &row_idx, const int &column_idx) const;
  void printGeneration(const std::string& filepath) const;
};