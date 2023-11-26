#pragma once
#include "./Matrix.hpp"
#include <memory>
#include <string>

class Generation
{

  Matrix generation;

public:
  /**
    Default Constructor

    @return Generation object
  */
  Generation(){};

  /**
    Constructs the Generation object from a Matrix object

    @param _generation Matrix object
    @return Generation object
  */
  explicit Generation(Matrix _generation);

  /**
     Gets the Matrix  of the Generation object

    @return Matrix
  */
  const Matrix &getGeneration() const;

  /**
     Non-const getter of the Matrix of the Generation object

    @return Matrix
  */
  Matrix &getGeneration();

  /**
    Getter of the underlying Matrix member of the Generation object

    @return Matrix member
  */
  const std::vector<std::vector<Cell>> &getGenerationCells() const;

  /**
    non-const Getter of the underlying Matrix member of the Generation object

    @return Matrix member
  */
  std::vector<std::vector<Cell>> &getGenerationCells();

  /**
    Counts Alive Neighbours of cell. Each cell has 8 neighbours in the Grid

    @param left_col_idx Left col index
    @param right_col_idx Right col index
    @param lower_row_idx Lower row index
    @param upper_row_idx Upper col index
    @param col_idx col index
    @param row_idx row index
    @return integer alive count
  */
  int countAliveNeighbours(const int &left_col_idx, const int &right_col_idx, const int &lower_row_idx, const int &upper_row_idx, const int &row_idx, const int &column_idx) const;

  /**
     Prints the Generation Object

    @param filename takes a filename without filetype ending
    @return void
  */
  void printGeneration(const std::string &filepath) const;
};