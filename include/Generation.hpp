#pragma once
#include <memory>
#include <string>
#include <vector>
#include "Cell.hpp"

class Generation
{
  std::vector<Cell> _generation;
  int _rows;
  int _cols;

public:
  /**
    Default Constructor

    @return Generation object
  */
  Generation();

  /**
    Constructs the Generation object from a Matrix object

    @param _generation Matrix object
    @return Generation object
  */
  explicit Generation(std::vector<Cell>&& generation, int rows, int cols);

  /**
    Constructs a Generation of size NxN and the grid of cells with a certain probability of being dead or alive for each cell

    @param rows row size of generation matrix
    @param cols col size of generation matrix
    @param prob_of_life Probability of individual cells being alive
    @return Generation Object
*/
  Generation(int rows, int cols, float prob_of_life);

  /**
    Set generation vector as well as rows and cols.

    @return void
  */
  void setGenerationAndProperties(std::vector<Cell>&& generation, int rows, int cols);

  /**
      Returns the Row Size of the Generation matrix
      @return Size of Row in Generation
  */
  int getRowSize() const;

  /**
      Returns the Col Size of the Generation matrix
      @return Size of Col in Generation
  */
  int getColSize() const;

  /**
    Getter of the underlying Matrix member of the Generation object

    @return Matrix member
  */
  const std::vector<Cell> &getGeneration() const;

  /**
    non-const Getter of the underlying Matrix member of the Generation object

    @return Matrix member
  */
  std::vector<Cell> &getGeneration();

  /**
    Return cell at specific index in the generation grid.

    @param i Row index
    @param j Column index

    @return Cell at index
  */
  const Cell &getCell(int i, int j) const;

    /**
    non-const Return cell at specific index in the generation grid.

    @param i Row index
    @param j Column index

    @return Cell at index
  */
  Cell &getCell(int i, int j);

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
  int countAliveNeighbours(const int &left_col_idx, const int &right_col_idx, const int &lower_row_idx,
                           const int &upper_row_idx, const int &row_idx, const int &column_idx) const;

  /**
     Prints the Generation Object

    @param filename takes a filename without filetype ending
    @return void
  */
  void printGeneration(const std::string &filepath);
};