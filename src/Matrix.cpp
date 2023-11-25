#include "../include/Matrix.hpp"
#include <random>
#include <memory>

/**
 Default Constructor
 Constructs a matrix given a std::vector<std::vector<Cell>>

 @param _matrix The Matrix member
 @return Matrix Object
*/
Matrix::Matrix(std::vector<std::vector<Cell>> _matrix)
{
    this->matrix = _matrix;
}

/**
 Parameter Constructor
 Constructs a matrix of size NxN and the grid of cells with a certain probability of being dead or alive for each cell

 @param N Matrix dimension NxN
 @param prob_of_life Probability of individual cells being alive
 @return Matrix Object
*/
Matrix::Matrix(int N, float prob_of_life)
{

    std::vector<std::vector<Cell>> _matrix;
    std::default_random_engine random_engine;
    std::uniform_real_distribution<> uniform_zero_to_one(0.0, 1.0);

    for (int i = 0; i < N; i++)
    {
        std::vector<Cell> row;

        for (int j = 0; j < N; j++)
        {
            if (uniform_zero_to_one(random_engine) > prob_of_life)
            {
                row.push_back(Cell{'d'});
            }
            else
            {
                row.push_back(Cell{'a'});
            }
        }

        _matrix.push_back(std::move(row));
    }

    this->matrix = _matrix;
}

/**
  Returns the Size of the Matrix
  @return Size of Matrix
*/
int Matrix::getSize() const
{
    return this->matrix.size();
}

/**
  Getter Function of the matrix class
  @return Matrix member
*/
const std::vector<std::vector<Cell>> &Matrix::getMatrix() const
{
    return matrix;
}