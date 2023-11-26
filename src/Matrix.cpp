#include "../include/Matrix.hpp"
#include <random>
#include <memory>

Matrix::Matrix() {};
Matrix::Matrix(std::vector<std::vector<Cell>> _matrix)
{
    this->matrix = _matrix;
}

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


int Matrix::getSize() const
{
    return this->matrix.size();
}

const std::vector<std::vector<Cell>> &Matrix::getMatrix() const
{
    return matrix;
}

std::vector<std::vector<Cell>> &Matrix::getMatrix()
{
    return matrix;
}