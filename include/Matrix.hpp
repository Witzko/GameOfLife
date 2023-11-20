#pragma once
#include <vector>
#include "./Cell.hpp"

class Matrix
{
    std::vector<std::vector<Cell>> matrix;

public:
    explicit Matrix(std::vector<std::vector<Cell>> _matrix); // default constructor
    Matrix(int N, float prob_of_life);               // parameter constructor: constructs matrix of size N x N
    int getSize() const;
    const std::vector<std::vector<Cell>>& getMatrix() const;
};