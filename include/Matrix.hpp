#pragma once
#include <vector>
#include "./Cell.hpp"

class Matrix
{
    std::vector<std::vector<Cell>> matrix;

public:
    Matrix(std::vector<std::vector<Cell>> _matrix); // default constructor
    Matrix(int N, float prob_of_life);               // parameter constructor: constructs matrix of size N x N
    void countAliveAndDeadCells(int &alive_count, int &dead_count);
    int getSize() const;
    const std::vector<std::vector<Cell>>& getMatrix() const;
};