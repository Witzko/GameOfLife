#pragma once
#include <vector>
#include "./Cell.hpp"

class Matrix
{
    std::vector<std::vector<Cell>> matrix;

public:
    /**
        Default constructor, constructs an empty matrix
        @return Matrix object
    */
    Matrix();
    
    /**
        Constructs a matrix given a std::vector<std::vector<Cell>>

        @param _matrix The Matrix member
        @return Matrix Object
    */
    explicit Matrix(std::vector<std::vector<Cell>> _matrix);

    /**
        Constructs a matrix of size NxN and the grid of cells with a certain probability of being dead or alive for each cell

        @param N Matrix dimension NxN
        @param prob_of_life Probability of individual cells being alive
        @return Matrix Object
    */
    Matrix(int N, float prob_of_life);

    /**
        Returns the Size of the Matrix
        @return Size of Matrix
    */
    int getSize() const;

    /**
        Getter Function of the matrix class
        @return Matrix member
    */
    const std::vector<std::vector<Cell>> &getMatrix() const;
};