#include "../include/Generation.hpp"
#include <fstream>
#include <iostream>
#include <random>
#include <memory>

Generation::Generation(std::vector<Cell> generation, int rows, int cols) : _generation(generation), _rows(rows), _cols(cols)
{
}

const std::vector<Cell> &Generation::getGeneration() const
{
    return _generation;
}

std::vector<Cell> &Generation::getGeneration()
{
    return _generation;
}

Cell &Generation::getCell(int i, int j) {
    return _generation[i * _cols + j];
}

const Cell &Generation::getCell(int i, int j) const {
    return _generation[i * _cols + j];
}

Generation::Generation(int rows, int cols, float prob_of_life)
{
    this->_rows = rows;
    this->_cols = cols;

    std::vector<Cell> generation(rows*cols, Cell('d'));
    this->_generation = generation;

    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<> uniform_zero_to_one(0.0, 1.0);

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            if (uniform_zero_to_one(gen) < prob_of_life)
            {
                this->getCell(i,j).setStateToAlive();
            }
        }
    }
}

int Generation::getRowSize() const
{
    return _rows;
}

int Generation::getColSize() const
{
    return _cols;
}

int Generation::countAliveNeighbours(const int &left_col_idx, const int &right_col_idx, const int &lower_row_idx,
                                     const int &upper_row_idx, const int &col_idx, const int &row_idx) const
{
    int alive_neighbours_count = 0;

    // anit-clockwise from bottom left
    this->getCell(lower_row_idx,left_col_idx).isAlive() && alive_neighbours_count++;
    this->getCell(lower_row_idx,col_idx).isAlive() && alive_neighbours_count++;
    this->getCell(lower_row_idx,right_col_idx).isAlive() && alive_neighbours_count++;
    this->getCell(row_idx,right_col_idx).isAlive() && alive_neighbours_count++;
    this->getCell(upper_row_idx,right_col_idx).isAlive() && alive_neighbours_count++;
    this->getCell(upper_row_idx,col_idx).isAlive() && alive_neighbours_count++;
    this->getCell(upper_row_idx,left_col_idx).isAlive() && alive_neighbours_count++;
    this->getCell(row_idx,left_col_idx).isAlive() && alive_neighbours_count++;

    return alive_neighbours_count;
}

void Generation::printGeneration(const std::string &filename)
{
    std::ofstream file("./debug/" + filename + ".csv");

    if (!file.is_open())
    {
        std::cout << "Error opening file." << std::endl;
        return;
    }

    file << "Generation  \n";

    for (int i = 0; i < _rows; i++)
    {
        for (int j = 0; j < _cols; j++)
        {
            if (this->getCell(i,j).isAlive())
            {
                file << "1 ";
            }
            else
            {
                file << "0 ";
            }
        }
        file << "\n";
    }

    file.close();
}
