#include "../include/Cell.hpp"

Cell::Cell(char s) : _state(s) {}

void Cell::setState(char s)
{
    this->_state = s;
}

char Cell::getState() const
{
    return this->_state;
}

void Cell::setStateToDead()
{
    this->_state = 'd';
}

void Cell::setStateToAlive()
{
    this->_state = 'a';
}

bool Cell::isAlive() const
{
    if (this->_state == 'a')
    {
        return true;
    }

    return false;
}

bool Cell::operator==(const Cell &other) const
{
    return _state == other._state;
}