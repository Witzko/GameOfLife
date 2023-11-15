#include "../include/Cell.hpp"

Cell::Cell(char s) : state(s) {}

char Cell::getState() const
{
    return this->state;
}

void Cell::setStateToDead()
{
    this->state = 'd';
}

void Cell::setStateToAlive()
{
    this->state = 'a';
}

bool Cell::isAlive() const
{
    if (this->state == 'a'){
        return true;
    }
    
    return false;
}
