#pragma once

class Cell
{
    char state;
    
public:
    explicit Cell(char s);

    char getState() const;
    void setStateToDead();
    void setStateToAlive();
    bool isAlive() const;
};