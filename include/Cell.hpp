#pragma once

class Cell
{
    char _state;

public:
    /**
        Constructs the Cell object

        @param s char: 'd' (dead) or 'a' (alive)
        @return Cell object
    */
    explicit Cell(char s);

    /**
        Gets the state member of Cell object

        @return Cell object member "state"
    */
    char getState() const;

    /**
        Sets the state member of Cell object

        @return void
    */
    void setState(char s);

    /**
        Sets the state member of Cell object to Dead

        @return void
    */
    void setStateToDead();

    /**
        Sets the state member of Cell object to Alive

        @return void
    */
    void setStateToAlive();

    /**
        Checks if the state member of Cell object is alive

        @return boolean
    */
    bool isAlive() const;

    /**
        Comparison oberator "==": equals. Returns true if two chars are equal

        @return boolean
    */
    bool operator==(const Cell &other) const;
};