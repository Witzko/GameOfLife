# General Concept and Ideas for the Project

In the following, we discuss general ideas on how to implement the project, e.g. which datatypes we use and which classes we define.

## Task description

We are implementing a two dimensional stencil computation on small integers, and for that purpose we choose *Game of Life*.
In the *Game of Life* a so-called **cell** can be either dead or alive. The cells are then distributed in a N x N matrix, which is 
illustrated in the following picture:

![](./pics/generic_matrix_scheme.PNG)

In the first generation, each cell is either dead or alive with a certain probability (input parameter in our program). In the *Game of Life* 
cells either live or die depending on the state of their neighbours. The following principles must be implemented in our case:

    - each cell has 8 neighbours
    - alive cell survives to next generation if:
        - 2 or 3 neigbours are alive
    - dead cell gets reborn to next generation if:
        - exactly 3 neighbours are alive
    - if not reborn or survive:
        - dead for next generation

These are the core principles to consider. Another important aspect of the task is, that the matrix is folded, so if we access an entry at a 
out-of-bounce position in an N x N matix:

    Matrix[N+1, 0] = Matrix[0,0]
    Matrix[0, N+1] = Matrix[0,0]
    Matrix[N+1, N+1] = Matrix[0,0]
    Matrix[N+1, 1] = Matrix[0, 1]
    (...... and so on)

The neighbour situation corresponding to these constraints, is illustrated below:

![](./pics/generic_matrix_neighbours.PNG)


So the goal is to calculate the state of the new generation, meaning calculating each state of each cell, a lot of times. In the task they say we should aim for at least
1000 generations. After that we should report the final state.

### MPI and processes

So how does MPI come into play here? We use independent sub-sections of the matrix and assign it to processes, which are by themselves "packed" into a grid:

![](./pics/generic_matrix_processes.PNG)

The calculations for these processes are only partly independent, because on the edges we need communication. However, if the block size for one process becomes big, we
can expect a almost linear speed up according to Prof. Träff.


## Task implementation

### Design ideas

There are two main considerations apart from the algorithms:

    - data structure for matrix
    - data type for cell state

We note the task description: "The cell contents should be represented by a short C data type, either char or a bit in a word. Your implementation must be space efficient
in the sense of maintaining at most two generations, the current and the next, that is at most two n×n matrices."

**Data Structure for Matrix**:

Some possibilities inlcude:

    - vector
    - array
    - list
    (...)

Vector is very often the correct choice because its feature rich, dynamic and provides data locality. Lets look at the advantages of array and list over a vector, and 
see if these advantages are important to us:

    - Vector vs. List: A list might be the better choice if our vector grows dynamically and the copying of that memory is expensive (list elements always stay at the same
                        memory location, vector elements don´t). Also list is better if we insert very often in the middle of the data struct. 
                        --> both advantages not relevant in our case

    - Vector vs. Array: An array might be the better choice if we want to allocate on the stack. 
                        --> not the case, our matrix size might lead to stack-overflow

Therefore *std::vector* is a good choice.

**Data Type for Cell State**:

Some possibilities include:

    - char: 8 Bit
    - short int (unsigned): 16 Bit
    - a single bit in a machine word, e.g. uint_64: 1 Bit

The advantage of the first two is obviously simplicity. However, the advantage of the last one is space efficiency for sure. 

*Char* for the beginning is a good choice, but bitwise operations on a machine word seem like a interesting (but more complicated) alternative.

**Do we need classes?**

Probably not, but they are convenient for the Cell. Look at this example:

class Cell
{
private:
    char state;
    
public:
    Cell(char s) : state(s) {}
    
    void setStateToDead() {
        this->state = 'd';
    }
    void setStateToAlive() {
        this->state = 'a';
    }
};

The member functions like getter, setter make things a little easier and safer. The alternative is to fill up the std::vector with chars and assign it directly. 
Also okay I guess... Memory concerns: Does the cell need more memory? No! member functions are not stored in the object but the class definition. Both are 1 byte!