CXX := g++
CXX_FLAGS := -Werror -Wall -Wextra -g

.DEFAULT_GOAL := sequential

# List of object files
OBJECTS := ./build/main.o ./build/Matrix.o ./build/Cell.o ./build/Generation.o

sequential: $(OBJECTS)
	@$(CXX) $(CXX_FLAGS) $(OBJECTS) -o ./build/sequential

./build/main.o: ./src/main.cpp ./include/Matrix.hpp ./include/Cell.hpp ./include/Generation.hpp
	@$(CXX) $(CXX_FLAGS) -c ./src/main.cpp -o ./build/main.o

./build/Generation.o: ./src/Generation.cpp ./include/Generation.hpp ./include/Matrix.hpp ./include/Cell.hpp
	@$(CXX) $(CXX_FLAGS) -c ./src/Generation.cpp -o ./build/Generation.o

./build/Matrix.o: ./src/Matrix.cpp ./include/Matrix.hpp ./include/Cell.hpp
	@$(CXX) $(CXX_FLAGS) -c ./src/Matrix.cpp -o ./build/Matrix.o

./build/Cell.o: ./src/Cell.cpp ./include/Cell.hpp
	@$(CXX) $(CXX_FLAGS) -c ./src/Cell.cpp -o ./build/Cell.o