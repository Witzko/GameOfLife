CXX := mpicxx
CXX_FLAGS := -Werror -Wall -Wextra -O3 -march=native -std=c++20

.DEFAULT_GOAL := sequential

# List of object files
OBJECTS_SEQUENTIAL := ./build/main_sequential.o ./build/Cell.o ./build/Generation.o ./build/functions.o
OBJECTS_PARALLEL := ./build/main_parallel.o ./build/Cell.o ./build/Generation.o ./build/functions.o

sequential: $(OBJECTS_SEQUENTIAL)
	@$(CXX) $(CXX_FLAGS) $(OBJECTS_SEQUENTIAL) -o ./build/sequential

parallel: $(OBJECTS_PARALLEL)
	@$(CXX) $(CXX_FLAGS) $(OBJECTS_PARALLEL) -o ./build/parallel

debug: CXX_FLAGS += -DDEBUG -g
debug: sequential parallel

./build/main_sequential.o: ./src/sequential/main.cpp ./include/Cell.hpp ./include/Generation.hpp ./include/functions.hpp
	@$(CXX) $(CXX_FLAGS) -c ./src/sequential/main.cpp -o ./build/main_sequential.o

./build/main_parallel.o: ./src/parallel/main.cpp ./include/Cell.hpp ./include/Generation.hpp ./include/functions.hpp
	@$(CXX) $(CXX_FLAGS) -c ./src/parallel/main.cpp -o ./build/main_parallel.o

./build/functions.o: ./src/functions.cpp ./include/Cell.hpp ./include/Generation.hpp
	@$(CXX) $(CXX_FLAGS) -c ./src/functions.cpp -o ./build/functions.o

./build/Generation.o: ./src/Generation.cpp ./include/Generation.hpp ./include/Cell.hpp
	@$(CXX) $(CXX_FLAGS) -c ./src/Generation.cpp -o ./build/Generation.o

./build/Cell.o: ./src/Cell.cpp ./include/Cell.hpp
	@$(CXX) $(CXX_FLAGS) -c ./src/Cell.cpp -o ./build/Cell.o

.PHONY: clean
clean:
	rm -rf build/*
