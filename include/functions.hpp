#pragma once
#include "Generation.hpp"
  
Generation calculateNextGen(const Generation& current_gen);

void countAliveAndDeadCells(const Generation& gen, int &alive_count, int &dead_count);