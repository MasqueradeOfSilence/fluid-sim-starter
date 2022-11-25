//============================================================================
// Name        : fluidsimstarter.cpp
// Author      : Sean Flynn
// Version     : 1.0
// Copyright   : Your copyright notice
// Description : 2D fluid simulator in C++, Ansi-style
//
// How to run me without using an IDE: g++ -o sim *.cpp (inside of the sim directory)
// ./sim
//============================================================================

#include <iostream>
using namespace std;

// vcpkg install eigen3:x64-windows
#include <eigen3/Eigen/Sparse>

#include "Simulator.h"
#include "MacGrid.h"
#include "GridCell.h"

const int NUM_FRAMES = 200;

int main()
{
	MacGrid *grid = new MacGrid(13, 26, 1);
	Simulator simulator(grid);
	simulator.run(NUM_FRAMES);
	return 0;
}
