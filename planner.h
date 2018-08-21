#ifndef myheader
#define myheader
#include <iostream>

//matrices
// #include "/usr/include/armadillo"
#include <vector>
#include <fstream>
#include <ctime>
#include <chrono>
#include <iomanip>
#include <string>
#include <utility>


struct state
{
	double* q;
	state *PREV;
	std::vector<state *> NEXT;
};
#endif