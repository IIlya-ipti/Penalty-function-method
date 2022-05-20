#pragma once
#include "includes.h"
#include "GradMethods.h"

/*
	func for random number
*/
double fRand(double fMin, double fMax);

/*
	implementation of the penalty function method
*/
std::vector<double> Penaltymethod(Func_constraint* f);