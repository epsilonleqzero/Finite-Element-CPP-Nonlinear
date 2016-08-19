/*
 * PdeFun.cpp
 *
 * This is an abstract parent class to be used by the finite element
 * method (FiniteElem.cpp) so that different functions can be called.
 *
 * There is a single method, which is abstract and overloaded by
 * the classes which extend this class.
 *
 *  Created on: Aug 6, 2016
 *      Author: Ted Kwan
 */

#include "PdeFun.h"

#include <cmath>

using namespace std;
using namespace arma;

/**
 * Constructor just in case the class gets called
 * directly.
 */
PdeFun::PdeFun() {
	vector<double> start(1);
	start[0]=1;
	params=start;
	dim =2;
}

/**
 * Constructor which sets a random vector to the parameters.
 * Not needed currently, but it may be needed later.
 *
 * @param params - vector containing parameters to be used.
 */
PdeFun::PdeFun(vector<double> params){
	this->params=params;
	dim=2;
}

/**
 * Evaluates the function value at the current values given.
 *
 * @param x - input matrix.
 * @return - vector containing the output.
 */
vec PdeFun::evalF(mat x){
	return x;
}

PdeFun::~PdeFun() {
	// TODO Auto-generated destructor stub
}
