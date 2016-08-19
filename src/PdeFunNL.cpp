/*
 * PdeFunNL.cpp
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

#include "PdeFunNL.h"

#include <cmath>

using namespace std;
using namespace arma;

/**
 * Constructor just in case the class gets called
 * directly.
 */
PdeFunNL::PdeFunNL() {
	vector<double> start(1);
	start[0]=1;
	params=start;
	dim =1;
}

/**
 * Constructor which sets a random vector to the parameters.
 * Not needed currently, but it may be needed later.
 *
 * @param params - vector containing parameters to be used.
 */
PdeFunNL::PdeFunNL(vector<double> params){
	this->params=params;
	dim=1;
}

/**
 * Evaluates the function value at the current values given.
 *
 * @param x - input matrix.
 * @return - vector containing the output.
 */
vec PdeFunNL::evalF(vec x){
	return x;
}

vec PdeFunNL::evalDF(vec x){
	return ones<vec>(x.n_rows);
}

PdeFunNL::~PdeFunNL() {
	// TODO Auto-generated destructor stub
}
