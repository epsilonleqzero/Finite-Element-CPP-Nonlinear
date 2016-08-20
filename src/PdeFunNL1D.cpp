/*
 * PdeFunNL1D.cpp
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

#include "PdeFunNL1D.h"


using namespace std;

/**
 * Constructor just in case the class gets called
 * directly.
 */
PdeFunNL1D::PdeFunNL1D() {
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
PdeFunNL1D::PdeFunNL1D(vector<double> params){
	this->params=params;
	dim=1;
}

/**
 * Evaluates the function value at the current values given.
 *
 * @param x - input matrix.
 * @return - vector containing the output.
 */
double PdeFunNL1D::evalF(double x){
	return x;
}

double PdeFunNL1D::evalDF(double x){
	return 1.0;
}

PdeFunNL1D::~PdeFunNL1D() {
	// TODO Auto-generated destructor stub
}
