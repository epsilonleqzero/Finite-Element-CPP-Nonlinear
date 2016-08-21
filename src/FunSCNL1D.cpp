/*
 * FunSCNL1D.cpp
 *
 * This function is a subclass of the class PdeFun.cpp
 *
 *  Created on: Aug 13, 2016
 *      Author: Ted Kwan
 */

#include "FunSCNL1D.h"
#include <vector>

using namespace std;
using namespace arma;

FunSCNL1D::FunSCNL1D() {
	// Using pi from armadillo.
	a=2.0*datum::pi;
	m=a;
	c=a;
}

FunSCNL1D::FunSCNL1D(vector<double> params) {
	// Using pi from armadillo.
	c=params[0];
	a=params[1];
	m=params[2];
}

/**
 * overloaded method from parent class.
 *
 * @param x - Input matrix.
 * @return vector of values calculated.
 */
double FunSCNL1D::evalF(double x){
	// Calculate vectors.
	return (a*x+m*sinh(x)+c);
}

double FunSCNL1D::evalDF(double x){
	// Calculate vectors.
	return (a+m*cosh(x));
}

FunSCNL1D::~FunSCNL1D() {
	// TODO Auto-generated destructor stub
}

