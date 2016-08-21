/*
 * FunSCNL.cpp
 *
 * This function is a subclass of the class PdeFun.cpp
 *
 *  Created on: Aug 13, 2016
 *      Author: Ted Kwan
 */

#include "FunSCNL.h"

using namespace std;
using namespace arma;

FunSCNL::FunSCNL() {
	// Using pi from armadillo.
	omega=2.0*datum::pi;
}

/**
 * overloaded method from parent class.
 *
 * @param x - Input matrix.
 * @return vector of values calculated.
 */
vec FunSCNL::evalF(vec x){
	// Calculate vectors.
	//vec y1=(sin(omega*x.col(1)));
	return (sinh(x));
}

vec FunSCNL::evalDF(vec x){
	// Calculate vectors.
	//vec y1=(sin(omega*x.col(1)));
	return (cosh(x));
}

FunSCNL::~FunSCNL() {
	// TODO Auto-generated destructor stub
}

