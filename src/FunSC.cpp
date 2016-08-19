/*
 * FunSC.cpp
 *
 * This function is a subclass of the class PdeFun.cpp
 *
 *  Created on: Aug 13, 2016
 *      Author: Ted Kwan
 */

#include "FunSC.h"
#include "Mesh.h"

using namespace std;
using namespace arma;

FunSC::FunSC() {
	// Using pi from armadillo.
	omega=2.0*datum::pi;
}

/**
 * overloaded method from parent class.
 *
 * @param x - Input matrix.
 * @return vector of values calculated.
 */
vec FunSC::evalF(mat x){
	// Calculate vectors.
	vec y1=(cos(omega*x.col(1)));
	vec y2=sin(omega*x.col(0));
	// Schur product.
	vec y=y1%y2;
	// Scale return value.
	return 2*pow(omega,2.0)*y;
}

FunSC::~FunSC() {
	// TODO Auto-generated destructor stub
}

