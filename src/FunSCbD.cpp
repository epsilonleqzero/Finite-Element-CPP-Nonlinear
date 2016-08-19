/*
 * FunSCbD.cpp
 *
 * This function is a subclass of the class PdeFun.cpp
 *
 *  Created on: Aug 13, 2016
 *      Author: Ted Kwan
 */

#include "FunSCbD.h"
#include "Mesh.h"

using namespace std;
using namespace arma;

FunSCbD::FunSCbD() {
	// Using pi from armadillo.
	omega=2*datum::pi;
}

/**
 * overloaded method from parent class.
 *
 * @param x - Input matrix.
 * @return vector of values calculated.
 */
vec FunSCbD::evalF(mat x){
	vec y=cos(omega*x.col(1))%sin(omega*x.col(0));
	return y;
}

FunSCbD::~FunSCbD() {
	// TODO Auto-generated destructor stub
}

