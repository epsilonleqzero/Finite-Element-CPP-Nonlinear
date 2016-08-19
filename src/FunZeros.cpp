/*
 * FunZeros.cpp
 *
 * This function is a subclass of the class PdeFun.cpp
 *
 *  Created on: Aug 13, 2016
 *      Author: Ted Kwan
 */

#include "FunZeros.h"
#include <cmath>

using namespace std;
using namespace arma;

FunZeros::FunZeros() {

}

/**
 * overloaded method from parent class.
 *
 * @param x - Input matrix.
 * @return vector of values calculated.
 */
vec FunZeros::evalF(mat x){
	vec y=zeros<vec>(x.n_rows);
	return y;
}

FunZeros::~FunZeros() {
	// TODO Auto-generated destructor stub
}

