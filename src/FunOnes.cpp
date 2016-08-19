/*
 * FunOnes.cpp
 *
 * This function is a subclass of the class PdeFun.cpp
 *
 *  Created on: Aug 13, 2016
 *      Author: Ted Kwan
 */

#include "FunOnes.h"
#include <cmath>

using namespace std;
using namespace arma;

FunOnes::FunOnes() {

}

/**
 * overloaded method from parent class.
 *
 * @param x - Input matrix.
 * @return vector of values calculated.
 */
vec FunOnes::evalF(mat x){
	vec y=ones<vec>(x.n_rows);
	return y;
}

FunOnes::~FunOnes() {

}

