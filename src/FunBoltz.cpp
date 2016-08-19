/*
 * FunBoltz.cpp
 *
 * This function is a subclass of the class PdeFun.cpp
 *
 *  Created on: Aug 13, 2016
 *      Author: Ted Kwan
 */

#include "FunBoltz.h"
#include "Mesh.h"

using namespace std;
using namespace arma;

FunBoltz::FunBoltz() {
	// Using pi from armadillo.
	omega=2*datum::pi;
}

/**
 * overloaded method from parent class.
 *
 * @param x - Input matrix.
 * @return vector of values calculated.
 */
vec FunBoltz::evalF(mat x){
	vec ub=(0.1*ones<vec>(x.n_rows)
			+((x.col(0)+2*x.col(1))/sqrt(5.0)));
	vec y=log((ones<vec>(x.n_rows)+cos(ub))/((ones<vec>(x.n_rows)-cos(ub))));
	return y;
}

FunBoltz::~FunBoltz() {
	// TODO Auto-generated destructor stub
}

