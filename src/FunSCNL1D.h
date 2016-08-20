/*
 * FunSCNL1D.h

 *
 *  Created on: Aug 13, 2016
 *      Author: Devils
 */

#include <armadillo>
#include <vector>
#include "PdeFunNL1D.h"

#ifndef FUNSCNL1D_H_
#define FUNSCNL1D_H_

class FunSCNL1D: public PdeFunNL1D {
public:
	FunSCNL1D();
	FunSCNL1D(std::vector<double>);
	virtual ~FunSCNL1D();
	double evalF(double x);
	double evalDF(double x);
private:
	double a;
	double m;
	double c;
};

#endif /* FUNSC_H_ */
