/*
 * PdeFunNL1D.h
 *
 *  Created on: Aug 14, 2016
 *      Author: devils
 */

#include <armadillo>
#include <vector>

#ifndef PDEFUNNL1D_H_
#define PDEFUNNL1D_H_

class PdeFunNL1D {
public:
	PdeFunNL1D();
	PdeFunNL1D(std::vector<double> params);
	virtual double evalF(double x);
	virtual double evalDF(double x);
	virtual ~PdeFunNL1D();
	int dim;
private:
    std::vector<double> params;
};


#endif /* PDEFUN_H_ */
