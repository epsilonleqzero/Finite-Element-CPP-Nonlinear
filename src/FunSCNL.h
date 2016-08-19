/*
 * FunSCNL.h

 *
 *  Created on: Aug 13, 2016
 *      Author: Devils
 */

#include <armadillo>
#include "PdeFunNL.h"

#ifndef FUNSCNL_H_
#define FUNSCNL_H_

class FunSCNL: public PdeFunNL {
public:
	FunSCNL();
	virtual ~FunSCNL();
	arma::vec evalF(arma::vec x);
	arma::vec evalDF(arma::vec x);
private:
	double omega;
};

#endif /* FUNSC_H_ */
