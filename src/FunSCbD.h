/*
 * FunSCbD.h

 *
 *  Created on: Aug 13, 2016
 *      Author: Devils
 */

#include <armadillo>
#include "PdeFun.h"

#ifndef FUNSCBD_H_
#define FUNSCBD_H_

class FunSCbD: public PdeFun {
public:
	FunSCbD();
	virtual ~FunSCbD();
	arma::vec evalF(arma::mat x);
private:
	double omega;
};

#endif /* FunSCbD_H_ */
