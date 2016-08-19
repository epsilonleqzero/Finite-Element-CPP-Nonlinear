/*
 * FunSC.h

 *
 *  Created on: Aug 13, 2016
 *      Author: Devils
 */

#include <armadillo>
#include "PdeFun.h"

#ifndef FUNSC_H_
#define FUNSC_H_

class FunSC: public PdeFun {
public:
	FunSC();
	virtual ~FunSC();
	arma::vec evalF(arma::mat x);
private:
	double omega;
};

#endif /* FUNSC_H_ */
