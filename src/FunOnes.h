/*
 * FunSC.h

 *
 *  Created on: Aug 13, 2016
 *      Author: Devils
 */

#include <armadillo>
#include "PdeFun.h"

#ifndef FUNONES_H_
#define FUNONES_H_

class FunOnes: public PdeFun {
public:
	FunOnes();
	virtual ~FunOnes();
	arma::vec evalF(arma::mat x);
};

#endif /* FUNONES_H_ */
