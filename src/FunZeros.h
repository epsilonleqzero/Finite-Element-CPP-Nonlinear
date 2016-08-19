/*
 * FunSC.h

 *
 *  Created on: Aug 13, 2016
 *      Author: Devils
 */

#include <armadillo>
#include "PdeFun.h"

#ifndef FUNZEROS_H_
#define FUNZEROS_H_

class FunZeros: public PdeFun {
public:
	FunZeros();
	virtual ~FunZeros();
	arma::vec evalF(arma::mat x);
};

#endif /* FUNSC_H_ */
