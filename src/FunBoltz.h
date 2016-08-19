/*
 * FunBoltz.h

 *
 *  Created on: Aug 13, 2016
 *      Author: Devils
 */

#include <armadillo>
#include "PdeFun.h"

#ifndef FUNBOLTZ_H_
#define FUNBOLTZ_H_

class FunBoltz: public PdeFun {
public:
	FunBoltz();
	virtual ~FunBoltz();
	arma::vec evalF(arma::mat x);
private:
	double omega;
};

#endif /* FunBoltz_H_ */
