/*
 * ArmaFuns.h
 *
 *  Created on: Aug 19, 2016
 *      Author: Devils
 */
#include <armadillo>

#ifndef ARMAFUNS_H_
#define ARMAFUNS_H_

class ArmaFuns {
public:
	ArmaFuns();
	arma::mat unique(arma:: mat);
	arma::umat uniqueu(arma:: umat);
	arma::uvec uniqueidx;
	virtual ~ArmaFuns();
};

#endif /* ARMAFUNS_H_ */
