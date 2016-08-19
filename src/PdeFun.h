/*
 * PdeFun.h
 *
 *  Created on: Aug 14, 2016
 *      Author: devils
 */

#include <armadillo>
#include <vector>

#ifndef PDEFUN_H_
#define PDEFUN_H_

class PdeFun {
public:
	PdeFun();
	PdeFun(std::vector<double> params);
	virtual arma::vec evalF(arma::mat x);
	virtual ~PdeFun();
	int dim;
private:
    std::vector<double> params;
};


#endif /* PDEFUN_H_ */
