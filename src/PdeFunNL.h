/*
 * PdeFunNL.h
 *
 *  Created on: Aug 14, 2016
 *      Author: devils
 */

#include <armadillo>
#include <vector>

#ifndef PDEFUNNL_H_
#define PDEFUNNL_H_

class PdeFunNL {
public:
	PdeFunNL();
	PdeFunNL(std::vector<double> params);
	virtual arma::vec evalF(arma::vec x);
	virtual arma::vec evalDF(arma::vec x);
	virtual ~PdeFunNL();
	int dim;
private:
    std::vector<double> params;
};


#endif /* PDEFUN_H_ */
