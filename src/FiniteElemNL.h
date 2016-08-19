/*
 * FiniteElemNL.h
 *
 *  Created on: Aug 13, 2016
 *      Author: Devils
 */
#include <armadillo>
#include <string>
#include <vector>
#include "Mesh.h"
#include "PdeFun.h"
#include "PdeFunNL.h"

#ifndef FINITEELEMNL_H_
#define FINITEELEMNL_H_

class FiniteElemNL {
public:
	FiniteElemNL();
	FiniteElemNL(std::vector<double> meshprops,std::string fun);
	virtual ~FiniteElemNL();
	Mesh mesh;
	arma::vec u;
private:
	arma::vec calcRHS();
	arma::vec NWTsolve(arma::vec b,arma::vec uf, arma::uword maxitr,
			arma::vec M, arma::mat A,double tol,arma::uvec freeNode);
	arma::vec accumArray(arma::uvec subs,arma::vec ar,arma::uword N);
	PdeFun * pde;
	PdeFunNL * pdenl;
	PdeFun * bdfun;
};

#endif /* FINITEELEM_H_ */
