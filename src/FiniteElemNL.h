/*
 * FiniteElemNL.h
 *
 *  Created on: Aug 13, 2016
 *      Author: Devils
 */
#include <armadillo>
#include <string>
#include <vector>
#include "MeshMG.h"
#include "PdeFun.h"
#include "PdeFunNL1D.h"
#include "PdeFunNL.h"
#include "FunSCNL1D.h"

#ifndef FINITEELEMNL_H_
#define FINITEELEMNL_H_

class FiniteElemNL {
public:
	FiniteElemNL();
	FiniteElemNL(std::vector<double> meshprops,std::string fun);
	virtual ~FiniteElemNL();
	MeshMG mesh;
	arma::vec u;
private:
	arma::vec calcRHS();
	arma::vec NWTsolve(arma::vec b,arma::vec uf, arma::uword maxitr,
			arma::vec M, arma::mat A,double tol,arma::uvec freeNode);
	arma::vec GSsolve(arma::vec b,arma::vec u, arma::uword maxitr,
				arma::vec M, arma::mat A,double tol,arma::uvec freeNode);
	arma::vec GSsolveb(arma::vec b,arma::vec u, arma::uword maxitr,
					arma::vec M, arma::mat A,double tol,arma::uvec freeNode);
	double NWTsolve1D(FunSCNL1D f, double x0, double tol, arma::uword maxitr);
	arma::vec accumArray(arma::uvec subs,arma::vec ar,arma::uword N);
	PdeFun * pde;
	PdeFunNL * pdenl;
	PdeFun * bdfun;
};

#endif /* FINITEELEM_H_ */
