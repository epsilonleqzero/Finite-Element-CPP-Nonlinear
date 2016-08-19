/*
 * Mesh.h
 *
 *  Created on: Aug 13, 2016
 *      Author: Devils
 */

#include <armadillo>
#include <vector>

#ifndef MESHMG_H_
#define MESHMG_H_

class MeshMG {
public:
	MeshMG();
	MeshMG(std::vector<double> meshprops);
	virtual ~MeshMG();
	void uniformrefine();
	arma::sp_mat stiffness;
	arma::sp_mat mass;
	arma::vec area;
	arma::mat node;
	arma::umat elem;
	arma::uvec bdNode;
	arma::uvec isbdNode;
	arma::uvec freeNode;
	arma::uword N;
	arma::umat HB;

private:

	std::vector<arma::sp_mat> assembleMatrix();
	arma::vec accumArrayM(arma::uvec subs,arma::vec ar,arma::uword N);
	void makeMesh(arma::vec xr,arma::vec yr);
	void findBoundary();
	arma::uvec coarseNodeFineIdx;
	double h;
	double n;
	arma::uword NT;
	};

#endif /* MESH_H_ */
