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
	std::vector<arma::sp_mat> assembleMatrix();
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
	std::vector<arma::umat> HBs;
	std::vector<arma::uvec> freeNodes;
	std::vector<arma::uvec> coarse2Fine;
	std::vector<arma::uword> Ns;
	std::vector<arma::sp_mat> stiffs;
	std::vector<arma::vec> masses;
	std::vector<arma::umat> elems;
	//std::vector<arma::vec> areas;
	std::vector<arma::mat> nodes;
	std::vector<arma::sp_mat> Pro;
	int refines;

private:

	arma::vec accumArrayM(arma::uvec subs,arma::vec ar,arma::uword N);
	void makeMesh(arma::vec xr,arma::vec yr);
	void ProHB();
	void findBoundary();
	arma::uvec coarseNodeFineIdx;
	double h;
	double n;
	arma::uword NT;
	};

#endif /* MESH_H_ */
