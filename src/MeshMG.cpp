/*
 * Mesh.cpp
 *
 * This class contains the mesh and mesh properties to be used
 * in the nonlinear finite element method (FiniteElemNL.cpp).
 *
 * It constructs the mesh, then constructs the stiffness matrix and
 * calculates the boundary. It can be refined by calling
 * uniformrefine().
 *
 *  Created on: Aug 13, 2016
 *      Author: Ted Kwan
 */

#include "MeshMG.h"
#include "ArmaFuns.h"
#include <cmath>
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

MeshMG::MeshMG() {

}

/**
 * Constructor for MeshMG. It takes one input, a vector containing the
 * properties of the mesh to be created.
 *
 * @param meshprops - double vector with all of the necessary properties.
 */
MeshMG::MeshMG(vector<double> meshprops) {
	// Extract the mesh size value.
	h = meshprops[4];
	n = round(1 / h);
	// Create vectors to be used for mesh creation.
	vec xr = linspace<vec>(meshprops[0], meshprops[1], n + 1);
	vec yr = linspace<vec>(meshprops[2], meshprops[3], n + 1);
	// Create the mesh.
	makeMesh(xr, yr);
	// Assemble stiffness matrix
	vector<sp_mat> stiffmass = assembleMatrix();
	stiffness=stiffmass[0];
	mass=stiffmass[1];
	// Calculate boundary nodes.
	findBoundary();
}

/**
 * makeMeshMG creates the mesh for the 2D Cartesian grid.
 *
 * @param xr - input vector containing the x values.
 * @param yr - input vector containing the y values.
 */
void MeshMG::makeMesh(vec xr, vec yr) {
	// Findthe lengths to use.
	uword xrl = xr.n_rows;
	uword yrl = xr.n_rows;
	mat x = zeros<mat>(xrl, yrl);
	mat y = zeros<mat>(xrl, yrl);
	// Construct the meshes. This is just the same
	// as calling meshgrid() in MATLAB.
	for (uword i = 0; i < yrl; i++) {
		x.row(i) = xr.t();
	}
	for (uword i = 0; i < xrl; i++) {
		y.col(i) = yr;
	}
	// Vectorise the values to use for the nodes matrix.
	vec xv = vectorise(x);
	N = xv.n_rows;
	Ns.push_back(N);
	uword ni = x.n_rows;
	// Initialize node matrix.
	node = zeros<mat>(N, 2);
	// put all the values in the node matrix
	node.col(0) = xv;
	node.col(1) = vectorise(y);
	// Map the indices for the nodes to create elem matrix.
	uvec t2nidxMapnz = regspace<uvec>(0, N - ni - 1);
	uvec topNode = regspace<uvec>(ni - 1, ni, N - ni - 1);
	// Set all of the doubled nodes to not be included.
	t2nidxMapnz(topNode) = zeros<uvec>(topNode.n_rows);
	// Find nonzeros in the original indices map.
	uvec nnz = nonzeros(t2nidxMapnz);
	uvec k = zeros<uvec>(nnz.n_rows + 1);
	// Calculate the different values used to create nodes.
	k(span(1, nnz.n_rows)) = nnz;
	uword NE = k.n_rows;
	HBs.push_back(zeros<umat>(3,3));
	coarse2Fine.push_back(zeros<uvec>(3));
	Pro.push_back(zeros<sp_mat>(3,3));
	// Create the elements which will differ by odd and
	// even elements.
	umat elemup = zeros<umat>(NE, 3);
	umat elemdown = zeros<umat>(NE, 3);
	uvec niv = ones<uvec>(k.n_rows);
	niv = niv * ni;
	uvec onek = ones<uvec>(k.n_rows);
	// Map elements to nodes in order.
	elemup.col(0) = k + niv;
	elemup.col(1) = k + niv + onek;
	elemup.col(2) = k;
	elemdown.col(0) = k + onek;
	elemdown.col(1) = k;
	elemdown.col(2) = k + niv + onek;
	// join all columns together and stack them properly.
	elem = join_cols(elemup, elemdown);
	NT = elem.n_rows;
	refines=0;
}

/**
 * assembleMatrix assembles the stiffness matrix using the quick
 * construction method from armadillo for sparse matrices.
 *
 * @return - Stiffness matrix as a sparse matrix.
 */
vector<sp_mat> MeshMG::assembleMatrix() {
	// Initialize index maps and value vector.
	uvec ii = zeros<uvec>(9 * NT);
	uvec jj = zeros<uvec>(9 * NT);
	vec sA = zeros<vec>(9 * NT);
	// Calculate the area of each node.
	cube ve(NT, 2, 3);
	ve.slice(0) = node.rows(elem.col(2)) - node.rows(elem.col(1));
	ve.slice(1) = node.rows(elem.col(0)) - node.rows(elem.col(2));
	ve.slice(2) = node.rows(elem.col(1)) - node.rows(elem.col(0));
	// Find the area using the dot product on the second dimension.
	area = 0.5*abs((ve.slice(2).col(0) % ve.slice(1).col(1))
					-(ve.slice(2).col(1) % ve.slice(1).col(0)));
	uword index = 0;
	// Loop to map values and indices.
	for (uword i = 0; i < 3; i++) {
		for (uword j = 0; j < 3; j++) {
			// Setup indices to map in this iteration.
			uvec inds = regspace<uvec>(index, index + NT - 1);
			// Setup element maps for indices.
			ii(inds) = elem.col(i);
			jj(inds) = elem.col(j);
			// Calculate values of stiffness matrix at
			// these points.
			mat prod = ve.slice(i) % ve.slice(j);
			// Store calculated value.
			sA(inds) = sum(prod, 1) / (4 * area);
			index = index + NT;
		}
	}
	// Setup index map to be a 2x9NT matrix.
	umat inds = join_horiz(ii, jj);
	// Create the sparse matrix using the same
	// method as sparse(row indices, col indices, values, size)
	sp_mat A(true, inds.t(), sA, N, N, true, true);
	vector<sp_mat> matvec(2);
	matvec[0]=A;
	// Create mass matrix.
	vec Mv=accumArrayM(join_vert(elem.col(0), join_vert(elem.col(1),elem.col(2)))
			,join_vert(area, join_vert(area,area))/(3.0), N);
	uvec inds2=regspace<uvec>(0,Mv.n_rows-1);
	umat subs2=join_horiz(inds2,inds2);
	sp_mat M(true,subs2.t(),Mv,N,N,true,true);
	matvec[1]=M;
	masses.push_back(Mv);
	stiffs.push_back(A);
	return matvec;
}

/**
 * Find the boundary nodes and elements.
 *
 * This method calculates the boundary for the given mesh to be used
 * to set the boundary condition for the finite element method.
 *
 * Values are stored so that they can be accesssed later.
 *
 */
void MeshMG::findBoundary() {

	// Setup as two column vectors to find the edges.
	umat e1 = join_vert(join_horiz(elem.col(2), elem.col(1)),
			join_horiz(elem.col(0), elem.col(2)));
	// Calculate all of the values of the edges.
	umat totalEdge = join_vert(e1, join_horiz(elem.col(1), elem.col(0)));
	totalEdge = sort(totalEdge, "ascend", 1);
	vec onev = ones<vec>(totalEdge.n_rows);
	// Create sparse matrix containing the edges which are not being double counted
	// and are thus exterior edges.
	sp_mat fndmat(true, totalEdge.t(), onev, totalEdge.n_rows, totalEdge.n_rows,
			true, true);
	// Initialize edge matrix.
	umat bdEdge = zeros<umat>(totalEdge.n_rows, 2);
	vec s = nonzeros(fndmat);
	uvec ii = zeros<uvec>(s.n_elem);
	uvec jj = zeros<uvec>(s.n_elem);
	uword k = 0;
	// Get the indices for the edges which are exterior edges.
	// The indices map back to nodes.
	for (uword i = 0; i < fndmat.n_rows; i++) {
		for (uword j = 0; j < fndmat.n_cols; j++) {
			double spot = fndmat(i, j);
			// Only find edges where there are two boundary nodes.
			if (spot == 1.0 && k < s.n_elem) {
				ii(k) = i;
				jj(k) = j;
				k++;
			}
		}
	}
	uvec i1 = nonzeros(ii);
	uvec j1 = nonzeros(jj);
	// Setup boolean vector which has whether or not a node is a boundary
	// node.
	isbdNode = zeros<uvec>(N);
	isbdNode.rows(i1) = ones<uvec>(i1.n_rows);
	isbdNode.rows(j1) = ones<uvec>(j1.n_rows);
	isbdNode(0) = 1.0;
	// Find indices of boundary nodes.
	bdNode = find(isbdNode);
	uvec onebd = ones<uvec>(isbdNode.n_rows);
	// Find indices of interior nodes.
	freeNode = find(onebd - isbdNode);
	freeNodes.push_back(freeNode);

}

/**
 * accumArrayM implements the accumarray method from MATLAB.
 *
 * This is a direct implementation and does the exact same thing
 * as the MATLAB function.
 *
 * @param subs - vector of indices.
 * @param ar - vector of values.
 * @param N - size of the array to be made.
 * @return accumulated array (same as MATLAB would).
 */
vec MeshMG::accumArrayM(uvec subs,vec ar,uword N){

	vec S=zeros<vec>(N);

	for(uword i=0;i<N;i++){
		// Get the subscripts.
		uvec q1=find(subs ==(i));
		vec spot;
		double thesum=0.0;
		if(!q1.is_empty()){
			// Find elements at indices q1
			spot=ar.elem(q1);
			// Sum elements.
			thesum=sum(spot);
			// Set at position i.
			S(i)=thesum;
		}
		// If it is empty, add 0.
		else{
			S(i)=0.0;
		}
	}

	return S;
}

/**
 * uniformrefine refines the mesh and re-calculates the boundary
 * as well as all of the needed matrices and vectors.
 *
 * This function ensures that the mesh is refined by adding nodes at
 * the midpoints of all edges in each element.
 *
 */
void MeshMG::uniformrefine(){
	// Find all of the edges.
	umat e23=join_horiz(elem.col(1),elem.col(2));
	umat e31=join_horiz(elem.col(2),elem.col(0));
	umat e12=join_horiz(elem.col(0),elem.col(1));
	// Matrix with all edge indices (n1 -> n2),
	umat totalEdge=sort(join_vert(e23,join_vert(e31,e12)),"ascend",1);
	ArmaFuns funs;
	// Implementation of unique in MATLAB.
	umat edge=funs.uniqueu(totalEdge);
	uvec js=funs.uniqueidx;
	uword NE=edge.n_rows;
	// Map elems to edges.
	umat elem2edge=zeros<umat>(NT*3,1);
	elem2edge.col(0)=js;
	elem2edge.reshape(NT,3);
	// New nodes to create.
	mat nnode=(node.rows(edge.col(0))+node.rows(edge.col(1)))/(2.0);
	// Concatinate the new nodes with the old.
	node=join_vert(node,nnode);
	// Re-map edges including new nodes.
	uvec edge2newNode=regspace<uvec>(N,N+NE-1);
	// Make hierarchical basis.
	HB=join_horiz(edge2newNode,edge);
	uvec t=regspace<uvec>(0,NT-1);
	// Indices for the new elements.
	umat p= join_horiz(elem,
			join_horiz(edge2newNode.rows(elem2edge.col(0)),
			join_horiz(edge2newNode.rows(elem2edge.col(1)),
					edge2newNode.rows(elem2edge.col(2)))));
	// Add elements to the elem matrix.
	elem.rows(t)=join_horiz(p.col(0),join_horiz(p.col(4),p.col(5)));
	// Create new elements.
	umat nelem=zeros<umat>(4*NT,3);
	nelem.rows(span(NT,2*NT-1))=join_horiz(p.col(5),join_horiz(p.col(1),p.col(3)));
	nelem.rows(span(2*NT,3*NT-1))=join_horiz(p.col(4),join_horiz(p.col(3),p.col(2)));
	nelem.rows(span(3*NT,4*NT-1))=join_horiz(p.col(3),join_horiz(p.col(4),p.col(5)));
	nelem.rows(span(0,NT-1))=elem;
	elem=nelem;
	// Save new information.
	NT=elem.n_rows;
	N=node.n_rows;
	refines=refines+1;
	Ns.push_back(N);
	HBs.push_back(HB);
	ProHB();
	vector<sp_mat> stiffmass = assembleMatrix();
	stiffness=stiffmass[0];
	mass=stiffmass[1];
	findBoundary();
}

/**
 * ProHB creates the prolongation and restriction matrices,
 * as well as maps the fine nodes to the coarse nodes for the FAS.
 *
 * This function ensures that the multigrid method works properly and
 * interpolates the coarse grid to the fine grid, as well as restricts
 * the coarse grid and the fine grid.
 *
 */
void MeshMG::ProHB(){
	// Get the number of coarse nodes.
	uword nCoarse = Ns[Ns.size()-2];
	uword nTotal=Ns[Ns.size()-1];
	// Get the number of fine grid only nodes.
	uword nFineNode=nTotal-nCoarse;
	uvec coarseNode=regspace<uvec>(0,nCoarse-1);
	// Save data for use by FEM.
	coarse2Fine.push_back(coarseNode);
	// Setup indices for the prolongation matrix.
	uvec ii=join_vert(coarseNode,join_vert(HB.col(0),HB.col(0)));
	uvec jj=join_vert(coarseNode,join_vert(HB.col(1),HB.col(2)));
	vec ss=join_vert(ones<vec>(nCoarse),
			join_vert(0.5*ones<vec>(nFineNode),0.5*ones<vec>(nFineNode)));
	// Quick sparse construction of the prolongation matrix.
	sp_mat Procurr(true,join_horiz(ii,jj).t(),ss,nTotal,nCoarse);
	Pro.push_back(Procurr);
}

MeshMG::~MeshMG() {

}

