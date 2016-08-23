//============================================================================
// Name        : cFEM-Nonlinear.cpp
// Author      : Ted Kwan
// Version     :
// Copyright   : Your copyright notice
// Description : Nonlinear Finite Element Method in C++
//============================================================================

#include <iostream>
#include "Mesh.h"
#include "MeshMG.h"
#include <vector>
#include <armadillo>
#include "FiniteElemNL.h"

using namespace std;
using namespace arma;

/**
* This is a testing program for the nonlinear finite element method
* in CPP.
*/
int main() {
	vector<double> meshvec(5);
	meshvec[0]=0;
	meshvec[1]=1.0;
	meshvec[2]=0;
	meshvec[3]=1.0;
	meshvec[4]=0.25;
	//Mesh msh(meshvec,"blahdiddyblah");
	clock_t start;
	start=clock();
	FiniteElemNL femp(meshvec,"blah man yo",3);
	double duration=clock()-start;
	cout << "This took: " << (duration/(double)(CLOCKS_PER_SEC))<< endl;
	//vec u=femp.u;
	//u.print("u: ");
//	MeshMG mshmg(meshvec);
//	mshmg.uniformrefine();
	//vector<sp_mat> stiffmass=mshmg.assembleMatrix();
	//sp_mat stiffness=stiffmass[0];
	//	mass.set_size(stiffmass[1].n_rows,stiffmass[1].n_cols);
	//sp_mat mass=stiffmass[1];
//	sp_mat stiffness=mshmg.stiffs.back();
//	stiffness.print("Stiffness: ");
	//cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!
	return 0;
}
