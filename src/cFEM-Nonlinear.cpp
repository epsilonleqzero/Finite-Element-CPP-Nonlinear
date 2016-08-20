//============================================================================
// Name        : cFEM-Nonlinear.cpp
// Author      : Ted Kwan
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include "Mesh.h"
#include "MeshMG.h"
#include <vector>
#include <armadillo>
#include "FiniteElemNL.h"

using namespace std;
using namespace arma;

int main() {
	vector<double> meshvec(5);
	meshvec[0]=0;
	meshvec[1]=1.0;
	meshvec[2]=0;
	meshvec[3]=1.0;
	meshvec[4]=0.25;
	//Mesh msh(meshvec,"blahdiddyblah");
	FiniteElemNL femp(meshvec,"blah man yo");
//	vec u=femp.u;
//	u.print("u: ");
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
