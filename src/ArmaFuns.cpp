/*
 * ArmaFuns.cpp
 *
 *  Created on: Aug 19, 2016
 *      Author: Devils
 */

#include "ArmaFuns.h"
#include <iostream>
#include <map>
#include <string>

using namespace arma;
using namespace std;


#include "ArmaFuns.h"

/**
 * Default construction, nothing to do.
 *
 */
ArmaFuns::ArmaFuns() {

}

/**
 * unique implements the unique method from MATLAB.
 *
 * This is a direct implementation and does the exact same thing
 * as the MATLAB function.
 *
 * @param inmat - matrix to run the method on.
 * @return matrix containing unique indices and values
 * (same as MATLAB would).
 */
mat ArmaFuns::unique(mat inmat){
	map<string,uword> unqmap;
	// Sort matrix by rows.
	inmat=sort(inmat,"ascend",1);
	// Setup map to use for checking for unique entries.
	string mapstr = to_string(inmat(0,0)) + "," + to_string(inmat(0,1));
	unqmap[mapstr]=0;
	mat outmat=zeros<mat>(1,inmat.n_cols);
	outmat.row(0)=inmat.row(0);
	uniqueidx=zeros<uvec>(inmat.n_rows);
	uniqueidx(0)=0;
	for(uword i=1;i<inmat.n_rows;i++){
		// Find current row as a string representation.
		mapstr = to_string(inmat(i,0)) + "," + to_string(inmat(i,1));
		// Check if it is unique.
		if(unqmap.count(mapstr)<1){
			// Store unique row.
			unqmap[mapstr]=i;
			uniqueidx(i)=i;
			outmat=join_vert(outmat,inmat.row(i));
		}
		else{
			// Map back non-unique row to the unique row.
			uniqueidx(i)=unqmap[mapstr];
		}
	}
	return outmat;
}

/**
 * unique implements the unique method from MATLAB.
 *
 * This is a direct implementation and does the exact same thing
 * as the MATLAB function. This method is for use on index matrices.
 *
 * @param inmat - matrix to run the method on.
 * @return matrix containing unique indices and values
 * (same as MATLAB would).
 */
umat ArmaFuns::uniqueu(umat inmat){
	map<string,uword> unqmap;
	// Sort matrix by rows.
	inmat=sort(inmat,"ascend",1);
	// Setup map to use for checking for unique entries.
	string mapstr = to_string(inmat(0,0)) + "," + to_string(inmat(0,1));
	//cout << mapstr << endl;
	unqmap[mapstr]=0;
	umat outmat=zeros<umat>(1,inmat.n_cols);
	outmat.row(0)=inmat.row(0);
	uniqueidx=zeros<uvec>(inmat.n_rows);
	uniqueidx(0)=0;
	uword k=1;
	for(uword i=1;i<inmat.n_rows;i++){
		// Find current row as a string representation.
		mapstr = to_string(inmat(i,0)) + "," + to_string(inmat(i,1));
		// Check if it is unique.
		if(unqmap.count(mapstr)<1){
			// Store unique row.
			unqmap[mapstr]=k;
			uniqueidx(i)=k++;
			outmat=join_vert(outmat,inmat.row(i));
		}else{
			// Map back non-unique row to the unique row.
			uniqueidx(i)=unqmap[mapstr];
		}
	}
	return outmat;
}

ArmaFuns::~ArmaFuns() {
	
}

