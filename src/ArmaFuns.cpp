/*
 * ArmaFuns.cpp
 *
 *  Created on: Aug 19, 2016
 *      Author: Devils
 */

#include "ArmaFuns.h"

using namespace arma;
using namespace std;
#include <map>
#include <string>


#include "ArmaFuns.h"

ArmaFuns::ArmaFuns() {


}

mat ArmaFuns::unique(mat inmat){
	map<string,uword> unqmap;
	inmat=sort(inmat,"ascend",1);
	string mapstr = to_string(inmat(0,0)) + "," + to_string(inmat(0,1));
	unqmap[mapstr]=0;
	mat outmat=zeros<mat>(1,inmat.n_cols);
	outmat.row(0)=inmat.row(0);
	uniqueidx=zeros<uvec>(inmat.n_rows);
	uniqueidx(0)=0;
	for(uword i=1;i<inmat.n_rows;i++){
		mapstr = to_string(inmat(i,0)) + "," + to_string(inmat(i,1));
		if(unqmap.count(mapstr)<1){
			unqmap[mapstr]=i;
			uniqueidx(i)=i;
			outmat=join_vert(outmat,inmat.row(i));
		}
		else{
			uniqueidx(i)=unqmap[mapstr];
		}
	}
	return outmat;
}

umat ArmaFuns::uniqueu(umat inmat){
	map<string,uword> unqmap;
	inmat=sort(inmat,"ascend",1);
	string mapstr = to_string(inmat(0,0)) + "," + to_string(inmat(0,1));
	unqmap[mapstr]=0;
	umat outmat=zeros<umat>(1,inmat.n_cols);
	outmat.row(0)=inmat.row(0);
	uniqueidx=zeros<uvec>(inmat.n_rows);
	uniqueidx(0)=0;
	for(uword i=1;i<inmat.n_rows;i++){
		mapstr = to_string(inmat(i,0)) + "," + to_string(inmat(i,1));
		if(unqmap.count(inmat.row(i).t())<1){
			unqmap[mapstr]=i;
			uniqueidx(i)=i;
			outmat=join_vert(outmat,inmat.row(i));
		}else{
			uniqueidx(i)=unqmap[mapstr];
		}
	}
	return outmat;
}

ArmaFuns::~ArmaFuns() {
	// TODO Auto-generated destructor stub
}

