/*
 * ArmaFuns.cpp
 *
 *  Created on: Aug 19, 2016
 *      Author: Devils
 */

using namespace arma;
using namespace std;
#include <map>


#include "ArmaFuns.h"

ArmaFuns::ArmaFuns() {


}

mat ArmaFuns::unique(mat inmat){
	map<vec,uword> unqmap;
	inmat=sort(inmat,"ascend",1);
	unqmap[inmat.row(0)]=0;
	mat outmat=zeros<mat>(1,inmat.n_cols);
	outmat.row(0)=inmat.row(0);
	uniqueidx=zeros<uvec>(inmat.n_rows);
	uniqueidx(0)=0;
	for(uword i=1;i<inmat.n_rows;i++){
		if(unqmap.count(inmat.row(i).t())<1){
			unqmap[inmat.row(i).t()]=i;
			uniqueidx(i)=i;
			outmat=join_vert(outmat,inmat.row(i));
		}
		else{
			uniqueidx(i)=unqmap[inmat.row(i).t()];
		}
	}
	return outmat;
}

umat ArmaFuns::uniqueu(umat inmat){
	map<uvec,uword> unqmap;
	inmat=sort(inmat,"ascend",1);
	unqmap[inmat.row(0)]=0;
	umat outmat=zeros<umat>(1,inmat.n_cols);
	outmat.row(0)=inmat.row(0);
	uniqueidx=zeros<uvec>(inmat.n_rows);
	uniqueidx(0)=0;
	for(uword i=1;i<inmat.n_rows;i++){
		if(unqmap.count(inmat.row(i).t())<1){
			unqmap[inmat.row(i).t()]=i;
			uniqueidx(i)=i;
			outmat=join_vert(outmat,inmat.row(i));
		}else{
			uniqueidx(i)=unqmap[inmat.row(i).t()];
		}
	}
	return outmat;
}

ArmaFuns::~ArmaFuns() {
	// TODO Auto-generated destructor stub
}

