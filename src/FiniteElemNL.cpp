/*
 * FiniteElemNL.cpp
 *
 * This class runs the finite element method and obtains a
 * numerical solution to a PDE on a domain given into the constructor.
 *
 * Currently, only two functions have been implemented, and square domains.
 *
 *  Created on: Aug 13, 2016
 *      Author: Ted Kwan
 */

#include "FiniteElemNL.h"
#include "FunZeros.h"
#include "FunOnes.h"
#include "FunSC.h"
#include "FunSCNL.h"
#include "FunSCbD.h"
#include "FunBoltz.h"
#include <iostream>

using namespace std;
using namespace arma;

FiniteElemNL::FiniteElemNL() {

}

/**
 * Constructor to create and run.
 *
 * This constructor creates and runs the finite element method using the properties
 * provided in meshprops and using the name provided in fun.
 *
 * @param meshprops - Vector containing the information on the square.
 * @param fun - string containing the name of the function to use.
 * @param j - integer representing the number of refinements and multigrid levels.
 *
 */
FiniteElemNL::FiniteElemNL(vector<double> meshprops,string fun, int j) {
	// Create mesh.
	MeshMG msh(meshprops);
	J=j;
	// Uniform refine the mesh J times.
	for(uword i=0;i<J;i++){
		msh.uniformrefine();
	}
	// Choose and apply function using polymorphism.
	pde=new FunZeros;
	bdfun=new FunBoltz;
	pdenl=new FunSCNL;
	// Set the mesh variable to the mesh created.
	mesh=msh;
	uword N=mesh.node.n_rows;
	// Initialize solution vector.
	u=zeros<vec>(N);
	// Calculate the right hand side of the equation.
	vec b=calcRHS();
	mat node=mesh.node;
	// Calculate the boundary condition.
	u.rows(mesh.bdNode)=bdfun->evalF(node.rows(mesh.bdNode));
	// Setup params to use for FAS.
	double tol=1e-6;
	// Global stiffness and mass matrices.
	mat A(mesh.stiffness);
	mat M(mesh.mass);
	vec Mv=M.diag();
	vec fu=pdenl->evalF(u);
	// Setup system to solve.
	vec residual=((A*u)+(Mv%fu)-b);
	uvec freeNode=mesh.freeNode;
	vec r=residual.rows(freeNode);
	uword k=0;
	uword maxitr =20;
	tol=tol*norm(r);
	double err=2*tol;
	vec u0=u;
	while(k<50 && err>tol){
		// Run the FAS method and solve.
		u=MgFAS(b,u,maxitr,1e-6,J);
		fu=pdenl->evalF(u);
		residual=((A*u)+(Mv%fu)-b);
		r=residual.rows(freeNode);
		err=norm(r);
		k++;
	}
	cout << k << " iterations." << endl;
	// Release the function classes created.
	delete pde;
	delete bdfun;
	delete pdenl;
}

/**
 * Constructor to create and run.
 *
 * This constructor creates and runs the finite element method using the properties
 * provided in meshprops and using the name provided in fun.
 *
 * @param meshprops - Vector containing the information on the square.
 * @param fun - string containing the name of the function to use.
 *
 */
FiniteElemNL::FiniteElemNL(vector<double> meshprops,string fun) {
	// Create mesh.
	MeshMG msh(meshprops);
	// Choose and apply function using polymorphism.
	pde=new FunZeros;
	bdfun=new FunBoltz;
	pdenl=new FunSCNL;
	// Set the mesh variable to the mesh created.
	mesh=msh;
	uword N=mesh.node.n_rows;
	// Initialize solution vector.
	u=zeros<vec>(N);
	// Calculate the right hand side of the equation.
	vec b=calcRHS();
	mat node=mesh.node;
	// Calculate the boundary condition.
	u.rows(mesh.bdNode)=bdfun->evalF(node.rows(mesh.bdNode));
	vec br=b-(mesh.stiffness*u);
	double tol=1e-6;
	mat A(mesh.stiffness);
	mat M(mesh.mass);
	vec Mv=M.diag();
	vec fu=pdenl->evalF(u);
	vec residual=((A*u)+(Mv%fu)-b);
	uvec freeNode=mesh.freeNode;
	vec r=residual.rows(freeNode);
	//uword k=0;
	//uword maxitr=20;
	tol=tol*norm(r);
	//double err=2*tol;
	u=NWTsolve(b, u, 10, Mv, A, 1e-6, mesh.freeNode);
	delete pde;
	delete bdfun;
	delete pdenl;
}

/**
 * calcRHS() is an internal method to find the right hand side.
 *
 * This is used above to set b for the system Au=b
 */
vec FiniteElemNL::calcRHS(){
	// Get needed Matrices and vectors.
	vec area=mesh.area; mat node=mesh.node; umat elem=mesh.elem;
	// Calculate Midpoints.
	mat mid1=(node.rows(elem.col(1))+node.rows(elem.col(2)))/2.0;
	mat mid2=(node.rows(elem.col(2))+node.rows(elem.col(0)))/2.0;
	mat mid3=(node.rows(elem.col(0))+node.rows(elem.col(1)))/2.0;
	// Calculate vectors to accumulate
	vec bt1=(area%((pde->evalF(mid2))+(pde->evalF(mid3))))/6.0;
	vec bt2=(area%((pde->evalF(mid3))+(pde->evalF(mid1))))/6.0;
	vec bt3=(area%((pde->evalF(mid1))+(pde->evalF(mid2))))/6.0;
	// Vertically concatenate vectors.
	vec bt12=join_vert(bt1,bt2);
	vec bts=join_vert(bt12,bt3);
	uvec elemv=vectorise(elem);
	// Accumulate the array.
	vec ret=accumArray(elemv,bts,node.n_rows);
	return ret;
}

/**
 * accumArray implements the accumarray method from MATLAB.
 *
 * This is a direct implementation and does the exact same thing
 * as the MATLAB function.
 *
 * @param subs - vector of indices.
 * @param ar - vector of values.
 * @param N - size of the array to be made.
 * @return accumulated array (same as MATLAB would).
 */
vec FiniteElemNL::accumArray(uvec subs,vec ar,uword N){

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
 * NWTsolve implements Newton's method for nonlinear equations.
 *
 * This method runs Newton's method to solve the nonlinear
 * system for the FAS method. It is used on the coarsest grid.
 *
 * @param b - vector for right hand side (or residual).
 * @param u - vector of initial guess, and return vector.
 * @param maxitr - maximum number of iterations.
 * @param M - vector of mass matrix.
 * @param A - stiffness matrix.
 * @param tol - double for the tolerance to solve to.
 * @param freeNode - vector containing non-boundary node indices.
 * @return u - solution to the nonlinear system.
 */
vec FiniteElemNL::NWTsolve(vec b,vec u, uword maxitr,
		vec M, mat A ,double tol,uvec freeNode){
	// Initial setup
	vec fu=pdenl->evalF(u);
	vec residual=((A*u)+(M%fu)-b);
	vec r=residual.rows(freeNode);
	uword k=0;
	tol=tol*norm(r);
	double err=2*tol;
	while(k<maxitr && err>tol){
		// Run Newton's method.
		vec Df=M%(pdenl->evalDF(u));
		mat Dfm=diagmat(Df);
		// Setup matrix.
		mat B=A+Dfm;
		// Solve the error equation.
		vec e=solve(B.submat(freeNode, freeNode), r);
		u.rows(freeNode)=u.rows(freeNode)-e;
		fu=pdenl->evalF(u);
		// Calculate residual.
		residual=((A*u)+(M%fu)-b);
		r=residual.rows(freeNode);
		err=norm(r);
		k=k+1;
//		cout << "Iteration: " << k << " error: " << err << endl;
	}
	return u;

}

/**
 * GSsolve implements Gauss-Seidel's method for nonlinear systems.
 *
 * This method runs the Gauss-Seidel method to solve the nonlinear
 * system for the FAS method. It is used as a pre-smoother.
 *
 * @param b - vector for right hand side (or residual).
 * @param u - vector of initial guess, and return vector.
 * @param maxitr - maximum number of iterations.
 * @param M - vector of mass matrix.
 * @param A - stiffness matrix.
 * @param tol - double for the tolerance to solve to.
 * @param freeNode - vector containing non-boundary node indices.
 * @return u - solution to the nonlinear system.
 */
vec FiniteElemNL::GSsolve(vec b,vec u, uword maxitr,
		vec M, mat A ,double tol,uvec freeNode){
	uword Nu=u.n_rows;
	for(uword i=1;i<Nu-1;i++){
		uvec isfree=find(freeNode==i);
		// Only solve on non-boundary nodes.
		if(!isfree.empty()){
			// Solve the nonlinear system.
			mat Ar=A.row(i);
			mat currA=(Ar.cols(span(0,i-1))*u.rows(span(0,i-1)))+
						(Ar.cols(span(i+1,Nu-1))*u.rows(span(i+1,Nu-1)));
			double ci=currA(0,0)-b(i);
			vector<double> pars(3);
			pars[0]=ci;
			pars[1]=A(i,i);
			pars[2]=M(i);
			FunSCNL1D nfun(pars);
			// Solve the 1-D equation.
			u(i)=NWTsolve1D(nfun,u(i),1e-6,10);
		}
	}
	return u;
}

/**
 * GSsolveb implements Gauss-Seidel's method (in reverse) 
 * for nonlinear systems.
 *
 * This method runs the backwards Gauss-Seidel method to solve 
 * the nonlinear system for the FAS method. It is used as a 
 * post-smoother.
 *
 * @param b - vector for right hand side (or residual).
 * @param u - vector of initial guess, and return vector.
 * @param maxitr - maximum number of iterations.
 * @param M - vector of mass matrix.
 * @param A - stiffness matrix.
 * @param tol - double for the tolerance to solve to.
 * @param freeNode - vector containing non-boundary node indices.
 * @return u - solution to the nonlinear system.
 */
vec FiniteElemNL::GSsolveb(vec b,vec u, uword maxitr,
		vec M, mat A ,double tol,uvec freeNode){
	uword Nu=u.n_rows;
	for(uword i=Nu-2;i>=1;i--){
		uvec isfree=find(freeNode==i);
		// Only solve on non-boundary nodes.
		if(!isfree.empty()){
			// Solve the nonlinear system.
			mat Ar=A.row(i);
			mat currA=(Ar.cols(span(0,i-1))*u.rows(span(0,i-1)))+
					  (Ar.cols(span(i+1,Nu-1))*u.rows(span(i+1,Nu-1)));
						//double ci=currA(0,0)+b(i);
			double ci=currA(0,0)-b(i);
			vector<double> pars(3);
			pars[0]=ci;
			pars[1]=A(i,i);
			pars[2]=M(i);
			FunSCNL1D nfun(pars);
			// Solve the 1-D equation.
			u(i)=NWTsolve1D(nfun,u(i),1e-6,10);
		}
	}
	return u;
}

/**
 * NWTsolve1D implements Newton's method for nonlinear equations.
 *
 * This method runs Newton's method for nonlinear equations in 1D.
 * It is used as part of the Gauss-Seidel methods.
 *
 * @param f - Function of equation to solve.
 * @param x0 - initial guess for the solution.
 * @param maxitr - maximum number of iterations.
 * @param tol - double for the tolerance to solve to.
 * @return x - solution to the nonlinear equation.
 */
double FiniteElemNL::NWTsolve1D(FunSCNL1D f, double x0, double tol, uword maxitr){
	double fp=f.evalF(x0);
	double df=f.evalDF(x0);
	uword n=0;
	double err=abs(fp);
	//cout << "Newton time" << endl;
	double x=x0;
	while(err>tol && n<maxitr){
		x=x-(fp/df);
		df=f.evalDF(x);
		fp=f.evalF(x);
		err=abs(fp);
		n++;
	}
	return x;
}

/**
 * MgFAS implements full approximation scheme for nonlinear finite
 * element systems. 
 *
 * This method runs the backwards FAS scheme to solve 
 * the nonlinear system for the Poisson-Boltzmann equation.
 *
 * @param b - vector for right hand side (or residual).
 * @param u - vector of current guess, and return vector.
 * @param maxitr - maximum number of iterations.
 * @param tol - double for the tolerance to solve to.
 * @param level - number of refinements on grid.
 * @return u - approximated solution to the PDE.
 */
vec FiniteElemNL::MgFAS(vec b,vec u, uword maxitr,double tol,uword level){

	if(level==0){
		// Base case.
		uvec freeNode=mesh.freeNodes[0];
		vec M=mesh.masses[0];
		mat A(mesh.stiffs[0]);
		// Directly solve on coarse grid.
		u=NWTsolve(b,u,15,M,A,1e-8,freeNode);
		return u;
	}
	else{
		uvec freeNode=mesh.freeNodes[level];
		vec M=mesh.masses[level];
		sp_mat stiff=(mesh.stiffs[level]);
		mat A(stiff);
		sp_mat Pro=(mesh.Pro[level]);
		sp_mat Res=Pro.t();
		// Pre-Smoothing
		u=GSsolve(b, u, 15, M, A, 1e-8, freeNode);

		// MG
		vec fp=pdenl->evalF(u);
		// Restrict to coarse grid.
		vec rHc=Res*(b-(A*u)-(M%fp));
		vec vc=u.rows(mesh.coarse2Fine[level]);
		vec fpvc=pdenl->evalF(vc);
		// Correction for restriction to coarse grid.
		vec rc=(mesh.stiffs[level-1]*vc)+(mesh.masses[level-1]%fpvc);
		// Recursion step
		vec uc=MgFAS(rHc+rc,vc,maxitr,tol,level-1);
		// Prolongate error equation to fine grid.
		vec eh=Pro*(uc-vc);
		// Correction step.
		u=u+eh;
		// Pre-Smoothing
		u=GSsolveb(b, u, 15, M, A, 1e-8, freeNode);
		return u;
	}
}

FiniteElemNL::~FiniteElemNL() {
}

