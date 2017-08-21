#ifndef _LINEARSOLVER_H
#define _LINEARSOLVER_H
#include "MatVec.h"

using namespace std;

// This file contains the functions to solve a linear system by GMRES method
//function to calculate the givens transform for given two numbers
//a, b are given numbers, c is cosine, s is sine, r is the nonzero
//entry after the Givens transform
void Givens(double a, double b, double& c, double& s, double& r);

//function to solve a linear system by GMRES method
//Arguments: 
//A: a reference of the matrix object (derived from class MatMulVec
//with function Multiply overrode)
//b: right hand side
//tol: the tolerance for convergence
//Return values:
// x: the solution
// flag: converge = 0, didn't converge = -1
// num_restart: number of restarts
// num_k: the number of q-vectors generated in the last restart

void GMRES(const MatMulVec& A, const valarray<double>& b, valarray<double>& x, double tol, int max_k,
           int& flag, double& res, int& num_restart, int& num_k);

//function to solve the Poisson equation with Neumann BC on a rectanble, the number of grid point
//in y direction is m+1, where m has to be power of 2, i.e., m = 2^(k+1), the number of grid points 
//in x direction is n, and the RHS is given as an object y of MultipleRHS, with col_num = m+1, and dim = n
//after the solve, the solutoin x overrides y
//Here the RHS must be multiplied by (dx)^2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

void Neumann_Poisson_Solver(int k, int n, MultipleRHS& y);

#endif
