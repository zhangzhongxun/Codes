#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <algorithm>
#include <valarray>
#include <ctime>
#include <cstdio>
#include <cstdlib>

#include "MatVec.h"
#include "IO.h"
#define PI 3.14159265359


//function to calculate the givens transform for given two numbers
//a, b are given numbers, c is cosine, s is sine, r is the nonzero
//entry after the Givens transform
void Givens(double a, double b, double& c, double& s, double& r){

  double tau;

  if( b == 0) 
    {
      c = 1.0; 
      s = 0.0;
    }
  else
    {
      if(fabs(b) > fabs(a))
	{
	  tau = -a/b;
	  s = 1.0/sqrt(1+tau*tau);
	  c = s*tau;
	}
      else
	{
	  tau = -b/a;
	  c = 1.0/sqrt(1+tau*tau);
	  s = c*tau;
	}
    }

  r = c*a - s*b;

}    


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

void GMRES(const MatMulVec& A, const valarray<double>& b, valarray<double>& x, double tol, int max_k,
           int& flag, double& res, int& num_restart , int& num_k){

  int m, n, m_tmp, i_restart, i, j, k, seed;
  double h10, h_kpk, h_tmp; //here h_kpk denotes h_(k+1,k)
  char debfile[80];

  if((n = A.get_dim()) != b.size()) // A is a nxn matrix, b is a vector of size n
    {
      cout << endl << "Dimension incompatible in GRMES() function,ERROR in LinearSolver.cc !!! " << endl;
      exit(-1);
    }

  m_tmp = n - 1; // = (int)(sqrt(1.0*n));
  m = m_tmp < max_k ? m_tmp : max_k; // do at most m steps before restart
  UpHessMat H(m);              // H is the upper Hessenberg Matrix 
   
  valarray<double> r(0.0, n); //the residual vector
//   valarray<double> x0(0.0, n); // initial value = 0 at this moment
  h10 = 0.0;

  //set initial value to be random
//   seed = time(0);
//   srand(seed);
//   for(i = 0; i < x0.size(); i++)
//     x0[i] = 1e-8*rand()/double(RAND_MAX);

  valarray<double> x0(b); // initial value, = RHS b at this moment
  A.get_diag(r); //get the diagonal vector of A, put it in r

//   if(Mini_abs_element(r) < 1e-15)  //get the smallest absolute value of r, prevent singular case
//     {
//       cout << endl << "The matrix A in GMRES solver has zero diagonal element, can not set x0 = diag(A)\b !!!" << endl;
//     }
//   else
//     {
//       x0 /= r;
//     }  

  if(Mini_abs_element(r) >= 1e-15)
    {
      x0 /= r;
    }
 
  //the vector arrays contain the q vectors in GMRES algorithm
//   valarray< valarray<double>* > q(m);
  valarray< valarray<double> > q(m);

  for(i = 0; i < m; i++)
      q[i].resize(n);
//     q[i] = new valarray<double>(n);

//    //***************************DEBUG**********************************
//     ofstream ofile("/var/tmp/2D/GMRESdeb.txt",ios::out);
//     ofile << endl << "__________________________________________________________________________________" 
//           << endl << endl << "In the GMRES solver, tol = " <<  tol << ", m = " << m << endl;
//     //***************************DEBUG**********************************

  i_restart = 0;
  while(i_restart < n) //we allow at most n times restart
    {
//       //***************************DEBUG**********************************
//       ofile << endl << " i_restart = " << i_restart << ", n = " << n << ", H.get_curr() = " << H.get_curr();
//       //***************************DEBUG**********************************

      if(i_restart > 0) //it is a restart, need to get x0
	{
	  valarray<double> y(0.0,m);
	  valarray<double> rho(0.0,m+1);
	  rho[0] = h10;
	  H.LS_solve(rho, y);
	  for(i = 0; i < m; i++)
	    x0 += (q[i])*y[i];

// 	  //***************************DEBUG**********************************
//            sprintf(debfile,"/var/tmp/2D/gmres_y_%d.txt",i_restart);
//           Write_Vector(debfile,y,y.size());
// 	  //***************************DEBUG**********************************
	}

      A.Multiply(x0,r); //r = Ax0
      r = b - r;
      h10 = sqrt(InnerProd(r,r));
      h_kpk = h10; //

//       //***************************DEBUG**********************************
//       ofile << ", h10 = " << h10 << endl;      
//       //***************************DEBUG**********************************

      k = 0;
      H.set_curr(k);
      while(k < m && h_kpk > tol) //doesn't reach m steps and doesn't converge
	{
// 	  //***************************DEBUG**********************************
//           ofile << endl << "k = " << k << ", curr = " << H.get_curr() 
//                 << ", h_kpk = " << h_kpk << endl;
// 	  //***************************DEBUG**********************************

	  q[k] = r/h_kpk;  //get the k-th column of Q

	  k += 1;

	  A.Multiply(q[k-1],r); //r = A*q[k-1]
        
	  for(i = 0; i < k; i++)
	    { 
 // 	      //***************************DEBUG**********************************
//               h_tmp  = sqrt(InnerProd(r,r));
// 	      ofile << endl << " i = " << i  << ", ||r|| = " << h_tmp;
// 	      //***************************DEBUG**********************************

	      h_tmp  = InnerProd(q[i],r);
	      H.get(i,k-1) = h_tmp;
	      r -= h_tmp*(q[i]);

// 	      //***************************DEBUG**********************************
// 	      ofile << ", H[" << i << "," << k -1 << "] = " << h_tmp << endl;
// 	      //***************************DEBUG**********************************
	    }

	  h_kpk = sqrt(InnerProd(r,r));
	  H.get(k,k-1) = h_kpk;
	  H.set_curr(k);

//           //***************************DEBUG**********************************
// 	  ofile << " h_kpk = " << h_kpk << endl;
//           //***************************DEBUG**********************************
	}

      if(h_kpk <= tol) //converged !
	{
	  flag = 0;
	  break;
	}
      //else, it doesn't converge, restart with updated x0
      i_restart += 1;
    }

  //solve for x anyway

  m_tmp = H.get_curr(); //get the size of the LS problem

  //only need to solve if the solution is non-trivial, i.e. x != x0 (i.e., H.get_curr() > 0)
  //also, if H.get_curr = m, k reaches its max-value, and x0 is already calcuted
  
  if((i_restart < n && m_tmp > 0) || (i_restart == n)) 
    {
      valarray<double> y(0.0,m_tmp);
      valarray<double> rho(0.0,m_tmp+1);
      rho[0] = h10;

      H.LS_solve(rho, y);

      for(i = 0; i < m_tmp; i++)
	x0 += (q[i])*y[i]; 


      A.Multiply(x0,r); //r = Ax0
      r = b - r;
      h10 = sqrt(InnerProd(r,r));
//       //***************************DEBUG**********************************
//       cout << endl << "At the end of GMRES y solve, i_restart = " << i_restart << ", h10 = " << h10 << endl;
//       sprintf(debfile,"/var/tmp/2D/gmres_y_%d",i_restart);
//       Write_Vector(debfile,y,y.size());
//       //***************************DEBUG**********************************
    }

  x = x0;       
   
  if(i_restart == n && h10 > tol)
    {
      flag = -1;
    }

  num_restart = i_restart;
  res = h10;
  num_k = H.get_curr();
 
//   delete the q vectors
//   for(i = 0; i < m; i++)
//     delete q[i];

//   //***************************DEBUG**********************************
//   ofile.close();
//   //***************************DEBUG**********************************
}


//function to solve the Poisson equation with Neumann BC on a unit rectangle, the number of grid point
//in y direction is m+1, where m has to be power of 2, i.e., m = 2^(k+1), the number of grid points 
//in x direction is n, and the RHS is given as an object y of MultipleRHS, with col_num = m+1, and dim = n
//after the solve, the solutoin x overrides y
//Here the RHS must be multiplied by (dx)^2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

void Neumann_Poisson_Solver(int k, int n, MultipleRHS& y){

  int i, j, m, r, shift;

  m = 1 << (k+1); //m = 2^(k+1)

  if(m+1 != y.get_col_num() || n != y.get_dim()) 
    {
      cout << endl << "Dimension Incompatible in Neumann_Poisson_Solver()"
	   << endl << "ERROR in LinearSolver.cc!! " << endl;
      exit(-1);
    }

  //define p and q vectors, note the initial value is 0
  MultipleRHS p(m+1,n);
  MultipleRHS q(m+1,n);

  //define matrix A
  valarray<double> vd(-4.0,n);
  valarray<double> vsd(1,n-1);
 
  valarray<double> T(1.0, n);

  TriDiag A(n,0); //A is  non-singular, with size n

  A.get_d() = vd;
  A.get_e() = vsd;
  A.get_f() = vsd;

  A.get_e()[n-2] = 2.0;
  A.get_f()[0] = 2.0;

  TriDiag Ar(n, 0); //non-singular tridiagonal, each only the diagonal change
  Ar.get_e() = A.get_e();
  Ar.get_f() = A.get_f();

  //compute p and q vectors
  //compute p^0 and q^0, note they all ocupy the same p and q
  //p^0_j = 0; q^0_j = y_j
  for(j = 0; j <=m; j++)
    q[j] = y[j];

  //compute p^1 and q^1, put in { } to make the variable tmp_p local
  {
    MultipleRHS tmp_p((1<<k)+1, n); // size = 2^k + 1
    for(j = 0; j <= m; j += 2)
      tmp_p[j/2] = y[j];

    //solve the tridiagonal systems
    A.ChaseSolve(tmp_p);

    p[0] = tmp_p[0];
    q[0] = 2.0*(y[1] - p[0]);

    for(j = 2; j <= m-2; j += 2)
      {
	p[j] = tmp_p[j/2];
        q[j] = y[j-1] + y[j+1] - 2.0*p[j];
      }

    p[m] = tmp_p[m/2];
    q[m] = 2.0*(y[m-1] - p[m]);

  }


  for(r = 1; r <= k-1; r++)
    {
      shift = (1 << r);        //shift = 2^r
      MultipleRHS tmp_p((1 << (k-r) ) + 1, n); 
  
      tmp_p[0] = 2.0*p[shift] - q[0];
      for(i = 1; i <= (1 << (k-r)) - 1; i++)
        {
	  j = i*(1 << (r+1));
 
          tmp_p[i] = p[j-shift] + p[j+shift] - q[j];
	}
      tmp_p[1 << (k-r)] = 2.0*p[m-shift] - q[m];

      //solve for A^(r)^(-1) tmp_p
      for(j = 1; j <= (1 << r); j++)
	{
	  double theta = (2.0*j - 1.0)*(PI)/(1 << (r+1));
          
          Ar.get_d() = A.get_d() + 2.0*cos(theta)*T;

          Ar.ChaseSolve(tmp_p);
	}

      //now update p and q
      p[0] += tmp_p[0];
      q[0]  = 2.0*(q[shift] - p[0]);
      for(i = 1; i <= (1 << (k-r)) - 1; i++)
        {
	  j = i*(1 << (r+1));
          p[j] += tmp_p[i];
          q[j]  = q[j-shift] + q[j+shift] - 2.0*p[j];
	}
      p[m] += tmp_p[1 << (k-r)];
      q[m]  = 2.0*(q[m-shift] - p[m]);
    }

  //r = k
  {
    MultipleRHS tmp_p(1,n);

    tmp_p[0] = p[0] + p[m] - q[1 << k];

    for(j = 1; j <= (1 << k); j++)
      {
       double theta = (2.0*j - 1.0)*(PI)/(1 << (k+1));
       Ar.get_d() = A.get_d() + 2.0*cos(theta)*T;

       Ar.ChaseSolve(tmp_p);
      }

    p[1 << k] += tmp_p[0];
    q[1 << k]  = q[0] + q[m] - 4.0*p[1 << k];

  }


  //now solve for x(2^k)
  {
    MultipleRHS tmp_p(1,n);
    tmp_p[0] = -q[1 << k];

    for(j = 1; j <= (1<< (k+1))-1; j++)
      {
	double theta = (1.0*j)*(PI)/(1 << k); 
             
        Ar.get_d() = A.get_d() + 2.0*cos(theta)*T; //Ar is still non-singular

        Ar.ChaseSolve(tmp_p);
      }

    // j = 2^(k+1)
    TriDiag As(n,1); //As is singular

    As.get_e() = A.get_e();
    As.get_f() = A.get_f();

    As.get_d() = A.get_d() + 2.0*T;

    As.ChaseSolve(tmp_p);
 
    //get x(2^k)
    y[1 << k] = p[1 << k] + tmp_p[0];
  }



  //now do the back-substitution
  //r = k
  {
    MultipleRHS tmp_p(2,n);
 
    tmp_p[0] = 2.0*y[1 << k] - q[0];
    tmp_p[1] = 2.0*y[1 << k] - q[m];

    for(j = 1; j <= (1 << k); j++)
      {
	double theta = (2.0*j - 1.0)*(PI)/(1 << (k+1));
	Ar.get_d() = A.get_d() + 2.0*cos(theta)*T;

	Ar.ChaseSolve(tmp_p);
      }

    y[0] = p[0] + tmp_p[0];
    y[m] = p[m] + tmp_p[1];

  }

  //for r = k -1 : -1 : 1
  for(r = k - 1; r >= 1; r--)
    {
      shift = (1 << r); // = 2^r
      MultipleRHS tmp_p(1 << (k-r), n);

      for(i = 1; i <= ( 1 << (k + 1 - r)) - 1; i += 2)
	{
	  j = i*shift;
          tmp_p[(i-1)/2] = y[j-shift] + y[j+shift] - q[j];
	}

      //solve for A^(r)^(-1) tmp_p
      for(j = 1; j <= (1 << r); j++)
	{
	  double theta = (2.0*j - 1.0)*(PI)/(1 << (r+1));
          
          Ar.get_d() = A.get_d() + 2.0*cos(theta)*T;

          Ar.ChaseSolve(tmp_p);
	}

      //get x at level r
      for(i = 1; i <= ( 1 << (k + 1 -r)) - 1; i += 2)
	{
	  j = i*shift;
          y[j] = p[j] + tmp_p[(i-1)/2];
	}
      
    }

   //for r = 0, same as before, treat it specially for efficiency
  {
    MultipleRHS tmp_p(1 << k, n);

    for(i = 1; i <= ( 1 << (k+1) ) - 1; i += 2)
       tmp_p[(i-1)/2] = y[i] - y[i-1] -  y[i+1];

    A.ChaseSolve(tmp_p);

    for(i = 1; i <= ( 1 << (k+1) ) - 1; i += 2)
      y[i] =  tmp_p[(i-1)/2];
  }

}
      
        



      

     


  

  
  
