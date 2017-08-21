#include "fish_subroutine.h"
#include <cmath>

//*****************************************************************************************************************
 int my_hwscrt(double a, double b, int m, int mbdcnd, double *bda, double *bdb, 
             double c, double d, int n, int nbdcnd, double *bdc, double *bdd, 
             double elmbda, double *f, int idimf, double& pertrb, int& ierror){


   int lw;
   lw = 4*(n+1) + (13 + int(log(n+1)/log(2.0)))*(m+1);

   double w[lw];
  

   hwscrt_(&a, &b, &m, &mbdcnd, bda, bdb, &c, &d, &n, &nbdcnd, bdc, bdd, &elmbda, f, &idimf, &pertrb, &ierror, w);

 }

//*****************************************************************************************************************   
int my_hstcrt(double a, double b, int m, int mbdcnd, double *bda, double *bdb, 
              double c, double d, int n, int nbdcnd, double *bdc, double *bdd, 
              double elmbda, double *f, int idimf, double& pertrb, int& ierror){


   int lw;
   lw = 13*m + 4*n + m*int(log(n)/log(2.0)); 

   double w[lw];

   hstcrt_(&a, &b, &m, &mbdcnd, bda, bdb, &c, &d, &n, &nbdcnd, bdc, bdd, &elmbda, f, &idimf, &pertrb, &ierror, w);

}

//*****************************************************************************************************************   

   

 
