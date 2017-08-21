// This file contains the C++ interface of the fishpack routine directly callable from my biofilm code
// Main change is that the scalar arguments are passed by value instead of by pointer, but notice that the vector
// arguments are still passed by pointer (would be better if C++ container can be used)

#ifndef _FISH_SUBROUTINE_H
#define _FISH_SUBROUTINE_H

#include "fishpack_db_cpp.h"

int my_hwscrt(double a, double b, int m, int mbdcnd, double *bda, double *bdb, 
              double c, double d, int n, int nbdcnd, double *bdc, double *bdd, 
              double elmbda, double *f, int idimf, double& pertrb, int& ierror);


int my_hstcrt(double a, double b, int m, int mbdcnd, double *bda, double *bdb, 
              double c, double d, int n, int nbdcnd, double *bdc, double *bdd, 
              double elmbda, double *f, int idimf, double& pertrb, int& ierror);


#endif
