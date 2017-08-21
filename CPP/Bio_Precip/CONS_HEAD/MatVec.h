//This file defines the data structure for the program, i.e., how to represent the 
//matrix and vector in the program
#ifndef _MATVEC_H
#define _MATVEC_H

#include <iostream>
#include <valarray>

using namespace std;

//************************************************************************************************************************

//this class represents an upper Hessenberg matrix with size (k+1)*k by
//k+1 vectors as its rows, with varying size of each row

class UpHessMat {

  valarray< valarray<double> > v;

  unsigned size; // value of k

  unsigned curr; // size will be used in the least square solve 

 public:

  UpHessMat(unsigned s = 0, unsigned c=0); //constructor

  //destructor
  ~UpHessMat() {}

  //the function to get the (i,j) entry of the matrix, here j is the natural column index
  double& get(int i, int j) {
    if(i > size || j > size - 1 || j < i - 1)
      {
  	cout << endl << "Out of range error in get() of class UpHessMat, in MatVec.h"
             << endl << "i = " << i << ", j = " << j << ", size = " << size << ", size -1 = " << size - 1
             << ", i -1 = " << i - 1 << endl; 
        exit(-1);
      }
    if(i == 0)
      return  v[i][j];
    else 
      return  v[i][j-i+1];
  }
  //the function to get the (i,j) entry of the matrix, here j is the natural column index
  //return const reference
  const double& get(int i, int j)const {
    if(i > size || j > size - 1 || j < i - 1)
      {
	cout << endl << "Out of range error in get() const of class UpHessMat, in MatVec.h" 
             << endl << "i = " << i << ", j = " << j << ", size = " << size << ", size -1 = " << size - 1
             << ", i -1 = " << i - 1 << endl; 

        exit(-1);
      }
    if(i == 0)
      return  v[i][j];
    else 
      return  v[i][j-i+1];
  }
  // function to display the matrix on screen, for debug purpose
  void display() const;

  //function to set the value of curr
  void set_curr(unsigned c) {
    curr = c;
    }

  unsigned get_curr() const {
    return curr;
  }

  //function to get the size of the matrix
  unsigned get_size() const {
    return size;
  }

  //Given curr and a vector b of size curr+1, this function solve for the solution y of the least square 
  //problem min || b - v*y||, note here after the solve, v is changed , so it can be used only once
  void LS_solve(valarray<double>& b, valarray<double>& y);
  
};

//************************************************************************************************************************

//The abstract base class describing a matrix vector multimplication
class MatMulVec{

 public:

  MatMulVec(){};

  virtual ~MatMulVec() {}

  virtual size_t get_dim() const = 0;

  virtual void Multiply(const valarray<double>& x, valarray<double>& y) const = 0;

  virtual double get_max_element() const = 0; //function to get the maximum magnitude of elements

  virtual void Scale() = 0;        //divide each element by its maximum magnitude of elements

  virtual void get_diag(valarray<double>&) const =0; //get the vector containg the diagonal of the matrix

  virtual void Write_to_file(char* filename) const = 0;   //write the matrix into a file

  virtual void Read_from_file(char* ) = 0;   //Read the matrix from a file
};

//************************************************************************************************************************

//the class represents the matrix obtained by discretizating the phi-transportation equation
//Boundary condition is indicated by Flag, it determines the x, y index range for the unknows
//Doesn NOT include the periodic case
//flag == 0, N_D_x_N_y (out flow bc)
//flag == 1, N_x_N_y (no-flux)
//flag == 2, N_x_N_D_y (no-flux at x and y = 0, Dirichlet at y = 1)
//iL, iR, jB, jT are the index ranged determined by Flag
//the matrix is represented by 13*Ny vectors
class PhiMat : public MatMulVec{

  size_t N_x, N_y, Nv;  

  int Flag, iL, iR, jB, jT, LX;     //flag identifies the type of BCs, index range for unknowns, depends on flag (BC)

  valarray< valarray<double> > v;

 public:

  PhiMat(size_t nx=0, size_t ny = 0, int flag = 0);

  ~PhiMat () { };

  size_t get_dim() const {

    return (iR - iL + 1)*(jT - jB + 1);
  }

  size_t get_Nx() const {
    return N_x;
  }

  size_t get_Ny() const {
    return N_y;
  }

  int get_LX() const {
    return LX;
  }

  void Multiply(const valarray<double>& x, valarray<double>& y) const;

  double get_max_element() const; //function to get the maximum magnitude of elements

  void Scale();        //divide each element by its maximum magnitude of elements

  void get_diag(valarray<double>&) const; //get the vector containg the diagonal of the matrix

  //overload operator[]
  valarray<double>& operator[](int i) {
    return v[i];
  }

  //overload operator[]
  const valarray<double> operator[](int i) const {
    return v[i];
  }

  //write the matrix into a file
  void Write_to_file(char* filename) const;

  //Read the matrix from a file
  void Read_from_file(char* );
};    


//************************************************************************************************************************

//the class to represent matrix obtained by discretizating the transportation equation 
//Boundary condition is indicated by Flag, it determines the x, y index range for the unknows
//Doesn NOT include the periodic case
//flag == 0, D_x_N_y (in-out flow bc)
//flag == 1, N_x_N_y (no-flux)
//flag == 2, N_x_N_D_y (no-flux at x and y = 0, Dirichlet at y = 1)
//iL, iR, jB, jT are the index ranged determined by Flag
//for the substrate, with in-out flow BC in x direction
class c_Mat : public MatMulVec{

  size_t N_x, N_y, Nv;

  int Flag, iL, iR, jB, jT, LX;     //flag identifies the type of BCs, index range for unknowns, depends on flag (BC)

  valarray< valarray<double> > v;

 public:

  c_Mat(size_t nx=0, size_t ny = 0, int flag = 0);

  ~c_Mat() { };

  size_t get_dim() const {

    return LX*(jT - jB + 1);
  }

  size_t get_Nx() const {
    return N_x;
  }

  size_t get_Ny() const {
    return N_y;
  }

  int get_LX() const {
    return LX;
  }

  void Multiply(const valarray<double>& x, valarray<double>& y) const;

  double get_max_element() const; //function to get the maximum magnitude of elements

  void Scale();        //divide each element by its maximum magnitude of elements

  void get_diag(valarray<double>&) const; //get the vector containg the diagonal of the matrix

  //overload operator[]
  valarray<double>& operator[](int i) {
    return v[i];
  }

  //overload operator[]
  const valarray<double> operator[](int i) const {
    return v[i];
  }

  //write the matrix into a file
  void Write_to_file(char* filename) const;

  //Read the matrix from a file
  void Read_from_file(char* );
};

//************************************************************************************************************************

//the class to represent matrix obtained by discretizating the transportation equation 
//for the substrate, with no-flux BC in both x and y direction
class c_NoFlux_x_NoFlux_y_mat : public MatMulVec{

  size_t N_x, N_y;
  valarray< valarray<double> > v;

 public:

  c_NoFlux_x_NoFlux_y_mat(size_t nx=0, size_t ny = 0);

  ~c_NoFlux_x_NoFlux_y_mat() { };

  size_t get_dim() const {

    return (N_x+2)*(N_y + 1);
  }

  size_t get_Nx() const {
    return N_x;
  }

  size_t get_Ny() const {
    return N_y;
  }


  void Multiply(const valarray<double>& x, valarray<double>& y) const;

  double get_max_element() const; //function to get the maximum magnitude of elements

  void Scale();        //divide each element by its maximum magnitude of elements

  void get_diag(valarray<double>&) const; //get the vector containg the diagonal of the matrix

  //overload operator[]
  valarray<double>& operator[](int i) {
    return v[i];
  }

  //overload operator[]
  const valarray<double> operator[](int i) const {
    return v[i];
  }

  //write the matrix into a file
  void Write_to_file(char* filename) const;

  //Read the matrix from a file
  void Read_from_file(char* );
};

//************************************************************************************************************************

//the class representing the RHS of a linear system, it can have multiple columns
//for our problem, dim is Nx+2 and col_num is Ny + 1
class MultipleRHS {

  int dim;      //dimension of the RHS vector
  int col_num;  //number of column on RHS

  valarray< valarray<double> > v;

 public:

  MultipleRHS(int coln, int d);
  MultipleRHS(const MultipleRHS& M); //copy constructor

  ~MultipleRHS() { }

  int get_dim() const {

    return dim;
  }

  int get_col_num() const {
 
    return col_num;
  }

  //overload operator[]
  valarray<double>& operator[](int i) {
    return v[i];
  }

  //overload operator[]
  valarray<double> operator[](int i) const {
    return v[i];
  }

  //function to calcute the gradient of a MultipleRHS variable (with size (Nx+2)x(Ny+1),
  //the result is a vector with size (Nx+2)*(Ny-1), P_x and P_y are derivative w.r.t. x and y respectively
  //bc_flag = 1 : periodic BC in x, = 0 : Dirichelt for u,v in x
  void To_Gradient(double dx, double dy, valarray<double>& P_x, valarray<double>& P_y, int bc_flag) const;

  //*************************DEBUG**************************** 
  void display();
  //**********************************************************
};

//************************************************************************************************************************

//overload operator* which multiply a double variable to a MultipleRHS object 
MultipleRHS operator*(double a, const MultipleRHS& v);

//************************************************************************************************************************

//the class representing a tridiagonal matrix as 3 vectors
// n is the number of unknowns
// d is the main diagonal of the coefficient matrix
// e is the first subdiagonal of the coefficient matrix
// f is the first superdiagonal of the coefficient matrix
// x is the right hand side of the equations, x can have multiple
//   column, and is represented by another class
// We also consdier a possible singular matrix arise from Neumann BC from the
// Poisson equation, we use the flag to indicate the case: 
// flag == 0, non-singular, flag ==1, singular

class TriDiag {

  size_t n;
  valarray<double> d;
  valarray<double> e;
  valarray<double> f;
  int flag;

 public:

  TriDiag(size_t _n, int _flag);

  ~TriDiag() {};

  TriDiag(const TriDiag& T);

  size_t get_n() const {
    return n;
  }

  int get_flag() const {
    return flag;
  }

  valarray<double>& get_d() {
    return d;
  }

  const valarray<double> get_d() const {
    return d;
  }

  valarray<double>& get_e() {
    return e;
  }

  const valarray<double> get_e() const {
    return e;
  }

  valarray<double>& get_f() {
    return f;
  } 

  const valarray<double> get_f() const {
    return f;
  } 


  //function to solve the tridiagonal system with given right hand side y 
  //if flag == 0, non-singular, regular chase method,
  //if flag == 1, set x[0] = 1, and modify y[1] -= e[0]*1, and solve the
  //resulting system with order n-1
  //after the solve, y is overriden by the solution
  void ChaseSolve(MultipleRHS& y) const ;

};

//************************************************************************************************************************

//The function to calculate the inner product of two vector
double InnerProd(const valarray<double>& u, const valarray<double>& v);

//************************************************************************************************************************

//function to calculate the smallest absolute value of a vector
double Mini_abs_element(const valarray<double>& v);


#endif
