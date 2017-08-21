#include "MatVec.h"
#include "LinearSolver.h"
#include <iomanip>
#include <fstream>

//************************************************************************************************************************

UpHessMat::UpHessMat(unsigned s, unsigned c) : size(s), curr(c) {

  v.resize(size+1);

  v[0].resize(size);


  for(int i = 1; i < size+1; i++)
    {
      v[i].resize(size + 1 - i);
    }
}

//************************************************************************************************************************

void UpHessMat::display() const {

  unsigned i, j;

  for(i = 0; i < size + 1; i++)
    {
      if(i == 0)
	{
	  for(j = 0; j < size; j++)
	    cout << v[i][j] << "   ";
	}
      else
        {
          for(j = 0; j < i-1; j++)
            cout << 0 << "  ";
          for(j = i-1; j < size; j++)
            cout << v[i][j-i+1] << "  ";
	}         
      cout << endl;
    }
}

//************************************************************************************************************************

void UpHessMat::LS_solve(valarray<double>& b, valarray<double>& y){

  if(b.size() != curr+1 || y.size() != curr)
    {
      cout << endl << "Dimension incompatible in the Least Square solver! "
           << " Error in UpHessMat::LS_solve, MatVec.cc !! " << endl;
      //XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX DEBUG XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      cout << endl << " b.size() = " << b.size() << ", curr + 1 = " << curr + 1 << ", and  y.size() = " << y.size() << ", curr = " << curr << endl;
      //XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX DEBUG XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      exit(-1);
    }

  int j, k;

  double c, s, r, tmp1, tmp2; 

  //do Givens for each column
  for(j = 0; j < curr; j++)
    {
      //get the Givens transform
      Givens(get(j,j), get(j+1,j), c, s, r);
      //update the diagonal element
      get(j,j) = r; 

      //apply Givens to row j and j+1
      for(k = j+1; k < curr; k++)
	{
	  tmp1 = get(j,k);
          tmp2 = get(j+1,k);

          get(j,k) = c*tmp1 - s*tmp2;
          get(j+1,k) = s*tmp1 + c*tmp2;
        }

      //apply Givens to vector b[j] and b[j+1]
      tmp1 = b[j];
      tmp2 = b[j+1];

      b[j] = c*tmp1 - s*tmp2;
      b[j+1] = s*tmp1 + c*tmp2;
    }

  //Now the matrix is upper triangular, solve for y by backward substitution
  y[curr-1] = b[curr-1]/get(curr-1,curr-1);
  for(j = curr-2; j >= 0; j--)
    {
      for(k = j+1; k < curr; k++)
        b[j] -= get(j,k)*y[k];

      y[j] = b[j]/get(j,j);
    }
}

//************************************************************************************************************************
//***************************class PhiMat*******************************************

PhiMat::PhiMat(size_t nx, size_t ny, int flag) : N_x(nx), N_y(ny), Flag(flag) {

  if(Flag == 0)
    {
      iL = 0;
      iR = N_x;
      jB = 0;
      jT = N_y;
    }
  else if(Flag == 1)
    {
      iL = 0; 
      iR = N_x + 1;
      jB = 0;
      jT = N_y;
    }
  else if(Flag == 2)
    {
      iL = 0;
      iR = N_x + 1;
      jB = 0;
      jT = N_y - 1;
    }
  else
    {
      cout << endl << "The flag for PhiMat is " << Flag << ", it must be 0, 1 or 2 !! Error at line " << __LINE__ << " of file " << __FILE__ << endl;
      exit(-1);
    }

  Nv = 13*(jT - jB +1); //number of vectors containing nonzero elements of the matrix
  v.resize(Nv); // v contains 13*N_y  valarray<double>

  LX = iR - iL + 1; //number of unknow grid points in x-direction
  
  for(int i = 0; i <  Nv; i++)
    v[i].resize(LX);
}     

//************************************************************************************************************************

void PhiMat::Multiply(const valarray<double>& x, valarray<double>& y) const{

  if(x.size() != (jT - jB +1)*(iR - iL +1) || x.size() != y.size())
    {
      cout << endl << " Dimension incompatible in PhiMat::Multiply(), error at line " << __LINE__ << " of file " << __FILE__ << " !! " 
           << endl;
      exit(-1);
    }

  int i, j;  

  int center, curr, shift;

  //treat j = 0, 1, jT - jB - 1, jT - jB differently, for each j, treat i = 0, 1, iR - iL - 1, iR - iL differently

  j = 0; 
 
  center = j*13 + 6; 
  shift  = j*LX;

  i = 0;
  curr = shift + i;
  y[curr] = v[center][i]*x[curr] + v[center+1][i]*x[curr+1] +
    v[center+2][i]*x[curr+2] +  v[center+4][i]*x[curr+LX] + v[center+5][i]*x[curr+1+LX] + 
    v[center+6][i]*x[curr+2*LX];

  i = 1;
  curr = shift + i;

  y[curr] =   v[center-1][i]*x[curr-1] + 
    v[center][i]*x[curr] + v[center+1][i]*x[curr+1] + v[center+2][i]*x[curr+2] + 
    v[center+3][i]*x[curr-1+LX] + v[center+4][i]*x[curr+LX] + v[center+5][i]*x[curr+1+LX] + 
    v[center+6][i]*x[curr+2*LX];

  for(i = 2; i <= iR - iL - 2; i++)
    {
      curr  =  shift + i;
  
      y[curr] =  v[center-2][i]*x[curr-2] + v[center-1][i]*x[curr-1] + 
	v[center][i]*x[curr] + v[center+1][i]*x[curr+1] + v[center+2][i]*x[curr+2] + 
	v[center+3][i]*x[curr-1+LX] + v[center+4][i]*x[curr+LX] + v[center+5][i]*x[curr+1+LX] + 
	v[center+6][i]*x[curr+2*LX];
    }

  i = iR - iL  - 1;
  curr  =  shift + i;

  y[curr] = v[center-2][i]*x[curr-2] + v[center-1][i]*x[curr-1] + 
    v[center][i]*x[curr] + v[center+1][i]*x[curr+1] + 
    v[center+3][i]*x[curr-1+LX] + v[center+4][i]*x[curr+LX] + v[center+5][i]*x[curr+1+LX] + 
    v[center+6][i]*x[curr+2*LX];
 
  i = iR - iL;
  curr  =  shift + i;

  y[curr] =  v[center-2][i]*x[curr-2] + v[center-1][i]*x[curr-1] + v[center][i]*x[curr] +
      v[center+3][i]*x[curr-1+LX] + v[center+4][i]*x[curr+LX] + v[center+6][i]*x[curr+2*LX];

  j = 1; 
 
  center = j*13 + 6; 
  shift =  j*LX;

  i = 0;
  curr = shift + i;
  y[curr] = v[center-4][i]*x[curr-LX]+ 
    v[center-3][i]*x[curr+1-LX] + v[center][i]*x[curr] + v[center+1][i]*x[curr+1] +
    v[center+2][i]*x[curr+2] +  v[center+4][i]*x[curr+LX] + v[center+5][i]*x[curr+1+LX] + 
    v[center+6][i]*x[curr+2*LX];

  i = 1;
  curr = shift + i;

  y[curr] =  v[center-5][i]*x[curr-1-LX] + v[center-4][i]*x[curr-LX]+ 
    v[center-3][i]*x[curr+1-LX] + v[center-1][i]*x[curr-1] + 
    v[center][i]*x[curr] + v[center+1][i]*x[curr+1] + v[center+2][i]*x[curr+2] + 
    v[center+3][i]*x[curr-1+LX] + v[center+4][i]*x[curr+LX] + v[center+5][i]*x[curr+1+LX] + 
    v[center+6][i]*x[curr+2*LX];

  for(i = 2; i <= iR - iL - 2; i++)
    {
      curr  =  shift + i;
  
      y[curr] = v[center-5][i]*x[curr-1-LX] + v[center-4][i]*x[curr-LX]+ 
	v[center-3][i]*x[curr+1-LX] + v[center-2][i]*x[curr-2] + v[center-1][i]*x[curr-1] + 
	v[center][i]*x[curr] + v[center+1][i]*x[curr+1] + v[center+2][i]*x[curr+2] + 
	v[center+3][i]*x[curr-1+LX] + v[center+4][i]*x[curr+LX] + v[center+5][i]*x[curr+1+LX] + 
	v[center+6][i]*x[curr+2*LX];
    }

  i = iR - iL - 1;
  curr  =  shift + i;

  y[curr] = v[center-5][i]*x[curr-1-LX] + v[center-4][i]*x[curr-LX]+ 
    v[center-3][i]*x[curr+1-LX] + v[center-2][i]*x[curr-2] + v[center-1][i]*x[curr-1] + 
    v[center][i]*x[curr] + v[center+1][i]*x[curr+1] + 
    v[center+3][i]*x[curr-1+LX] + v[center+4][i]*x[curr+LX] + v[center+5][i]*x[curr+1+LX] + 
    v[center+6][i]*x[curr+2*LX];
 
  i = iR - iL;
  curr  =  shift + i;

  y[curr] = v[center-5][i]*x[curr-1-LX] + v[center-4][i]*x[curr-LX]+ 
    v[center-2][i]*x[curr-2] + v[center-1][i]*x[curr-1] + v[center][i]*x[curr] +
    v[center+3][i]*x[curr-1+LX] + v[center+4][i]*x[curr+LX] + v[center+6][i]*x[curr+2*LX];


  for(j = 2; j <= jT - jB - 2; j++)
    {
      center = j*13 + 6; 
      shift  = j*LX;

      i = 0;
      curr = shift + i;
      y[curr] = v[center-6][i]*x[curr-2*LX] + v[center-4][i]*x[curr-LX]+ 
         	v[center-3][i]*x[curr+1-LX] + v[center][i]*x[curr] + v[center+1][i]*x[curr+1] +
                v[center+2][i]*x[curr+2] +  v[center+4][i]*x[curr+LX] + v[center+5][i]*x[curr+1+LX] + 
        	v[center+6][i]*x[curr+2*LX];

      i = 1;
      curr = shift + i;

      y[curr] = v[center-6][i]*x[curr-2*LX] + v[center-5][i]*x[curr-1-LX] + v[center-4][i]*x[curr-LX]+ 
	v[center-3][i]*x[curr+1-LX] + v[center-1][i]*x[curr-1] + 
	v[center][i]*x[curr] + v[center+1][i]*x[curr+1] + v[center+2][i]*x[curr+2] + 
	v[center+3][i]*x[curr-1+LX] + v[center+4][i]*x[curr+LX] + v[center+5][i]*x[curr+1+LX] + 
	v[center+6][i]*x[curr+2*LX];

      for(i = 2; i <= iR - iL - 2; i++)
	{
	  curr  =  shift + i;
  
	  y[curr] = v[center-6][i]*x[curr-2*LX] + v[center-5][i]*x[curr-1-LX] + v[center-4][i]*x[curr-LX]+ 
            v[center-3][i]*x[curr+1-LX] + v[center-2][i]*x[curr-2] + v[center-1][i]*x[curr-1] + 
            v[center][i]*x[curr] + v[center+1][i]*x[curr+1] + v[center+2][i]*x[curr+2] + 
            v[center+3][i]*x[curr-1+LX] + v[center+4][i]*x[curr+LX] + v[center+5][i]*x[curr+1+LX] + 
            v[center+6][i]*x[curr+2*LX];
	}

      i = iR - iL - 1;
      curr  =  shift + i;

      y[curr] = v[center-6][i]*x[curr-2*LX] + v[center-5][i]*x[curr-1-LX] + v[center-4][i]*x[curr-LX]+ 
	v[center-3][i]*x[curr+1-LX] + v[center-2][i]*x[curr-2] + v[center-1][i]*x[curr-1] + 
	v[center][i]*x[curr] + v[center+1][i]*x[curr+1] + 
	v[center+3][i]*x[curr-1+LX] + v[center+4][i]*x[curr+LX] + v[center+5][i]*x[curr+1+LX] + 
	v[center+6][i]*x[curr+2*LX];
 
      i = iR - iL;
      curr  =  shift + i;

      y[curr] = v[center-6][i]*x[curr-2*LX] + v[center-5][i]*x[curr-1-LX] + v[center-4][i]*x[curr-LX]+ 
	  v[center-2][i]*x[curr-2] + v[center-1][i]*x[curr-1] + v[center][i]*x[curr] +
  	  v[center+3][i]*x[curr-1+LX] + v[center+4][i]*x[curr+LX] + v[center+6][i]*x[curr+2*LX];
    }

  j = jT - jB - 1; 
 
  center = j*13 + 6; 
  shift =  j*LX;

  i = 0;
  curr = shift + i;
  y[curr] = v[center-6][i]*x[curr-2*LX] + v[center-4][i]*x[curr-LX]+ 
    v[center-3][i]*x[curr+1-LX] + v[center][i]*x[curr] + v[center+1][i]*x[curr+1] +
    v[center+2][i]*x[curr+2] +  v[center+4][i]*x[curr+LX] + v[center+5][i]*x[curr+1+LX];

  i = 1;
  curr = shift + i;

  y[curr] = v[center-6][i]*x[curr-2*LX] + v[center-5][i]*x[curr-1-LX] + v[center-4][i]*x[curr-LX]+ 
    v[center-3][i]*x[curr+1-LX] + v[center-1][i]*x[curr-1] + 
    v[center][i]*x[curr] + v[center+1][i]*x[curr+1] + v[center+2][i]*x[curr+2] + 
    v[center+3][i]*x[curr-1+LX] + v[center+4][i]*x[curr+LX] + v[center+5][i]*x[curr+1+LX];

  for(i = 2; i <= iR - iL - 2; i++)
    {
      curr  =  shift + i;
  
      y[curr] = v[center-6][i]*x[curr-2*LX] + v[center-5][i]*x[curr-1-LX] + v[center-4][i]*x[curr-LX]+ 
	v[center-3][i]*x[curr+1-LX] + v[center-2][i]*x[curr-2] + v[center-1][i]*x[curr-1] + 
	v[center][i]*x[curr] + v[center+1][i]*x[curr+1] + v[center+2][i]*x[curr+2] + 
	v[center+3][i]*x[curr-1+LX] + v[center+4][i]*x[curr+LX] + v[center+5][i]*x[curr+1+LX];
    }

  i = iR - iL - 1;
  curr  =  shift + i;

  y[curr] = v[center-6][i]*x[curr-2*LX] + v[center-5][i]*x[curr-1-LX] + v[center-4][i]*x[curr-LX]+ 
    v[center-3][i]*x[curr+1-LX] + v[center-2][i]*x[curr-2] + v[center-1][i]*x[curr-1] + 
    v[center][i]*x[curr] + v[center+1][i]*x[curr+1] + 
    v[center+3][i]*x[curr-1+LX] + v[center+4][i]*x[curr+LX] + v[center+5][i]*x[curr+1+LX]; 
 
  i = iR - iL;
  curr  =  shift + i;

  y[curr] = v[center-6][i]*x[curr-2*LX] + v[center-5][i]*x[curr-1-LX] + v[center-4][i]*x[curr-LX]+ 
    v[center-2][i]*x[curr-2] + v[center-1][i]*x[curr-1] + v[center][i]*x[curr] +
    v[center+3][i]*x[curr-1+LX] + v[center+4][i]*x[curr+LX];

  j = jT - jB; 
 
  center = j*13 + 6; 
  shift =  j*LX;

  i = 0;
  curr = shift + i;
  y[curr] = v[center-6][i]*x[curr-2*LX] + v[center-4][i]*x[curr-LX]+ 
    v[center-3][i]*x[curr+1-LX] + v[center][i]*x[curr] + v[center+1][i]*x[curr+1] +
    v[center+2][i]*x[curr+2];

  i = 1;
  curr = shift + i;

  y[curr] = v[center-6][i]*x[curr-2*LX] + v[center-5][i]*x[curr-1-LX] + v[center-4][i]*x[curr-LX]+ 
    v[center-3][i]*x[curr+1-LX] + v[center-1][i]*x[curr-1] + 
    v[center][i]*x[curr] + v[center+1][i]*x[curr+1] + v[center+2][i]*x[curr+2];

  for(i = 2; i <= iR - iL - 2; i++)
    {
      curr  =  shift + i;
  
      y[curr] = v[center-6][i]*x[curr-2*LX] + v[center-5][i]*x[curr-1-LX] + v[center-4][i]*x[curr-LX]+ 
	v[center-3][i]*x[curr+1-LX] + v[center-2][i]*x[curr-2] + v[center-1][i]*x[curr-1] + 
	v[center][i]*x[curr] + v[center+1][i]*x[curr+1] + v[center+2][i]*x[curr+2]; 
    }

  i = iR - iL - 1;
  curr  =  shift + i;

  y[curr] = v[center-6][i]*x[curr-2*LX] + v[center-5][i]*x[curr-1-LX] + v[center-4][i]*x[curr-LX]+ 
    v[center-3][i]*x[curr+1-LX] + v[center-2][i]*x[curr-2] + v[center-1][i]*x[curr-1] + 
    v[center][i]*x[curr] + v[center+1][i]*x[curr+1];
 
  i = iR - iL;
  curr  =  shift + i;

  y[curr] = v[center-6][i]*x[curr-2*LX] + v[center-5][i]*x[curr-1-LX] + v[center-4][i]*x[curr-LX]+ 
    v[center-2][i]*x[curr-2] + v[center-1][i]*x[curr-1] + v[center][i]*x[curr];
}

//************************************************************************************************************************

double PhiMat::get_max_element() const {

  int i;

  double tmin, tmax, result;

  result = 0.0;

  for(i = 0; i < Nv; i++)
    {
      tmin = fabs(v[i].min());
      tmax = fabs(v[i].max());
 
      tmax = (tmax > tmin ? tmax : tmin);

      if(tmax > result)
        result = tmax;
    }

  return result;
}

//************************************************************************************************************************

void PhiMat::Scale() {

  double maxele = get_max_element();

  if(maxele <= 1e-10)
    {
      cout << endl << "The matrix is zero, can not be scaled !!!!!!!! " << endl;
      exit(-1);
    }
  if(maxele <= 1.0)
    {
      cout << endl << " The maximum entry of PhiMat matrix is less than 1, no need to scale " << endl;
      return;
    }

  int i;

  for(i = 0; i < Nv; i++)
    v[i] /= maxele;
}

//************************************************************************************************************************

void PhiMat::get_diag(valarray<double>& d) const {

  if(d.size() != get_dim())
    {
      cout << endl << "Dimension incompatible in PhiMat::get_diag(), "
           << " ERROR in MatVec.cc() !!!! " << endl;
      exit(-1);
    }

  int i, j, center;

  for(j = 0; j <= (jT - jB); j++)
    {
      center = j*13 + 6;

      for(i = 0; i <= (iR - iL); i++)
	d[j*LX + i] = v[center][i];
    }
}

//************************************************************************************************************************

void PhiMat::Write_to_file(char* filename) const {

  ofstream ofile(filename,ios::out);

  if(!ofile)
    {
      cout << endl << " Can not open output file " << filename << ", ERROR in PhiMat::Write_to_file " 
           << " at MatVec.cpp !!! " << endl << endl;
      exit(-1);
    }

  int i, j;

  for(i = 0; i < Nv; i++)
    {
      for(j = 0; j < v[i].size() - 1; j++)
	ofile << v[i][j] << " ";
      ofile << v[i][v[i].size()-1] << endl;
    }

  ofile.close();
}

//************************************************************************************************************************

void PhiMat::Read_from_file(char* filename) {

  ifstream ifile(filename,ios::in);

  if(!ifile)
    {
      cout << endl << " Can not open input file " << filename << ", ERROR in PhiMat::Read_from_file " 
           << " at MatVec.cpp !!! " << endl << endl;
      exit(-1);
    }

  int i, j;

  for(i = 0; i < Nv; i++)
    {
      for(j = 0; j < v[i].size(); j++)
	ifile >> v[i][j];
    }

  ifile.close();
}


//************************************************************************************************************************

//class c_Mat 

c_Mat::c_Mat(size_t nx, size_t ny, int flag) : N_x(nx), N_y(ny) , Flag(flag) {
  
  if(Flag == 0)
    {
      iL = 1;
      iR = N_x;
      jB = 0;
      jT = N_y;
    }
  else if(Flag == 1)
    {
      iL = 0; 
      iR = N_x + 1;
      jB = 0;
      jT = N_y;
    }
  else if(Flag == 2)
    {
      iL = 0;
      iR = N_x + 1;
      jB = 0;
      jT = N_y - 1;
    }
  else
    {
      cout << endl << "The flag for c_Mat is " << Flag << ", it must be 0, 1 or 2 !! Error at line " << __LINE__ << " of file " << __FILE__ << endl;
      exit(-1);
    }

  Nv = 5*(jT - jB + 1);
  LX = iR - iL + 1;

  v.resize(Nv); // v contains 5*(N_y + 1)  valarray<double>
  
  for(int i = 0; i < Nv; i++)
    v[i].resize(LX);
}     

//************************************************************************************************************************

void c_Mat::Multiply(const valarray<double>& x, valarray<double>& y) const {

  if(x.size() != get_dim() || x.size() != y.size() )
    {
      cout << endl << "Dimension Incompatible in c_Mat::Multiply(), ERROR in MatVec.cc!!!" << endl;
      exit(-1);
    }

  int i, j;
  int center, curr, shift;

  //at the bottom boundary
  j = 0;

  center = j*5 + 2;
  shift  = j*LX;

  i = 0;
  curr = shift + i;

  y[curr] = v[center][i]*x[curr] + v[center+1][i]*x[curr+1] + 
            v[center+2][i]*x[curr+LX];

  for(i = 1; i <= iR - iL - 1; i++)
    {
      curr = shift + i;

      y[curr] =	v[center-1][i]*x[curr-1] + v[center][i]*x[curr] + v[center+1][i]*x[curr+1] + 
          	v[center+2][i]*x[curr+LX];
    }

  i = iR - iL;
  curr = shift + i;

  y[curr] = v[center-1][i]*x[curr-1] + v[center][i]*x[curr] + 
            v[center+2][i]*x[curr+LX];

  //points between boundary
  for(j = 1; j <= jT - jB - 1; j++)
    {
      center = j*5 + 2;
      shift  = j*LX;

      i = 0;
      curr = shift + i;

      y[curr] = v[center-2][i]*x[curr-LX] + 
                v[center][i]*x[curr] + v[center+1][i]*x[curr+1] + 
           	v[center+2][i]*x[curr+LX];

      for(i = 1; i <= iR - iL - 1; i++)
	{
	  curr = shift + i;

	  y[curr] = v[center-2][i]*x[curr-LX] + 
                    v[center-1][i]*x[curr-1] + v[center][i]*x[curr] + v[center+1][i]*x[curr+1] + 
	            v[center+2][i]*x[curr+LX];
	}

      i = iR - iL;
      curr = shift + i;

      y[curr] = v[center-2][i]*x[curr-LX] + 
        	v[center-1][i]*x[curr-1] + v[center][i]*x[curr] + 
         	v[center+2][i]*x[curr+LX];
    }

  //at the top boundary
  j = jT - jB;

  center = j*5 + 2;
  shift  = j*LX;

  i = 0;
  curr = shift + i;

  y[curr] = v[center-2][i]*x[curr-LX] + 
    v[center][i]*x[curr] + v[center+1][i]*x[curr+1];

  for(i = 1; i <= iR - iL - 1; i++)
    {
      curr = shift + i;

      y[curr] = v[center-2][i]*x[curr-LX] + 
	v[center-1][i]*x[curr-1] + v[center][i]*x[curr] + v[center+1][i]*x[curr+1];
    }

  i = iR - iL;
  curr = shift + i;

  y[curr] = v[center-2][i]*x[curr-LX] + 
    v[center-1][i]*x[curr-1] + v[center][i]*x[curr];

}  

//************************************************************************************************************************

double c_Mat::get_max_element() const {

  int i;

  double tmin, tmax, result;

  result = 0.0;

  for(i = 0; i < Nv; i++)
    {
      tmin = fabs(v[i].min());
      tmax = fabs(v[i].max());
 
      tmax = (tmax > tmin ? tmax : tmin);

      if(tmax > result)
        result = tmax;
    }

  return result;
}

//************************************************************************************************************************

void c_Mat::Scale() {

  double maxele = get_max_element();

  if(maxele == 0.0)
    {
      cout << endl << "The matrix is zero, can not be scaled !!!!!!!! " << endl;
      exit(-1);
    }

  int i;

  for(i = 0; i < Nv; i++)
    v[i] /= maxele;
}

//************************************************************************************************************************

void c_Mat::get_diag(valarray<double>& d) const {

  if(d.size() != get_dim())
    {
      cout << endl << "Dimension incompatible in c_Mat::get_diag(), ERROR in MatVec.cc() !!!! " << endl;
      exit(-1);
    }

  int i, j, center;

  for(j = 0; j <= jT - jB; j++)
    {
      center = j*5 + 2;

      for(i = 0; i <= iR - iL; i++)
	d[j*LX + i] = v[center][i];
    }
}

//************************************************************************************************************************

void c_Mat::Write_to_file(char* filename) const {

  ofstream ofile(filename,ios::out);

  if(!ofile)
    {
      cout << endl << " Can not open output file " << filename << ", ERROR in c_Mat::Write_to_file " 
           << " at MatVec.cpp !!! " << endl << endl;
      exit(-1);
    }

  int i, j;

  for(i = 0; i < Nv; i++)
    {
      for(j = 0; j < v[i].size() - 1; j++)
	ofile << v[i][j] << " ";
      ofile << v[i][v[i].size()-1] << endl;
    }

  ofile.close();
}

//************************************************************************************************************************

void c_Mat::Read_from_file(char* filename) {

  ifstream ifile(filename,ios::in);

  if(!ifile)
    {
      cout << endl << " Can not open input file " << filename << ", ERROR in c_InOutFlow_x_N_y_mat::Read_from_file " 
           << " at MatVec.cpp !!! " << endl << endl;
      exit(-1);
    }

  int i, j;

  for(i = 0; i < Nv; i++)
    {
      for(j = 0; j < v[i].size(); j++)
	ifile >> v[i][j];
    }

  ifile.close();
}


//************************************************************************************************************************

//************************************************************************************************************************

//class c_NoFlux_x_NoFlux_y_mat 

c_NoFlux_x_NoFlux_y_mat::c_NoFlux_x_NoFlux_y_mat(size_t nx, size_t ny) : N_x(nx), N_y(ny) {
  
  int Nv = 5*(N_y+1);
  v.resize(Nv); // v contains 5*(N_y + 1)  valarray<double>
  
  for(int i = 0; i <  Nv; i++)
    v[i].resize(N_x+2);
}     

//************************************************************************************************************************

void  c_NoFlux_x_NoFlux_y_mat::Multiply(const valarray<double>& x, valarray<double>& y) const {

  if( x.size() != get_dim() || x.size() != y.size() )
    {
      cout << endl << "Dimension Incompatible in c_NoFlux_x_NoFlux_y_mat::Multiply(), ERROR in MatVec.cc!!!" << endl;
      exit(-1);
    }

  int i, j;
  int center, curr, shift, LX;

  LX = N_x + 2;
 
  //at the bottom boundary
  j = 0;

  center = j*5 + 2;
  shift  = j*(N_x+2);

  i = 0;
  curr = shift + i;

  y[curr] = v[center][i]*x[curr] + v[center+1][i]*x[curr+1] + 
            v[center+2][i]*x[curr+LX];

  for(i = 1; i <= N_x; i++)
    {
      curr = shift + i;

      y[curr] =	v[center-1][i]*x[curr-1] + v[center][i]*x[curr] + v[center+1][i]*x[curr+1] + 
          	v[center+2][i]*x[curr+LX];
    }

  i = N_x + 1;
  curr = shift + i;

  y[curr] = v[center-1][i]*x[curr-1] + v[center][i]*x[curr] + 
            v[center+2][i]*x[curr+LX];

  //points between boundary
  for(j = 1; j <= N_y - 1; j++)
    {
      center = j*5 + 2;
      shift  = j*(N_x+2);

      i = 0;
      curr = shift + i;

      y[curr] = v[center-2][i]*x[curr-LX] + 
                v[center][i]*x[curr] + v[center+1][i]*x[curr+1] + 
           	v[center+2][i]*x[curr+LX];

      for(i = 1; i <= N_x; i++)
	{
	  curr = shift + i;

	  y[curr] = v[center-2][i]*x[curr-LX] + 
                    v[center-1][i]*x[curr-1] + v[center][i]*x[curr] + v[center+1][i]*x[curr+1] + 
	            v[center+2][i]*x[curr+LX];
	}

      i = N_x + 1;
      curr = shift + i;

      y[curr] = v[center-2][i]*x[curr-LX] + 
        	v[center-1][i]*x[curr-1] + v[center][i]*x[curr] + 
         	v[center+2][i]*x[curr+LX];
    }

  //at the top boundary
  j = N_y;

  center = j*5 + 2;
  shift  = j*(N_x+2);

  i = 0;
  curr = shift + i;

  y[curr] = v[center-2][i]*x[curr-LX] + 
    v[center][i]*x[curr] + v[center+1][i]*x[curr+1];

  for(i = 1; i <= N_x; i++)
    {
      curr = shift + i;

      y[curr] = v[center-2][i]*x[curr-LX] + 
	v[center-1][i]*x[curr-1] + v[center][i]*x[curr] + v[center+1][i]*x[curr+1];
    }

  i = N_x + 1;
  curr = shift + i;

  y[curr] = v[center-2][i]*x[curr-LX] + 
    v[center-1][i]*x[curr-1] + v[center][i]*x[curr];

}  

//************************************************************************************************************************

double c_NoFlux_x_NoFlux_y_mat::get_max_element() const {

  int i, Nv;
 
  Nv = 5*(N_y+1);

  double tmin, tmax, result;

  result = 0.0;

  for(i = 0; i < Nv; i++)
    {
      tmin = fabs(v[i].min());
      tmax = fabs(v[i].max());
 
      tmax = (tmax > tmin ? tmax : tmin);

      if(tmax > result)
        result = tmax;
    }

  return result;
}

//************************************************************************************************************************

void c_NoFlux_x_NoFlux_y_mat::Scale() {

  double maxele = get_max_element();

  if(maxele == 0.0)
    {
      cout << endl << "The matrix is zero, can not be scaled !!!!!!!! " << endl;
      exit(-1);
    }

  int i, Nv;

  Nv = 5*(N_y+1);

  for(i = 0; i < Nv; i++)
    v[i] /= maxele;
}

//************************************************************************************************************************

void c_NoFlux_x_NoFlux_y_mat::get_diag(valarray<double>& d) const {

  if(d.size() != get_dim())
    {
      cout << endl << "Dimension incompatible in c_NoFlux_x_NoFlux_y_mat::get_diag(), ERROR in MatVec.cc() !!!! " << endl;
      exit(-1);
    }

  int i, j, center;

  for(j = 0; j <= N_y; j++)
    {
      center = j*5 + 2;

      for(i = 0; i <= N_x+1; i++)
	d[j*(N_x + 2) + i] = v[center][i];
    }
}

//************************************************************************************************************************

void c_NoFlux_x_NoFlux_y_mat::Write_to_file(char* filename) const {

  ofstream ofile(filename,ios::out);

  if(!ofile)
    {
      cout << endl << " Can not open output file " << filename << ", ERROR in c_NoFlux_x_NoFlux_y_mat::Write_to_file " 
           << " at MatVec.cpp !!! " << endl << endl;
      exit(-1);
    }

  int Nv = 5*(N_y+1);
  int i, j;

  for(i = 0; i < Nv; i++)
    {
      for(j = 0; j < v[i].size() - 1; j++)
	ofile << v[i][j] << " ";
      ofile << v[i][v[i].size()-1] << endl;
    }

  ofile.close();
}

//************************************************************************************************************************

void c_NoFlux_x_NoFlux_y_mat::Read_from_file(char* filename) {

  ifstream ifile(filename,ios::in);

  if(!ifile)
    {
      cout << endl << " Can not open input file " << filename << ", ERROR in c_NoFlux_x_NoFlux_y_mat::Read_from_file " 
           << " at MatVec.cpp !!! " << endl << endl;
      exit(-1);
    }

  int Nv = 5*(N_y+1);
  int i, j;

  for(i = 0; i < Nv; i++)
    {
      for(j = 0; j < v[i].size(); j++)
	ifile >> v[i][j];
    }

  ifile.close();
}

//************************************************************************************************************************

  
//class MultipleRHS***********************************************************************

MultipleRHS::MultipleRHS(int coln, int d) : col_num(coln), dim(d) {
 
  v.resize(col_num);
  for(int i = 0; i < col_num; i++)
    v[i].resize(dim);
}

//************************************************************************************************************************

MultipleRHS::MultipleRHS(const MultipleRHS& M){

  dim  = M.get_dim();
  col_num = M.get_col_num();

  v.resize(col_num);
  for(int i = 0; i < col_num; i++)
    {
      v[i].resize(dim);
      v[i] = M[i];
    }
}
    

//********************************DEBUG******************************
void MultipleRHS::display() {
  cout << endl;
  for(int i = 0; i < col_num; i++)
    {
      cout << "y[" << i << "] = " ;
      for(int j = 0; j < dim; j++)
        cout << setw(8) << v[i][j] << "   ";
      cout << endl;
    }
}

//************************************************************************************************************************

//function to calcute the gradient of a MultipleRHS variable (with size (Nx+2)x(Ny+1),
//the result is a vector with size (Nx+2)*(Ny-1), but we only use its first (Nx+1) columns
void MultipleRHS::To_Gradient(double dx, double dy, valarray<double>& P_x, valarray<double>& P_y, int bc_flag) const {

  int i, j, Nx, Ny, T;

  Nx = dim - 2;
  Ny = col_num - 1;

  if(P_x.size() != P_y.size() || P_x.size() != (Nx + 2)*(Ny - 1)  )
    {
      cout << endl << "Dimension Incompatible in MultipleRHS::To_Gradient(), ERROR in MatVec.cc!!!" << endl;
      exit(-1);
    }

  if(bc_flag == 1)     //periodic BC in x-direction, pay attention at the boundary!!!
    {
      T = Nx + 1;

      for(j = 1; j <= Ny-1; j++)
	{
	  {
	    i = 0;
	    P_x[(j-1)*(Nx+2)+i] = .5*(v[j][1] - v[j][Nx])/dx; 
	    P_y[(j-1)*(Nx+2)+i] = .5*(v[j+1][i] - v[j-1][i])/dy;
	  }
       
	  for(i = 1; i <= Nx; i++)
	    {
	      P_x[(j-1)*(Nx+2)+i] = .5*(v[j][i+1] - v[j][i-1])/dx;
	      P_y[(j-1)*(Nx+2)+i] = .5*(v[j+1][i] - v[j-1][i])/dy;
	    }
	  {
	    i = Nx + 1;
	    P_x[(j-1)*(Nx+2)+i] = .5*(v[j][1] - v[j][Nx])/dx;
	    P_y[(j-1)*(Nx+2)+i] = .5*(v[j+1][i] - v[j-1][i])/dy;
	  }
	}
    }

  else if(bc_flag == 0) //No-flux BC (or Dirichelet BC for u, v) in x-direction
    {
      for(j = 1; j <= Ny-1; j++)
	{
	  for(i = 1; i <= Nx; i++)
	    {
	      P_x[(j-1)*(Nx+2)+i] = .5*(v[j][i+1] - v[j][i-1])/dx;
	      P_y[(j-1)*(Nx+2)+i] = .5*(v[j+1][i] - v[j-1][i])/dy;
	    }
         
	  // p_x, p_y at i = 0 and i = Nx + 1 are 0
	  P_x[(j-1)*(Nx+2)+0] = 0.0;
	  P_y[(j-1)*(Nx+2)+0] = 0.0;

	  P_x[(j-1)*(Nx+2)+Nx+1] = 0.0;
	  P_y[(j-1)*(Nx+2)+Nx+1] = 0.0;
	}
    }

//   else if(bc_flag == 2) //out-flow BC for u, v in x-direction
//     {
//       for(j = 1; j <= Ny-1; j++)
// 	{
	   
// 	  i = 0;
// 	  P_x[(j-1)*(Nx+2)+i] = (v[j][i+1] - v[j][i])/dx;
// 	  P_y[(j-1)*(Nx+2)+i] = .5*(v[j+1][i] - v[j-1][i])/dy;
	     

// 	  for(i = 1; i <= Nx; i++) 
// 	    {
// 	      P_x[(j-1)*(Nx+2)+i] = .5*(v[j][i+1] - v[j][i-1])/dx;
// 	      P_y[(j-1)*(Nx+2)+i] = .5*(v[j+1][i] - v[j-1][i])/dy;
// 	    }

// 	  {
// 	    i = Nx + 1;

// 	    P_x[(j-1)*(Nx+2)+i] = (v[j][i] - v[j][i-1])/dx;
// 	    P_y[(j-1)*(Nx+2)+i] = .5*(v[j+1][i] - v[j-1][i])/dy;
// 	  }
// 	}
//     }

  else
    {
      cout << endl << " bc_flag must be 0, 1, ERROR in MultipleRHS::To_Gradient() at MatVec.cc() !!! " << endl ;
      exit(-1);
    }
    
}
//************************************************************************************************************************
	    
//overload operator* which multiply a double variable to a MultipleRHS object 
MultipleRHS operator*(double a, const MultipleRHS& v){

  MultipleRHS w(v);

  int colnum = v.get_col_num();

  for(int i = 0; i < colnum; i++)
    w[i] *= a;

  return w;
}

//************************************************************************************************************************

//class TriDiag

TriDiag::TriDiag(size_t _n, int _flag) : n(_n), flag(_flag) {

  d.resize(n);
  e.resize(n-1);
  f.resize(n-1);
}


//************************************************************************************************************************

TriDiag::TriDiag(const TriDiag& T) {

  n = T.get_n();
  flag = T.get_flag();

  d = T.get_d();
  e = T.get_e();
  f = T.get_f();
}
 

//************************************************************************************************************************
 
//function to solve the tridiagonal system with given right hand side y 
//if flag == 0, non-singular, regular chase method,
//if flag == 1, set x[0] = 0, then x[1] = f[0]/f[0]
//and the modified system can be solved by direct substitution 
//after the solve, y is overriden by the solution
void TriDiag::ChaseSolve(MultipleRHS& y) const {

    int ic, j, col_num;

    if(y.get_dim() != n)
      {
       cout << endl << "Dimension incompatible in Chase solve of the tridiagonal system." 
            << endl << "Error in TriDiag.ChaseSolve(), MatVec.cc !" << endl;
       exit(-1);
      }

    col_num = y.get_col_num();

    if(flag == 0) //non-singular case
      {
       valarray<double> b(n);
       valarray<double> a(n-1);

       b[0] = d[0];
       for(j = 1; j < n; j++)
         {
	   a[j-1] = e[j-1]/b[j-1];
           b[j] = d[j] - a[j-1]*f[j-1];
	 }

       for(j = 1; j < n; j++)
	 for(ic = 0; ic < col_num; ic++)
           y[ic][j] -= a[j-1]*y[ic][j-1];

       for(ic = 0; ic < col_num; ic++)
	 y[ic][n-1] /= b[n-1];
 
       for(j = n-2; j >= 0; j--)
         for(ic = 0; ic < col_num; ic++)
	   y[ic][j] = (y[ic][j] - f[j]*y[ic][j+1])/b[j];
      }

    else if(flag == 1) // the tridiagonal system is singular, set y[ic][0] = 0
      {
        MultipleRHS solu(y);

	for(ic = 0; ic < col_num; ic++)
	  {
	    solu[ic][0] = 0.0;
	    solu[ic][1] = y[ic][0]/f[0];

	    for(j = 2; j <= n-1; j++)
	      solu[ic][j] = (y[ic][j-1] - e[j-2]*solu[ic][j-2] - d[j-1]*solu[ic][j-1])/f[j-1];

// 	    solu[ic][n-1] = (y[ic][n-1] - e[n-2]*solu[ic][n-2])/f[n-1];
	  }

	y = solu;
      }

             
//         valarray<double> b(n-1);
//         valarray<double> a(n-2);

//         b[0] = d[1];
//         for(j = 1; j < n-1; j++)
//           {
// 	    a[j-1] = e[j]/b[j-1];
//             b[j] = d[j+1] - a[j-1]*f[j];
// 	  }

//         for(j = 2; j < n; j++)
//           for(ic = 0; ic < col_num; ic++)
// 	    y[ic][j] -= a[j-2]*y[ic][j-1];

//         for(ic = 0; ic < col_num; ic++)
//           y[ic][n-1] /= b[n-2];

//         for(j = n-2; j >= 1; j--)
//           for(ic = 0; ic < col_num; ic++)
// 	    y[ic][j] = (y[ic][j] - f[j]*y[ic][j+1])/b[j-1];
//       }

    else
      {
	cout << endl << "Wrong flag in TriDiag.ChaseSolve(), flag must be 0 or 1 ! " 
             << endl << "Error in MatVec.cc " << endl;
        exit(-1);
      }
  }


//************************************************************************************************************************

//The function to calculate the inner product of two vector
double InnerProd(const valarray<double>& u, const valarray<double>& v){

  if(u.size() != v.size())
    {
      cout << endl << "Take Inner Product of two vectors with different size!! ERROR in LinearSolver.cc, InnerProd() !! " 
           << endl;
      exit(-1);
    }

  double res = 0.0;

  for(int i = 0; i < u.size(); i++)
    res += u[i]*v[i];

  return res;
}

//************************************************************************************************************************

//function to calculate the smallest absolute value of a vector
double Mini_abs_element(const valarray<double>& v){
  int i;

  double result = fabs(v[0]);

  for(i = 1; i < v.size(); i++)
{
   if(result > fabs(v[i]))
           result = fabs(v[i]);
}

  return result;
}

//********************************************* END of MatVec.cc *********************************************************
