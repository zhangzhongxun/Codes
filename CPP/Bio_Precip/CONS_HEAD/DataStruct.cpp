#include "DataStruct.h"

//************************************************************************************************************************
//the global parameters are defined here

double Gamma_0;   //surface of the CaCO3
double Gamma_1;   //used in the chemican potential
double Gamma_2;   //used in the chemican potential
double Kai_CH;  //used in the chemican potential
double Lambda_c;    //mobility factor of CaCO3
double Lambda_b;    //mobility factor of biofilm
double Lambda_s;    //mobility factor of solvent
double sub_epsilon; //the small parameter control the growth rate of the network
double sub_mu;    //used in definition of g_b
double sub_K_c;   //used in definition of g_b
double sub_A;     //used in definition of g_b
double Ds;        //the diffusion coefficient for the nutrient for biofilm
double D_ur;      //the diffusion coefficient for urea
double D_nh3;     //the diffusion coefficient for [NH3]
double D_nh4;     //the diffusion coefficient for [NH4+]
double D_h2co3;   //the diffusion coefficient for [H2CO3]
double D_hco3;    //the diffusion coefficient for [HCO3 -]
double D_co3;     //the diffusion coefficient for [CO3 2-]
double D_ca;      //the diffusion coefficient for [Ca2+]
double D_cl;      //the diffusion coefficient for [Cl-]
double D_h;       //the diffusion coefficient for [H+]
double D_oh;      //the diffusion coefficient for [OH-]
double D_cahco3;  //the diffusion coefficient for [CaHCO3+]
double D_caco3;  //the diffusion coefficient for aqueous [CaCO3]
const double Z_ur = 0.0, Z_nh3 = 0.0, Z_nh4 = 1.0, Z_h2co3 = 0.0, Z_hco3 = -1.0, Z_co3 = -2.0, Z_ca = 2.0, 
                           Z_cl = -1.0, Z_h = 1.0, Z_oh = -1.0, Z_electron = 1.602e-19, Z_c = 0.0, Z_cahco3 = 1.0,  Z_caco3 = 0.0; //charge of variaous particles
double eta_cac;   //viscosity of CaCO3  (should be very large, virtually infinity
double eta_bio;   //viscosity of biofilm
double eta_sol;   //viscosity of solvent
double rho_cac;   //density of CaCO3 (variable density of fluid?)
double rho_bio;   //density of biofilm
double rho_sol;   //density of solvent
double NP_inv;    //the 1/(N_n*V_n/V_s), used in the mixed energy
double PHI_reg;   //the parameter to regularize log(phi) at phi = 0, i.e., use log(phi+PHI_reg)
double PHIS_reg;  //the parameter to regularize the transportation equation when phis is zero
double eta_ave;   //average viscosity, used in the LHS of the discretized NSE
double k_ur;      //the rate for urea hydrolysis
double k_p;       //the rate of CaCO3 precipitation
double S_crit;    //the critial value of saturation state for CaCO3 to begin precipitate
double K_acid_1;  //acidity constant 1
double K_acid_2;  //acidity constant 2
double K_nh;      //equilibrium rate constant for NH3 hydrolysis
double K_w;       //Equilibrium constant for water
double K_cahco3;  // Equilibrium constant for [CaHCO3-] <-> [Ca2+] + [HCO3-], K4 in the function
double K_caco3;  // Equilibrium constant for [CaCO3] <-> [Ca2+] + [CO3^2-], K5 in the function
double K_so;      //equilibrium calcite solubility product
double A_Davies, B_Davies, alpha_Davies;  //constants used in Davies equation for activity coefficient
double Ca_phase;  //the phase separation parameter for CaCO3, 1/epsilon, should be big
double Ca_vol_coef; //coefficient of volume change from [Ca2+] to CaCO3
double Vol_Penal; //the penalty cosntant for phi_c + phi_b + phi_s = 1
char   SaveDir[100];//Director to save data
int chemc_flag;   //0: double wells (4-th order polynomial), 1: quadratic
double a_EPU;     //the unit constant for electric potential
double L_theta;   //the ratio of computation domain length to the whole domain length (< 1)
double Q_flux;    //the flux from the left boundary, induced by given pressure drop
double U_ref;     //the maximum velocity at the centerline of the left boundary induced by given pressure drop
double G_pres;    //the linear pressure gradient induced by the 
double Re_s;      //Reynolds of the solvent
int Newton_max_k; // maximum iteratoin number for Newton's method
double Newton_tol; //Tolerance for Newton's method
int Activity_max_k; //maximum iteration number for activity coefficients update
double Activity_tol; //Tolerance for successive update of activity coefficients

//************************************************************************************************************************

Node_Data::Node_Data(size_t nx, size_t ny) : N_x(nx), N_y(ny) {

  v.resize(N_x + 2);
  
  for(int i = 0; i < N_x + 2; i++)
    v[i].resize(N_y + 1);
}

Node_Data::Node_Data(const Node_Data& u){

  N_x = u.get_Nx();
  N_y = u.get_Ny();

  v.resize(N_x + 2);
  
  for(int i = 0; i < N_x + 2; i++)
    {
      v[i].resize(N_y+1);
      v[i] = u[i];
    }
}


//************************************************************************************************************************

Node_Data& Node_Data::operator=(const Node_Data& u){

  int y_differ = 0;

  if(this != &u) //avoid self-assignment
    {
      if(N_x != u.get_Nx())
	{
	  N_x = u.get_Nx();
	  v.resize(N_x + 2);
	}
  
      if(N_y != u.get_Ny())
	{
          y_differ = 1;
	  N_y = u.get_Ny();
	}

      for(int i = 0; i < N_x + 2; i++)
	{
	  if(y_differ == 1)
	    {
	      v[i].resize(N_y+1);
	    }
	  v[i] = u[i];
	}
    }

  return *this;
}


//************************************************************************************************************************
Node_Data&  Node_Data::operator=(double a){


    for(int i = 0; i < N_x + 2; i++)
	{
	  v[i] = a;
	}

    return *this;
}

//************************************************************************************************************************

Node_Data& Node_Data::operator+=(const Node_Data& u){

  if(u.get_Nx() != N_x || u.get_Ny() != N_y)
    {
      cout << endl << "Dimension Incompatible in Node_Data::operator+=, ERROR in DataStruct.cc !!!!!!" << endl;
      exit(-1);
    }

  for(int i = 0; i <= N_x + 1; i++)
    v[i] += u[i];

  return *this;
}

//************************************************************************************************************************

double Node_Data::L_inf_norm() const {
 
  double norm = 0.0;

  int i, j;

  for(j = 0; j <= N_y; j++)
    for(i = 0; i <= N_x + 1; i++)
      {
        if(fabs(v[i][j]) > norm)
	   norm = fabs(v[i][j]);
      }

  return norm;
}

//************************************************************************************************************************

double  Node_Data::sum() const {

  double result = 0.0;

  for(int i = 0; i < N_x + 2; i++)
    result += v[i].sum();

  return result;
}


//************************************************************************************************************************

double  Node_Data::max() const {

  double result = v[0].max();

  for(int i = 1; i < N_x + 2; i++)
    if(v[i].max() > result)
      result = v[i].max();

  return result;
}

//************************************************************************************************************************

double  Node_Data::min() const {

  double result = v[0].min();

  for(int i = 1; i < N_x + 2; i++)
    if(v[i].min() < result)
      result = v[i].min();

  return result;
}

//************************************************************************************************************************

Node_Data_Ext_N_x_N_y::Node_Data_Ext_N_x_N_y(size_t nx, size_t ny) : N_x(nx), N_y(ny) {

  v.resize(N_x + 6);
  
  for(int i = 0; i < v.size(); i++)
    v[i].resize(N_y + 5);
} 

Node_Data_Ext_N_x_N_y::Node_Data_Ext_N_x_N_y(const Node_Data& u){

  int i, j;

  N_x = u.get_Nx();
  N_y = u.get_Ny();

  v.resize(N_x + 6);
  
  for(i = 0; i < v.size(); i++)
    {
      v[i].resize(N_y+5);
    }

  for(i = 2; i <= N_x + 3; i++)
    for(j = 2; j <= N_y + 2; j++)
      v[i][j] = u[i-2][j-2];

  //now extend it outside the boundary
  for(i = 2; i <= N_x + 3; i++)
    {
      j = 0;
      v[i][j] = u[i-2][2];
     
      j = 1; 
      v[i][j] = u[i-2][1];

      j = N_y + 3;
      v[i][j] = u[i-2][N_y-1];

      j = N_y + 4;
      v[i][j] = u[i-2][N_y-2];
    }

  for( j = 2; j <= N_y + 2; j++)
    {
      i = 0; 
      v[i][j] = u[2][j-2];

      i = 1;
      v[i][j] = u[1][j-2];
 
      i = N_x + 4;
      v[i][j] = u[N_x][j-2];
     
      i = N_x + 5;
      v[i][j] = u[N_x-1][j-2];
    }

  v[1][1] = u[1][1];
  v[0][1] = u[2][1];
  v[1][0] = u[1][2];
  v[0][0] = u[2][2];

  v[N_x+4][1] = u[N_x][1];
  v[N_x+4][0] = u[N_x][2];
  v[N_x+5][1] = u[N_x-1][1];
  v[N_x+5][0] = u[N_x-1][2];

  v[1][N_y+3] = u[1][N_y-1];
  v[1][N_y+4] = u[1][N_y-2];
  v[0][N_y+3] = u[2][N_y-1];
  v[0][N_y+4] = u[2][N_y-2];

  v[N_x+4][N_y+3] = u[N_x][N_y-1];
  v[N_x+4][N_y+4] = u[N_x][N_y-2];
  v[N_x+5][N_y+3] = u[N_x-1][N_y-1];
  v[N_x+5][N_y+4] = u[N_x-1][N_y-2];
}


//************************************************************************************************************************

Node_Data_Ext_IO_x_N_y::Node_Data_Ext_IO_x_N_y(size_t nx, size_t ny) : N_x(nx), N_y(ny) {

  v.resize(N_x + 4);
  
  for(int i = 0; i < v.size(); i++)
    v[i].resize(N_y + 3);
}

Node_Data_Ext_IO_x_N_y::Node_Data_Ext_IO_x_N_y(const Node_Data& u, const valarray<double>& uR){

  int i, j;

  N_x = u.get_Nx();
  N_y = u.get_Ny();

  v.resize(N_x + 4);
  
  for(i = 0; i < v.size(); i++)
    {
      v[i].resize(N_y+3);
    }

  for(i = 1; i <= N_x+2; i++)
    for(j = 1; j <= N_y+1; j++)
      v[i][j] = u[i-1][j-1];

  //no-flux in y-direction
  for(i = 1; i <= N_x+2; i++)
    {
      j = 0;
      v[i][j] = u[i-1][1];

      j = N_y+2;
      v[i][j] = u[i-1][N_y-1];
    }

  for(j = 1; j <= N_y+1; j++)
    {
      //In-flow at left boundary
      i = 0;
      v[i][j] = u[0][j-1];

      //Out-flow at right boundary
      i = N_x + 3;
      v[i][j] = uR[j-1];
    }

  //take care the corner, though may not be used at all
  v[0][0] = v[1][0];
  v[0][N_y+2] = v[1][N_y+2];
  v[N_x+3][0] = v[N_x+2][0];
  v[N_x+3][N_y+2] = v[N_x+2][N_y+2];        
}

//************************************************************************************************************************

Node_Data_Ext_P_x_N_y::Node_Data_Ext_P_x_N_y(size_t nx, size_t ny) : N_x(nx), N_y(ny) {

  v.resize(N_x + 5);
  
  for(int i = 0; i < v.size(); i++)
    v[i].resize(N_y + 5);
}

Node_Data_Ext_P_x_N_y::Node_Data_Ext_P_x_N_y(const Node_Data& u){

  int i, j;

  N_x = u.get_Nx();
  N_y = u.get_Ny();

  v.resize(N_x + 5);
  
  for(i = 0; i < v.size(); i++)
    {
      v[i].resize(N_y+5);
    }

  for(i = 2; i <= N_x + 3; i++)
    for(j = 2; j <= N_y + 2; j++)
      v[i][j] = u[i-2][j-2];

  //now extend it outside the boundary
  for(i = 2; i <= N_x + 3; i++)
    {
      j = 0;
      v[i][j] = u[i-2][2];
     
      j = 1; 
      v[i][j] = u[i-2][1];

      j = N_y + 3;
      v[i][j] = u[i-2][N_y-1];

      j = N_y + 4;
      v[i][j] = u[i-2][N_y-2];
    }

  for( j = 2; j <= N_y + 2; j++)
    {
      i = 0; 
      v[i][j] = u[N_x-1][j-2];

      i = 1;
      v[i][j] = u[N_x][j-2];
 
      i = N_x + 4;
      v[i][j] = u[1][j-2];
    }

  v[1][1] = v[1][3]; 
  v[0][1] = v[0][3];
  v[1][0] = v[1][4];
  v[0][0] = v[0][4];

  v[N_x+4][1] = v[N_x+4][3];
  v[N_x+4][0] = v[N_x+4][4];

  v[1][N_y+3] = v[1][N_y+1];
  v[1][N_y+4] = v[1][N_y];
  v[0][N_y+3] = v[0][N_y+1];
  v[0][N_y+4] = v[0][N_y];

  v[N_x+4][N_y+3] = v[N_x+4][N_y+1];
  v[N_x+4][N_y+4] = v[N_x+4][N_y];
}

//************************************************************************************************************************

//overload operator+ for Node_Data
Node_Data operator+(const Node_Data&u, const Node_Data& v){

  if( u.get_Nx() != v.get_Nx() || u.get_Ny() != v.get_Ny() )
    {
     cout << endl << "Dimension Incompatible in operator+(Node_Data,Node_Data), ERROR in DataStruct.cc !" << endl;
     exit(-1);
    }

  Node_Data w(u);

  int i, Nx;

  Nx = u.get_Nx();
  for(i = 0; i <= Nx + 1; i++)
    w[i] += v[i];

  return w;
}

//************************************************************************************************************************

//overload operator+ which add a Node_Data to a vector with size (Nx+2)*(Ny+1)
Node_Data operator+(const Node_Data& u, const valarray<double>& v){

  int i, j, Nx, Ny;

  Nx = u.get_Nx();
  Ny = u.get_Ny();

  if(v.size() != (Nx+2)*(Ny+1))
    {
     cout << endl << "Dimension Incompatible in operator+(Node_Data,valarray<double>), ERROR in DataStruct.cc !" << endl;
     exit(-1);
    }
 
  Node_Data w(u);

  for(j = 0; j <= Ny; j++)
    for(i = 0; i <= Nx+1; i++)
      w[i][j] += v[j*(Nx+2)+i];

  return w;
}

//************************************************************************************************************************

//overload operator+ which add a Node_Data to a MultipleRHS object with dim = Nx+2 and col_num = Ny+1
Node_Data operator+(const Node_Data& u, const MultipleRHS& v){

  int i, j, Nx, Ny;

  Nx = u.get_Nx();
  Ny = u.get_Ny();

  if(v.get_dim() != (Nx+2) || v.get_col_num() !=(Ny+1))
    {
     cout << endl << "Dimension Incompatible in operator+(Node_Data,MultipleRHS), ERROR in DataStruct.cc !" << endl;
     exit(-1);
    }
 
  Node_Data w(u);

  for(j = 0; j <= Ny; j++)
    for(i = 0; i <= Nx+1; i++)
      w[i][j] += v[j][i];

  return w;
}


//************************************************************************************************************************
//overload operator+ for Node_Data, which calculate double + Node_Data
Node_Data operator+(double a, const Node_Data& u){

  Node_Data w(u);

  int i, Nx;

  Nx = u.get_Nx();
  for(i = 0; i <= Nx + 1; i++)
    w[i] = a + w[i];

  return w;
}

//************************************************************************************************************************
//overload operator+ for Node_Data, which calculate Node_Data + double
Node_Data operator+(const Node_Data& u, double a){

  return a+u;
}

//************************************************************************************************************************
//overload operator- for Node_Data, which calculate double - Node_Data
Node_Data operator-(double a, const Node_Data& u){

  Node_Data w(u);

  int i, Nx;

  Nx = u.get_Nx();
  for(i = 0; i <= Nx + 1; i++)
    w[i] = a - w[i];

  return w;
}

//************************************************************************************************************************

//overload operator* for Node_Data and double variable
Node_Data operator*(double a, const Node_Data& u){
 
  Node_Data w(u);

  int i, Nx;

  Nx = u.get_Nx();
  for(i = 0; i <= Nx + 1; i++)
    w[i] *= a;

  return w;
}


//************************************************************************************************************************

//overload operator/ for double variable and Node_Data
Node_Data operator/(double a, const Node_Data& u){
 
  Node_Data w(u.get_Nx(), u.get_Ny());

  int i, Nx;

  Nx = u.get_Nx();
  for(i = 0; i <= Nx + 1; i++)
    w[i] = a/u[i];

  return w;
}
//************************************************************************************************************************
void Add_Array(double* u, double* v, int n){

  for(int i = 0; i < n; i++)
    u[i] += v[i];
}

//************************************************************************************************************************

//function to point-wise multiply two Node_Data, the result is stored in the SECOND argument
void PointWise_Multiply_Node_Data(const Node_Data& u, Node_Data& v){

  int i, j, Nx, Ny;

  if(u.get_Nx() != v.get_Nx() || u.get_Ny() != v.get_Ny())
    {
      cout << endl << "Dimension Incompatible PointWise_Multiply_Node_Data(), ERROR in DataStruct.cc !!! "<< endl;
      exit(-1);
    }

  Nx = u.get_Nx();
  Ny = u.get_Ny();

  for(j = 0; j <= Ny; j++)
    for(i = 0; i <= Nx + 1; i++)
      v[i][j] *= u[i][j];
}

//************************************************************************************************************************

//function to point-wise multiply a Node_Data and a MultipleRHS, the result is stored in the SECOND argument
//which is the MultipleRHS 
void PointWise_Multiply_Node_Data_MultipleRHS(const Node_Data& u, MultipleRHS& v){

  int i, j, Nx, Ny;

  if(u.get_Nx() != v.get_dim() - 2 || u.get_Ny() != v.get_col_num() - 1)
    {
      cout << endl << "Dimension Incompatible PointWise_Multiply_Node_Data_MultipleRHS(), "
           << " ERROR in DataStruct.cc !!! "<< endl;
      exit(-1);
    }

  Nx = u.get_Nx();
  Ny = u.get_Ny();

  for(j = 0; j <= Ny; j++)
    for(i = 0; i <= Nx + 1; i++)
      v[j][i] *= u[i][j];
}

//************************************************************************************************************************
void Node_Data_Array_Conversion(Node_Data& u, double* v, int flag){

  int i, j, Nx, Ny;

  Nx = u.get_Nx();
  Ny = u.get_Ny();

  if(flag == 1) // Node_Data to Array
    {
      for(j = 0; j <= Ny; j++)
	for(i = 0; i <= Nx+1; i++)
          v[j*(Nx+2)+i] = u[i][j];
    }
  else if(flag == 2) //Array to Node_Data
    {
      for(j = 0; j <= Ny; j++)
	for(i = 0; i <= Nx+1; i++)
	  u[i][j] = v[j*(Nx+2)+i];
    }
  else
    {
      cout << endl << "flag in Node_Data_Array_Conversion() must be 1 or 2, ERROR in DataStruct.cpp !! " << endl;
      exit(-1);
    }
}

//************************************************************************************************************************
void MultipleRHS_Vector_Conversion(valarray<double>& v, MultipleRHS& mrhs, int flag){

  int i, j, Nx, Ny;

  Nx = mrhs.get_dim() - 2;
  Ny = mrhs.get_col_num() - 1;

  if(v.size() != (Nx+2)*(Ny+1) - 1)
    {
           cout << endl << "Dimension Incompatible MultipleRHS_Vector_Conversion(), "
           << " ERROR in DataStruct.cc !!! "<< endl;
      exit(-1);
    }

  if(flag == 1) //mrhs to v, remove the first element of mrhs
    {
      {
      j = 0;
 
      for(i = 1; i <= Nx+1; i++)
        v[i-1] = mrhs[j][i];
      }

      for(j = 1; j <= Ny; j++)
        for(i = 0; i <= Nx+1; i++)
          v[j*(Nx+2)+i-1] = mrhs[j][i];
    }

  else if(flag == 2) // v to mrhs, set the first element of mrhs to ZERO or the average of v[0][1] and v[1][0]
    {
      mrhs[0][0] = .5*(v[0] + v[Nx+1]); // 0.0;

      {
      j = 0;
  
      for(i = 1; i <= Nx + 1; i++)
     	mrhs[j][i] = v[i-1];
      }

      for(j = 1; j <= Ny; j++)
	for(i = 0; i <= Nx+1; i++)
          mrhs[j][i] = v[j*(Nx+2)+i-1];
    }

  else
    {
      cout << endl << "The flag in  MultipleRHS_Vector_Conversion() must be 1 or 2, ERROR in DataStruct.cc !! " << endl;
      exit(-1);
    }
      
}

//************************************************************************************************************************
void MultipleRHS_Array_Conversion(double* v, MultipleRHS& mrhs, int flag){


  int i, j, Nx, Ny;

  Nx = mrhs.get_dim() - 2;
  Ny = mrhs.get_col_num() - 1;

  //have to make sure the dimension of v and mrhs are compatible

  if(flag == 1) //mrhs to v
    {
      for(j = 0; j <= Ny; j++)
        for(i = 0; i <= Nx+1; i++)
          v[j*(Nx+2)+i] = mrhs[j][i];
    }

  else if(flag == 2) // v to mrhs
    {
      for(j = 0; j <= Ny; j++)
        for(i = 0; i <= Nx+1; i++)
	  mrhs[j][i] = v[j*(Nx+2)+i];
    }

  else
    {
      cout << endl << "The flag in  MultipleRHS_Vector_Conversion() must be 1 or 2, ERROR in DataStruct.cc !! " << endl;
      exit(-1);
    }
}

//************************************************************************************************************************
//we have to make sure the dimension of p, px, py are compatible
void STG_pressure_gradient(int Nx, int Ny, double dx, double dy, double* p, valarray<double>& P_x, valarray<double>& P_y, int bc_flag){

  int i, j;

  for(j = 1; j <= Ny-1; j++)
    for(i = 1; i <= Nx; i++)
      {
	P_x[(j-1)*(Nx+2)+i] = .5*((p[(j-1)*(Nx+1)+i] + p[j*(Nx+1)+i]) - (p[(j-1)*(Nx+1)+i-1] + p[j*(Nx+1)+i-1]))/dx;
	P_y[(j-1)*(Nx+2)+i] = .5*((p[j*(Nx+1)+i-1] + p[j*(Nx+1)+i]) - (p[(j-1)*(Nx+1)+i-1] + p[(j-1)*(Nx+1)+i]))/dy;
      }

  if(bc_flag == 1)     //periodic BC in x-direction, pay attention at the boundary!!!
    {
      for(j = 1; j <= Ny-1; j++)
	{
	  {
	    i = 0;
	    P_x[(j-1)*(Nx+2)+i] = .5*((p[(j-1)*(Nx+1)+i] + p[j*(Nx+1)+i]) - (p[(j-1)*(Nx+1)+Nx] + p[j*(Nx+1)+Nx]))/dx;
	    P_y[(j-1)*(Nx+2)+i] = .5*((p[j*(Nx+1)+Nx] + p[j*(Nx+1)+i]) - (p[(j-1)*(Nx+1)+Nx] + p[(j-1)*(Nx+1)+i]))/dy;
	  }
 
	  {
	    i = Nx + 1;
	    P_x[(j-1)*(Nx+2)+i] = P_x[(j-1)*(Nx+2)];
	    P_y[(j-1)*(Nx+2)+i] = P_y[(j-1)*(Nx+2)];
	  }
	}
    }

  else if(bc_flag == 0) //No-flux BC (or Dirichelet BC for u, v) in x-direction
    {
      for(j = 1; j <= Ny-1; j++)
	{
	  // p_x, p_y at i = 0 and i = Nx + 1
          i = 0;
	  P_x[(j-1)*(Nx+2)+i] = .5*((p[(j-1)*(Nx+1)+i+1] + p[j*(Nx+1)+i+1]) - (p[(j-1)*(Nx+1)+i] + p[j*(Nx+1)+i]))/dx;
	  P_y[(j-1)*(Nx+2)+i] = .5*((3*p[j*(Nx+1)+i] - p[j*(Nx+1)+i+1]) - (3*p[(j-1)*(Nx+1)+i] - p[(j-1)*(Nx+1)+i+1]))/dy;

          i = Nx + 1;
	  P_x[(j-1)*(Nx+2)+i] = .5*((p[(j-1)*(Nx+1)+i-1] + p[j*(Nx+1)+i-1]) - (p[(j-1)*(Nx+1)+i-2] + p[j*(Nx+1)+i-2]))/dx;
	  P_y[(j-1)*(Nx+2)+i] = .5*((3*p[j*(Nx+1)+i-1] - p[j*(Nx+1)+i-2]) - (3*p[(j-1)*(Nx+1)+i-1] - p[(j-1)*(Nx+1)+i-2]))/dy;
	}
    }

  else
    {
      cout << endl << " bc_flag must be 0, 1, ERROR in MultipleRHS::To_Gradient() at MatVec.cc() !!! " << endl ;
      exit(-1);
    }

}


//************************************************************************************************************************
//overload the output operator << for Node_Data
ostream& operator<<(ostream& _ostr, const Node_Data& u){

  int i, j, Nx, Ny;

  Nx = u.get_Nx();
  Ny = u.get_Ny();

  //std::_Ios_Fmtflags oldflags = _ostr.flags(ios::scientific | ios::showpos | ios::showpoint );
  std::ios::fmtflags oldflags = _ostr.flags(ios::scientific | ios::showpos | ios::showpoint );

  int oldprecision =  _ostr.precision(16);   

  for(j = 0; j <= Ny; j++)
    {
      for(i = 0; i <= Nx + 1; i++)
        _ostr << u[i][j] << "  ";

      _ostr << endl;
    }

  _ostr << endl;

  _ostr.precision(oldprecision);

  _ostr.flags(oldflags);

  return _ostr;

}

//************************************************************************************************************************
 
//overload the input operator >> for Node_Data
istream& operator>>(istream& _istr, Node_Data& u){

  int i, j, Nx, Ny;

  Nx = u.get_Nx();  
  Ny = u.get_Ny();

  for(j = 0; j <= Ny; j++)
    for(i = 0; i <= Nx+1; i++)
      _istr >> u[i][j];

  return _istr;
}

//************************************************************************************************************************

//function to extrapolate the Node_Data at time n+1 by data at time n and n-1
void Extrapolate_Node_Data(const Node_Data& dn, const Node_Data& dn_1, double dt_n, double dt_n_1, double alpha, Node_Data& d) {

if(dn.get_Nx() != dn_1.get_Nx() || dn.get_Nx() != d.get_Nx() || dn.get_Ny() != dn_1.get_Ny() || dn.get_Ny() != d.get_Ny())
  {
   cout << endl << "Dimension Incompatible in Extrapolate_Node_Data(), ERROR in DataStruct.cc !!! " << endl;
   exit(-1);
  }

int i, j, Nx, Ny;

Nx = d.get_Nx();
Ny = d.get_Ny();

for(i = 0; i <= Nx + 1; i++)
  for(j = 0; j <= Ny; j++)
    d[i][j] = dn[i][j] + alpha*dt_n*(dn[i][j] - dn_1[i][j])/dt_n_1;

}

//************************************************************************************************************************

void Interpolate_Node_Data(const Node_Data& dnp1, const Node_Data& dn, double alpha, Node_Data& d) {

if(dnp1.get_Nx() != dn.get_Nx() || dn.get_Nx() != d.get_Nx() || dnp1.get_Ny() != dn.get_Ny() || dn.get_Ny() != d.get_Ny())
  {
   cout << endl << "Dimension Incompatible in Extrapolate_Node_Data(), ERROR in DataStruct.cc !!! " << endl;
   exit(-1);
  }

int i, j, Nx, Ny;

Nx = d.get_Nx();
Ny = d.get_Ny();

for(i = 0; i <= Nx + 1; i++)
  for(j = 0; j <= Ny; j++)
    d[i][j] = alpha*dnp1[i][j] + (1.0 - alpha)*dn[i][j];

}

//************************************************************************************************************************
 
//function to confine every entry of a Node_Data variable within [a,b], mainly used to confine phi and c
void Confine_Node_Data(Node_Data& u, double a, double b){

  int i, j, Nx, Ny;

  Nx = u.get_Nx();
  Ny = u.get_Ny();

  if(a > b)
    {
      cout << endl << " a > b in  Confine_Node_Data() at line " << __LINE__ << " of " << __FILE__ << "!! a mus be <= b !! " << endl;
      exit(-1);
    }

  for(j = 0; j <= Ny; j++)
    for(i = 0; i <= Nx + 1; i++)
      {
        if(u[i][j] < a)
	  u[i][j] = a;
        else if(u[i][j] > b)
          u[i][j] = b;
      }
}

//************************************************************************************************************************

Node_Data_x_Average::Node_Data_x_Average(size_t nx, size_t ny) : N_x(nx), N_y(ny) {
 
  v.resize(N_x + 1);

  for(int i = 0; i < N_x + 1; i++)
    v[i].resize(N_y+1);
}

//************************************************************************************************************************

Node_Data_x_Average::Node_Data_x_Average(const Node_Data& D) {

  N_x = D.get_Nx();
  N_y = D.get_Ny();
 
  v.resize(N_x + 1);

  for(int i = 0; i < N_x + 1; i++)
    v[i].resize(N_y+1);

  int i, j;

  for(i = 0; i <= N_x; i++)
    for(j = 0; j <= N_y; j++)
      v[i][j] = .5*(D[i][j] + D[i+1][j]);
}


//************************************************************************************************************************

Node_Data_x_Ave_Ext_N::Node_Data_x_Ave_Ext_N(size_t nx, size_t ny) : N_x(nx), N_y(ny) {
 
  v.resize(N_x + 3);

  for(int i = 0; i < v.size(); i++)
    v[i].resize(N_y+1);
}

//************************************************************************************************************************

Node_Data_x_Ave_Ext_N::Node_Data_x_Ave_Ext_N(const Node_Data& D) {

  N_x = D.get_Nx();
  N_y = D.get_Ny();
 
  v.resize(N_x + 3);

  for(int i = 0; i < v.size(); i++)
    v[i].resize(N_y+1);

  int i, j;

  for(i = 1; i <= N_x+1; i++)
    for(j = 0; j <= N_y; j++)
      v[i][j] = .5*(D[i-1][j] + D[i][j]);

  for(j = 0; j <= N_y; j++)
    {
      i = 0;
      v[i][j] = v[1][j];

      i = N_x + 2;
      v[i][j] = v[i-1][j];
    }
}


//************************************************************************************************************************

Node_Data_x_Ave_Ext_P::Node_Data_x_Ave_Ext_P(size_t nx, size_t ny) : N_x(nx), N_y(ny) {
 
  v.resize(N_x + 2);

  for(int i = 0; i < v.size(); i++)
    v[i].resize(N_y+1);
}

//************************************************************************************************************************

Node_Data_x_Ave_Ext_P::Node_Data_x_Ave_Ext_P(const Node_Data& D) {

  N_x = D.get_Nx();
  N_y = D.get_Ny();
 
  v.resize(N_x + 2);

  for(int i = 0; i < v.size(); i++)
    v[i].resize(N_y+1);

  int i, j;

  for(i = 1; i <= N_x+1; i++)
    for(j = 0; j <= N_y; j++)
      v[i][j] = .5*(D[i-1][j] + D[i][j]);

  for(j = 0; j <= N_y; j++)
    {
      i = 0;
      v[i][j] = v[N_x+1][j];
    }
}

//************************************************************************************************************************

Node_Data_x_Ave_Ext_IO::Node_Data_x_Ave_Ext_IO(size_t nx, size_t ny) : N_x(nx), N_y(ny) {
 
  v.resize(N_x + 3);

  for(int i = 0; i < v.size(); i++)
    v[i].resize(N_y+1);
}

//************************************************************************************************************************

Node_Data_x_Ave_Ext_IO::Node_Data_x_Ave_Ext_IO(const Node_Data_Ext_IO_x_N_y& D) {

  N_x = D.get_Nx();
  N_y = D.get_Ny();
 
  v.resize(N_x + 3);

  for(int i = 0; i < v.size(); i++)
    v[i].resize(N_y+1);

  int i, j;

  for(i = 0; i <= N_x+2; i++)
    for(j = 0; j <= N_y; j++)
      v[i][j] = .5*(D.get(i-1,j) + D.get(i,j));
}

//************************************************************************************************************************

Node_Data_y_Average::Node_Data_y_Average(size_t nx, size_t ny) : N_x(nx), N_y(ny) {
 
  v.resize(N_x + 2);

  for(int i = 0; i < v.size(); i++)
    v[i].resize(N_y);
}

//************************************************************************************************************************

Node_Data_y_Average::Node_Data_y_Average(const Node_Data& D) {

  N_x = D.get_Nx();
  N_y = D.get_Ny();
 
  v.resize(N_x + 2);

  for(int i = 0; i < v.size(); i++)
    v[i].resize(N_y);

  int i, j;
   
  for(i = 0; i <= N_x+1; i++)
    for(j = 0; j < N_y; j++)
      v[i][j] = .5*(D[i][j] + D[i][j+1]);
}

//************************************************************************************************************************

Node_Data_y_Ave_Ext_N::Node_Data_y_Ave_Ext_N(size_t nx, size_t ny) : N_x(nx), N_y(ny) {
 
  v.resize(N_x + 2);

  for(int i = 0; i < v.size(); i++)
    v[i].resize(N_y + 2);
}

//************************************************************************************************************************

Node_Data_y_Ave_Ext_N::Node_Data_y_Ave_Ext_N(const Node_Data& D) {

  N_x = D.get_Nx();
  N_y = D.get_Ny();
 
  v.resize(N_x + 2);

  for(int i = 0; i < v.size(); i++)
    v[i].resize(N_y+2);

  int i, j;
  
  for(j = 1; j <= N_y; j++) 
    for(i = 0; i <= N_x+1; i++)
       v[i][j] = .5*(D[i][j-1] + D[i][j]);

  for(i = 0; i <= N_x + 1; i++)
    {
      j = 0;
      v[i][j] = v[i][j+1];

      j = N_y + 1;
      v[i][j] = v[i][j-1];
    }

}

//************************************************************************************************************************

Node_Data_y_Ave_Ext_N::Node_Data_y_Ave_Ext_N(const Node_Data_Ext_IO_x_N_y& D){

  N_x = D.get_Nx();
  N_y = D.get_Ny();
 
  v.resize(N_x + 2);

  for(int i = 0; i < v.size(); i++)
    v[i].resize(N_y+2);

  int i, j;

  for(i = 0; i <= N_x+1; i++)
    for(j = 0; j <= N_y+1; j++)
      v[i][j] = .5*(D.get(i,j-1) + D.get(i,j));
} 

//************************************************************************************************************************

Edge_Data::Edge_Data(size_t nx, size_t ny) : N_x(nx), N_y(ny) {

  int i;

  vx.resize(N_x + 1);
  
  for(i = 0; i < N_x + 1; i++)
    vx[i].resize(N_y + 1);

  vy.resize(N_x + 2);
  for(i = 0; i < N_x + 2; i++)
    vy[i].resize(N_y);
}


