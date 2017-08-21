//the file defines the data structure for the model, and how to compute the 
//elements in each datastructure

#ifndef _DATASTRUCT_H
#define _DATASTRUCT_H

#include <valarray>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "MatVec.h"

using namespace std;

//decleration  of global variables
extern double Gamma_0;   //surface of the CaCO3
extern double Gamma_1; //used in the chemican potential
extern double Gamma_2; //used in the chemican potential
extern double Kai_CH;  //used in the chemican potential
extern double Lambda_c;    //mobility factor of CaCO3
extern double Lambda_b;    //mobility factor of biofilm
extern double Lambda_s;    //mobility factor of solvent
extern double sub_epsilon; //the small parameter control the growth rate of the network
extern double sub_mu;  //used in definition of g_b
extern double sub_K_c; //used in definition of g_b
extern double sub_A;   //used in definition of g_b
extern double Ds;      //the diffusion coefficient for the nutrient for biofilm
extern double D_ur;      //the diffusion coefficient for urea
extern double D_nh3;     //the diffusion coefficient for [NH3]
extern double D_nh4;     //the diffusion coefficient for [NH4+]
extern double D_h2co3;   //the diffusion coefficient for [H2CO3]
extern double D_hco3;    //the diffusion coefficient for [HCO3 -]
extern double D_co3;     //the diffusion coefficient for [CO3 2-]
extern double D_ca;      //the diffusion coefficient for [Ca2+]
extern double D_cl;      //the diffusion coefficient for [Cl-]
extern double D_h;       //the diffusion coefficient for [H+]
extern double D_oh;      //the diffusion coefficient for [OH-]
extern double eta_cac;   //viscosity of CaCO3  (should be very large, virtually infty)
extern double eta_bio;   //viscosity of network
extern double eta_sol;   //viscosity of solvent
extern double rho_cac;   //density of CaCO3 (variable density of fluid?)
extern double rho_bio;   //density of network
extern double rho_sol;   //density of solvent
extern double NP_inv;    //the 1/(N_n*V_n/V_s), used in the mixed energy
extern double PHI_reg;   //the parameter to regularize log(phi) at phi = 0, i.e., use log(phi+PHI_reg)
extern double PHIS_reg;  //the parameter to regularize the transportation equation when phis is zero
extern double eta_ave;   //average viscosity, used in the LHS of the discretized NSE
extern double k_ur;      //the rate for urea hydrolysis
extern double k_p;       //the rate of CaCO3 precipitation
extern double S_crit;    //the critial value of saturation state for CaCO3 to begin precipitate
extern double K_acid_1;  //Equilibrium rate constant of [H2CO3] <--> [HCO3-] + [H+],  K1 in the function
extern double K_acid_2;  //Equilibrium rate constant of [HCO3-] <--> [CO3 2-] + [H+], K2 in the function
extern double K_nh;      //Equilibrium rate constant of [NH3] + [H+] <--> [NH4 +], K3 in the function
extern double K_w;       //Equilibrium constant for water = 1e-14
extern double K_so;      //equilibrium calcite solubility product
extern double Ca_phase;  //the phase separation parameter for CaCO3, = 1/epsilon, should be big
extern double Ca_vol_coef; //coefficient of volume change from [Ca2+] to CaCO3
extern double Vol_Penal; //the penalty cosntant for phi_c + phi_b + phi_s = 1
extern char SaveDir[100];//Director to save data
extern int chemc_flag;   //0: double wells (4-th order polynomial), 1: quadratic
extern const double Z_ur, Z_nh3, Z_nh4, Z_h2co3, Z_hco3, Z_co3, Z_ca, Z_cl, Z_h, Z_oh, Z_c, Z_electron;
extern double a_EPU;     //the unit constant for electric potential

//************************************************************************************************************************

//define a class contain all the control parameters such as dx, dy , Nx, k (m = 2^(k+1), Ny = m + 1)
//total time T_final, CFL, tolerance of GMRES, ...
class Control_Parameter {

public:
  
    Control_Parameter() {}   //default constructor do nothing
    ~Control_Parameter() {}  //destructor do nothing

    double dx, dy, T_final, CFL, GMRES_tol, rho_0, c_0, h, t_0, tfac, alph, eps_dif, Shear_v; 

    double Ini_phc, Ini_phb, Ini_Ur, Ini_NH3, Ini_NH4, Ini_H2CO3, Ini_HCO3, Ini_CO3, Ini_Ca, Ini_Cl, Ini_H, Ini_OH,
           Ini_c, cutoff_c_time, D_phic_tol, D_delta, a_epu;

    size_t Nx, k;

    int Init_flag,  BC_flag, Grid_flag,  velo_flag, Init_N, GMRES_max_k, F_Stride, N_dif_steps, D_flag, Difu_Init_flag, Ave_Vis_flag, EP_flag;

    char savedir[100];
    
    void display(char*) const;
};

//************************************************************************************************************************

//define a class contain all the global parameters, these are the values before nondimensionalization
//it has a member function to convert them to non-dimensionalized parameters
class Global_Parameter {

public:
 
    Global_Parameter () {}
    ~Global_Parameter () {}

    double gam_0, gam_1, gam_2, KB, T, kai_CH, lam_c, lam_b, lam_s, s_eps, s_mu, s_K_c, s_A, ds, d_ur, d_nh3, d_nh4, d_h2co3, d_hco3, d_co3, d_ca, d_cl, d_h, d_oh,
           eta_b, eta_s, eta_c, rho_b, rho_s, rho_c, kur, kp, s_crit;
           
    double np_inv, phi_reg, phis_reg, kacid_1, kacid_2, k_nh, k_w, k_so, ca_phase, ca_vol_coef, vol_penal;

    //function to convert parameters to non-dimensionalized value
    void To_Nondimension(const Control_Parameter& CP);
    void display(char*) const;
};



//************************************************************************************************************************
  
//this class is for data defined at grid-node, it has
//size (N_x+2)x(N_y + 1), in x-direction we have N_x +2
//grid since we use the outflow boundary condition

class Node_Data {

  size_t N_x, N_y;

  valarray< valarray<double> > v;

 public:

  Node_Data(size_t nx=0, size_t ny=0);
  //copy constructor
  Node_Data(const Node_Data& u);
 
  ~Node_Data() { }

  size_t get_Nx() const {
    return N_x;
  }

  size_t get_Ny() const {
    return N_y;
  }

  //overload operator=
  Node_Data& operator=(const Node_Data& u);

  //overload operator= (assign constant value a to Node_Data)
  Node_Data& operator=(double a);

  //overload operator+=
  Node_Data& operator+=(const Node_Data& u);

  //overload operator[]
  valarray<double>& operator[](int i){
    return v[i];
  }

  //overload operator[]
  const valarray<double> operator[](int i) const {
    return v[i];
  }

  //function to return the L_infty norm of the variable, for debuging purpose
  double L_inf_norm() const;

  //function to return the sum of all the entries of a Node_Data variable
  double sum() const;

  double max() const;

  double min() const;
};

//************************************************************************************************************************

//this class extends the Node_Data to index i = -2, -1, Nx+2, Nx+3, j = -2, -1, Ny+1, Ny+2, using No-flux BC
//it has size (N_x+6)x(N_y + 5),
  
class Node_Data_Ext_N_x_N_y {

  size_t N_x, N_y;

  valarray< valarray<double> > v;

 public:

  Node_Data_Ext_N_x_N_y(size_t nx=0, size_t ny=0);
  //copy constructor
  Node_Data_Ext_N_x_N_y(const Node_Data& u);
 
  ~Node_Data_Ext_N_x_N_y() { }

  size_t get_Nx() const {
    return N_x;
  }

  size_t get_Ny() const {
    return N_y;
  }

  //overload operator[]
  valarray<double>& operator[](int i){
    return v[i];
  }

  //overload operator[]
  const valarray<double> operator[](int i) const {
    return v[i];
  }

  //function to return entry (i,j) where (i,j) is the index of a Node_Data variable, i.e., 0 <= i <= Nx+1, 0 <= j <= Ny
  double get(int i, int j) const {
    return v[i+2][j+2];
  }

};


//************************************************************************************************************************

//this class extends the Node_Data to index i =, -1, Nx+2, j = -1, Ny+1, using In-Out flow BC in x and No-flux BC in y
//it has size (N_x+4)x(N_y+3),

class Node_Data_Ext_IO_x_N_y {

  size_t N_x, N_y;

  valarray< valarray<double> > v;

 public:

  Node_Data_Ext_IO_x_N_y(size_t nx=0, size_t ny=0);
  //copy constructor
  Node_Data_Ext_IO_x_N_y(const Node_Data& u, const valarray<double>& uR);
 
  ~Node_Data_Ext_IO_x_N_y() { }

  size_t get_Nx() const {
    return N_x;
  }

  size_t get_Ny() const {
    return N_y;
  }

  //overload operator[]
  valarray<double>& operator[](int i){
    return v[i];
  }

  //overload operator[]
  const valarray<double> operator[](int i) const {
    return v[i];
  }

  //function to return entry (i,j) where (i,j) is the index of a Node_Data variable, i.e., 0 <= i <= Nx+1, 0 <= j <= Ny
  double get(int i, int j) const {
    return v[i+1][j+1];
  }

};


//************************************************************************************************************************

//this class extends the Node_Data to index i = -2, -1, Nx+2, j = -2, -1, Ny+1, Ny+2, using No-flux BC in y-direction, periodic BC in x-direction
//it has size (N_x+5)x(N_y + 5),

class Node_Data_Ext_P_x_N_y {

  size_t N_x, N_y;

  valarray< valarray<double> > v;

 public:

  Node_Data_Ext_P_x_N_y(size_t nx=0, size_t ny=0);
  //copy constructor
  Node_Data_Ext_P_x_N_y(const Node_Data& u);
 
  ~Node_Data_Ext_P_x_N_y() { }

  size_t get_Nx() const {
    return N_x;
  }

  size_t get_Ny() const {
    return N_y;
  }

  //overload operator[]
  valarray<double>& operator[](int i){
    return v[i];
  }

  //overload operator[]
  const valarray<double> operator[](int i) const {
    return v[i];
  }

  //function to return entry (i,j) where (i,j) is the index of a Node_Data variable periodic in x, i.e., 0 <= i <= Nx, 0 <= j <= Ny
  double get(int i, int j) const {
    return v[i+2][j+2];
  }

};
//************************************************************************************************************************

//overload the output operator << for Node_Data
ostream& operator<<(ostream& _ostr, const Node_Data& u);

//overload the input operator << for Node_Data
istream& operator>>(istream& _istr, Node_Data& u);

//overload operator+ for Node_Data, which add two Node_Data variable
Node_Data operator+(const Node_Data& u, const Node_Data& v);

//overload operator+ which add a Node_Data to a vector with size (Nx+2)*(Ny+1)
Node_Data operator+(const Node_Data& u, const valarray<double>& v);

//overload operator+ which add a Node_Data to a MultipleRHS object with dim = Nx+2 and col_num = Ny+1
Node_Data operator+(const Node_Data& u, const MultipleRHS& v);

//overload operator+ for Node_Data, which calculate double + Node_Data
Node_Data operator+(double a, const Node_Data& u);

//overload operator+ for Node_Data, which calculate Node_Data + double
Node_Data operator+(const Node_Data& u, double a);

//overload operator- for Node_Data, which calculate double - Node_Data
Node_Data operator-(double a, const Node_Data& u);

//overload operator* for Node_Data and double variable
Node_Data operator*(double a, const Node_Data& u);

//overload operator/ for double variable and Node_Data, need to make sure u has no zero element!!!!
Node_Data operator/(double a, const Node_Data& u);

//************************************************************************************************************************

//function to add two arrays, put the result in first argument
void Add_Array(double* u, double* v, int n);

//************************************************************************************************************************
//function to point-wise multiply two Node_Data, the result is stored in the SECOND argument
void PointWise_Multiply_Node_Data(const Node_Data& u, Node_Data& v);

//************************************************************************************************************************

//function to point-wise multiply a Node_Data and a MultipleRHS, the result is stored in the SECOND argument
//which is the MultipleRHS 
void PointWise_Multiply_Node_Data_MultipleRHS(const Node_Data& u, MultipleRHS& v);

//************************************************************************************************************************

//function to do conversion between a Node_Data and a double* array,
//flag = 1: Node_Data to double*, flag = 2: double* to Node_Data
void Node_Data_Array_Conversion(Node_Data& u, double* v, int flag);

//************************************************************************************************************************

//function to do conversion between a vector v and a MultipleRHS object mrhs, use flag to indicate the direction
//flag = 1: MultipleRHS to vector, but remove the first element
//flag = 2: vector to MultipleRHS, set the first element of MultipleRHS to ZERO
//note v.size = mrhs.get_col()*mrhs.get_dim() - 1
void MultipleRHS_Vector_Conversion(valarray<double>& v, MultipleRHS& mrhs, int flag);

//************************************************************************************************************************

//function to do conversion between an array v (given by double*) and a MultipleRHS object mrhs, use flag to indicate the direction
//flag = 1: MultipleRHS to array
//flag = 2: array to MultipleRHS
//note v.size = mrhs.get_col()*mrhs.get_dim()
void MultipleRHS_Array_Conversion(double* v, MultipleRHS& mrhs, int flag);

//************************************************************************************************************************

//function to calculate the pressure gradient for the staggered grid 
//input:  pressure p as an array of size (Nx+1)*Ny
//output: px and py as array of size (Nx+2)*(Ny-1), but only the first (Nx+1) columns are used
//flag == 0: In-Out flow BC, == 1: Period BC

void STG_pressure_gradient(int Nx, int Ny, double dx, double dy, double* p, valarray<double>& px, valarray<double>& py, int bc_flag);

//************************************************************************************************************************

//function to extrapolate the Node_Data at time n+alpha by data at time n and n-1, if alpha = 1, the extrapolation is at n+1
void Extrapolate_Node_Data(const Node_Data& dn, const Node_Data& dn_1, double dt_n, double dt_n_1, double alpha, Node_Data& d);

//************************************************************************************************************************

//function to interpolate the Node_Data at time n+alpha by data at time n and n+1, if alpha = 1, it is at n+1
void Interpolate_Node_Data(const Node_Data& dnp1, const Node_Data& dn, double alpha, Node_Data& d);

//************************************************************************************************************************

//function to confine every entry of a Node_Data variable within [a,b], mainly used to confine phi and c
void Confine_Node_Data(Node_Data& u, double a, double b);

//************************************************************************************************************************

//this class is the average of a Node_Data variable in x-direction,
//and it has size (N_x+1)x(N_y+1)

class Node_Data_x_Average {

  size_t N_x, N_y;

  valarray< valarray<double> > v;
 
public:

  Node_Data_x_Average(size_t nx, size_t ny);

  ~Node_Data_x_Average() { }

  Node_Data_x_Average(const Node_Data& D); //copy constructor

  size_t get_Nx() const {
    return N_x;
  }

  size_t get_Ny() const {
    return N_y;
  }

  //overload operator[]
  valarray<double>& operator[](int i){
    return v[i];
  }

  //overload operator[]
  const valarray<double> operator[](int i) const {
    return v[i];
  }
};

//************************************************************************************************************************
  
//this class is the average of a Node_Data variable in x-direction,
//and it has size (N_x+2)xN_y

class Node_Data_y_Average {

  size_t N_x, N_y;

  valarray< valarray<double> > v;

 public:

   Node_Data_y_Average(size_t nx, size_t ny);

  ~Node_Data_y_Average() { }

  Node_Data_y_Average(const Node_Data& D); //copy constructor

  size_t get_Nx() const {
    return N_x;
  }

  size_t get_Ny() const {
    return N_y;
  }

  //overload operator[]
  valarray<double>& operator[](int i){
    return v[i];
  }
 
  //overload operator[]
  const valarray<double> operator[](int i) const {
    return v[i];
  }
};

//************************************************************************************************************************

//this class is the average of a Node_Data variable in x-direction, 
//also extend it to index i = -1 and i = Nx+1 for the No-flux boundary condition
//and it has size (N_x+3)x(N_y+1)

class Node_Data_x_Ave_Ext_N {

  size_t N_x, N_y;

  valarray< valarray<double> > v;
 
public:

  Node_Data_x_Ave_Ext_N(size_t nx, size_t ny);

  ~Node_Data_x_Ave_Ext_N() { }

  Node_Data_x_Ave_Ext_N(const Node_Data& D); //copy constructor

  size_t get_Nx() const {
    return N_x;
  }

  size_t get_Ny() const {
    return N_y;
  }

  //overload operator[]
  valarray<double>& operator[](int i){
    return v[i];
  }

  //overload operator[]
  const valarray<double> operator[](int i) const {
    return v[i];
  }
};

//************************************************************************************************************************

//this class is the average of a Node_Data variable in x-direction, 
//also extend it to index i = -1 for the periodic boundary condition
//and it has size (N_x+2)x(N_y+1)

class Node_Data_x_Ave_Ext_P {

  size_t N_x, N_y;

  valarray< valarray<double> > v;
 
public:

  Node_Data_x_Ave_Ext_P(size_t nx, size_t ny);

  ~Node_Data_x_Ave_Ext_P() { }

  Node_Data_x_Ave_Ext_P(const Node_Data& D); //copy constructor

  size_t get_Nx() const {
    return N_x;
  }

  size_t get_Ny() const {
    return N_y;
  }

  //overload operator[]
  valarray<double>& operator[](int i){
    return v[i];
  }

  //overload operator[]
  const valarray<double> operator[](int i) const {
    return v[i];
  }
};

//************************************************************************************************************************

//this class is the average of a Node_Data variable in x-direction, 
//also extend it to index i = -1 and i = Nx+1 for the In-Out flow boundary condition
//and it has size (N_x+3)x(N_y+1)

class Node_Data_x_Ave_Ext_IO {

  size_t N_x, N_y;

  valarray< valarray<double> > v;
 
public:

  Node_Data_x_Ave_Ext_IO(size_t nx, size_t ny);

  ~Node_Data_x_Ave_Ext_IO() { }

  Node_Data_x_Ave_Ext_IO(const Node_Data_Ext_IO_x_N_y& D); //copy constructor

  size_t get_Nx() const {
    return N_x;
  }

  size_t get_Ny() const {
    return N_y;
  }

  //overload operator[]
  valarray<double>& operator[](int i){
    return v[i];
  }

  //overload operator[]
  const valarray<double> operator[](int i) const {
    return v[i];
  }
};

//************************************************************************************************************************
  
//this class is the average of a Node_Data variable in y-direction, 
//also extend it to index j = -1 and i = Ny for the No-flux boundary condition
//and it has size (N_x+2)x(N_y+2)

class Node_Data_y_Ave_Ext_N {

  size_t N_x, N_y;

  valarray< valarray<double> > v;

 public:

   Node_Data_y_Ave_Ext_N(size_t nx, size_t ny);

  ~Node_Data_y_Ave_Ext_N() { }

  Node_Data_y_Ave_Ext_N(const Node_Data& D); //copy constructor
  Node_Data_y_Ave_Ext_N(const Node_Data_Ext_IO_x_N_y& D); //copy constructor 

  size_t get_Nx() const {
    return N_x;
  }

  size_t get_Ny() const {
    return N_y;
  }

  //overload operator[]
  valarray<double>& operator[](int i){
    return v[i];
  }
 
  //overload operator[]
  const valarray<double> operator[](int i) const {
    return v[i];
  }
};

//************************************************************************************************************************

//the class for variables defined on the grid edge, contains two component,
//one for x-direction with size (N_x + 1)x(N_y + 1),
//one for y-direction with size (N_x + 2)x(N_y)
class Edge_Data {

  size_t N_x, N_y;

  valarray< valarray<double> > vx;
  valarray< valarray<double> > vy;

 public:

  Edge_Data(size_t nx, size_t ny);

  ~Edge_Data() { }   


  double& x(int i, int j) {
    return vx[i][j];
  }

  const double& x(int i, int j) const {
    return vx[i][j];
  }

  double& y(int i, int j){
    return vy[i][j];
  }

  const double& y(int i, int j) const {
    return vy[i][j];
  }

  size_t get_Nx() const {
    return N_x;
  }

  size_t get_Ny() const {
    return N_y;
  }

};



#endif
