//this file contains the driver routine for the biofilm model
#ifndef _DRIVER_H
#define _DRIVER_H

#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include "LinearSolver.h"
#include "Build.h"
#include "fish_subroutine.h"
#include "IO.h"

using namespace std;


//************************************************************************************************************************
  
//function to read in the control parameters such as dx, dy , Nx, k (m = 2^(k+1), Ny = m + 1)
//total time T_final, CFL, tolerance of GMRES, etc.
void Read_Control_Parameters(char* filename, Control_Parameter& CP);


//************************************************************************************************************************

//function to read in the global physical parameters, such Gamma_1, Gamma_2, Kai, eta_n, eta_s, rho_n, rho_s, ....
void Read_Global_Parameters(char* filename, Global_Parameter& GP,  const Control_Parameter& CP);

//************************************************************************************************************************

//function to save the intermediate data for later restart, so save the result at time step n and n-1,
//also need to save the time step size at n and n-1 for the use in extrapolation, and the current time
void Save_Inter_Data(char* filename, const Node_Data& un, const Node_Data& vn, const Node_Data& p, const Node_Data& s, double* s_stg, const Node_Data& phi_c, 
                     const Node_Data& phi_b, const Node_Data& Ur, const Node_Data& NH3, const Node_Data& NH4, const Node_Data& H2CO3, const Node_Data& HCO3, 
                     const Node_Data& CO3, const Node_Data& Ca, const Node_Data& Cl, const Node_Data& H, Node_Data& OH, const Node_Data& c, 
                     double T_current, int Inter_file_NO, int n);

//************************************************************************************************************************

//function to read in the saved intermediate data
void Read_Inter_Data(char* filename,  size_t Nx, size_t Ny, Node_Data& un,  Node_Data& vn,  Node_Data& p, Node_Data& s, double* s_stg, Node_Data& phi_c,
                     Node_Data& phi_b, Node_Data& Ur, Node_Data& NH3, Node_Data& NH4, Node_Data& H2CO3, Node_Data& HCO3, Node_Data& CO3, Node_Data& Ca,
		     Node_Data& Cl, Node_Data& H, Node_Data& OH, Node_Data& c, double& T_current, int& Inter_file_NO, int& n);

//************************************************************************************************************************

//function to setup the initial conditinos, use flag to indicate the method to set the initial condition
//flag >= 0, different values corresponding to different initial data, phi may be read in from a  file
void Set_Initial_Condition(char fname[][80], Node_Data& u, Node_Data& v, Node_Data& phi_c, Node_Data& phi_b, Node_Data& Ur, Node_Data& NH3, Node_Data& NH4,
                           Node_Data& H2CO3, Node_Data& HCO3, Node_Data& CO3, Node_Data& Ca,  Node_Data& Cl, Node_Data& H, Node_Data& OH, Node_Data& c, 
                           double dx, double dy, const Control_Parameter& CP);

//************************************************************************************************************************

//function to smooth the initial data by pure diffusion, the input includes the diffusion constant and how many time steps
//will the diffusion run
void Diffuse_Initial_Condition(Node_Data& phi, int N_steps, double dx, double dy, double dt, double alpha, double eps_dif);

//************************************************************************************************************************

//function to update a Node_Data variable with periodic BC in x-direction and, Dirichelet BC  in y-direction, 
//u is a Node_Data and x is solution of the linear system with size (Nx+1)*(Ny-1), including u, v, and phi
void Update_Node_Data_P_x_D_y(Node_Data& u, const valarray<double>& x);

//************************************************************************************************************************

//function to update a Node_Data variable with  in Dirichelet BC  in both x and y direction, 
//u is a Node_Data and x is solution of the linear system with size Nx*(Ny-1), including u, v
void Update_Node_Data_D_x_D_y(Node_Data& u, const valarray<double>& x);

//************************************************************************************************************************
//function to update a Node_Data variable with OutFlow BC in x-direction and, Dirichelet BC  in y-direction, 
//u is a Node_Data and x is solution of the linear system with size Nx*(Ny-1), including u, v
void Update_Node_Data_OutFlow_x_D_y(Node_Data& u, const valarray<double>& x);

//************************************************************************************************************************

//function to update a Node_Data variable with periodic BC in x-drection and no-flux BC in y-direction
//u is a Node_Data and x is solution of the linear system with size (Nx+1)*(Ny+1), including c
void Update_Node_Data_P_x_N_y(Node_Data& u, const valarray<double>& x);

//************************************************************************************************************************

//function to update a Node_Data variable with no-flux BC in x and y direction
//u is a Node_Data and x is solution of the linear system with size (Nx+2)*(Ny+1), including c
void Update_Node_Data_N_x_N_y(Node_Data& u, const valarray<double>& x);

//************************************************************************************************************************

//function to update a Node_Data variable with no-flux BC in x direction and Dirichelet BC in y-direction
//u is a Node_Data and x is solution of the linear system with size (Nx+2)*(Ny-1), including c
void Update_Node_Data_N_x_D_y(Node_Data& u, const valarray<double>& x);

//************************************************************************************************************************
//function to update a Node_Data variable with no-flux BC in x direction and y = 0, Dirichelet BC at y = 1
//u is a Node_Data and x is solution of the linear system with size (Nx+2)*Ny, primarily for c
void Update_Node_Data_N_x_D_N_y(Node_Data& u, const valarray<double>& x);

//************************************************************************************************************************
//function to update a Node_Data variable with periodic BC in x direction, no-flux BC and y = 0, Dirichelet BC at y = 1
//u is a Node_Data and x is solution of the linear system with size (Nx+1)*Ny, primarily for c
void Update_Node_Data_P_x_D_N_y(Node_Data& u, const valarray<double>& x);

//************************************************************************************************************************

//function to update a Node_Data variable with periodic BC in x-direction and, Dirichelet BC  in y-direction
//by adding a vector with size (Nx+2)*(Ny-1) to the original Node_Data variable, mainly used for the pressure 
//correction of the velocity, i.e. , u = u + grad(p)
//u is a Node_Data and x is a vector with size (Nx+2)*(Ny-1), i.e. p_x or p_y
void ADD_Update_Node_Data_P_x_D_y(Node_Data& u, const valarray<double>& x);

//************************************************************************************************************************

//function to update a Node_Data variable with  Dirichelet BC  in both x and y direction
//by adding a vector with size (Nx+2)*(Ny-1) to the original Node_Data variable, mainly used for the pressure 
//correction of the velocity, i.e. , u = u + grad(p), u is a Node_Data and x is a vector with size (Nx+2)*(Ny-1)
//we only update Nx*(Ny-1) entries inside the boundary
void ADD_Update_Node_Data_D_x_D_y(Node_Data& u, const valarray<double>& x);

//************************************************************************************************************************

//function to update a Node_Data variable with  out-flow BC  in both x and Dirichelet BC in y direction
//by adding a vector with size (Nx+2)*(Ny-1) to the original Node_Data variable, mainly used for the pressure 
//correction of the velocity, i.e. , u = u + grad(p), u is a Node_Data and x is a vector with size (Nx+2)*(Ny-1)
//we  update Nx*(Ny-1) entries inside the boundary, and extrapolate it to i = 0 and i = Nx + 1
void ADD_Update_Node_Data_OutFlow_x_D_y(Node_Data& u, const valarray<double>& x);

//************************************************************************************************************************

//fucntion to update the velocity from the solution of the momentum equation
//use bc_flag to distinguish which update_Node_Data function to call
void Update_Velo_By_Momentum(Node_Data& u, const valarray<double>& x, int bc_flag);

//************************************************************************************************************************

//fucntion to update the velocity from the solution of the pressure correction equation
//use bc_flag to distinguish which update_Node_Data function to call
void Update_Velo_By_Pressure(Node_Data& u, const valarray<double>& x, int bc_flag);

//************************************************************************************************************************

//function to update the variable diffusion coefficient of various species depending on the calcite volume fraction phi_c
void Update_Diffusion_Node(Node_Data& D_Ur_n, Node_Data& D_NH3_n, Node_Data& D_NH4_n, Node_Data& D_H2CO3_n,
                           Node_Data& D_HCO3_n, Node_Data& D_CO3_n, Node_Data& D_Ca_n, Node_Data& D_Cl_n, Node_Data& D_H_n,
			   Node_Data& D_OH_n, Node_Data& D_c_n, const Node_Data& phic, const Control_Parameter& CP);

//************************************************************************************************************************


//function to average two Node_Data, the result is put into the second argument, and the first argument
//is constant
void Average_Node_Data(const Node_Data& v, Node_Data& u);

//************************************************************************************************************************

//function to calculate the time step dt at a given node point
//given velocity u, v, sound of speed a, CFL number, dx, dy 
double Explicit_dt(double u, double v, double a, double dx, double dy, double cfl);

//************************************************************************************************************************

//function to calculate the time step for the time integration, 
//which is the smallest dt from all the node points, note here we treat speed of sound as a constant
double Calculate_dt(const Node_Data& u, const Node_Data& v, double a, double dx, double dy, double cfl);

//************************************************************************************************************************

//function to calculate the average viscosity of the biofilm, given the volume fraction of phi_c and phi_b
double Cal_Ave_Viscosity(const Node_Data& phi_c, const Node_Data& phi_b);

//************************************************************************************************************************

//function to determine phi_s from phi_c and phi_b, make sure they sum up to 1 and each is between 0 and 1
//adjust phi_c and phi_b if necessary
void Balance_phi(Node_Data& phi_c, Node_Data& phi_b, Node_Data& phi_s);


//************************************************************************************************************************

//function to calculate the effluent of substance c (integrate it at x = L from y = 0 to y = 1)
double Cal_Effluent(const Node_Data& c, double dy);

//************************************************************************************************************************

//function to calculate the influent of substance c (integrate it at x = 0 from y = 0 to y = 1)  by trapezoidal rule
double Cal_Influent(const Node_Data& c, double dy);

//***********************************************************************************************************************

//the main driver function, integrate over time, filename contains the names of the files containing the parameters
void Time_Evolve(char* filename);

#endif
