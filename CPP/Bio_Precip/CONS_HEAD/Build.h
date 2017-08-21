//This file contains the functions building the linear system corresponding to the discretized governing PDE system
#ifndef _BUILD_LINEAR_SYSTEM_H
#define _BUILD_LINEAR_SYSTEM_H

#include "DataStruct.h"
#include "IO.h"
#include <limits>

using namespace std;

//************************************************************************************************************************
//RHS of the NSE equation, No-flux BC in both x and y direction, regular grid
void Build_NSE_RHS_REG_GRID(double* Fx, double* Fy, const Node_Data& rho_n, const Node_Data& mu, const valarray<double>& Rmu,  const valarray<double>& Rphb, 
                            const valarray<double>& Rphc, const Node_Data& un, const Node_Data& vn, const Node_Data& s, const Node_Data& phi_b, const Node_Data& phi_c, 
                            const double& dx, const double& dy, const double& dt, const Control_Parameter& CP);

//************************************************************************************************************************
//RHS of the NSE equation, No-flux BC in both x and y direction, staggered grid

void Build_NSE_RHS_STG_GRID(double* Fx, double* Fy, const Node_Data& rho_n, const Node_Data& mu,  const valarray<double>& Rmu,  const valarray<double>& Rphb, 
                            const valarray<double>& Rphc, const Node_Data& un, const Node_Data& vn, double* s, const Node_Data& phi_b,  const Node_Data& phi_c, 
                            const double& dx, const double& dy, const double& dt, const Control_Parameter& CP);

//************************************************************************************************************************

//function to comput the RHS of the Pressure equation i.e., div(U)
//input:  u, v as the Node_Data, dx and dy are the mesh size,
//bc_flag is the boundary condition flay, calculation is different for different BC
//output: du/dx + dv/dy as an object of MultipleRHS
//here the object of MultipleRHS with dim = Nx + 2 and col_num = Ny + 1;

void Build_Poisson_RHS_REG_GRID(const Node_Data& u, const Node_Data& v, double dx, double dy, MultipleRHS& y, int bc_flag);


//************************************************************************************************************************

//function to comput the RHS of the Pressure equation i.e., div(U) for half-staggered grid (pressure at cell center and velocity at cell corner)
//input:  u, v as the Node_Data, dx and dy are the mesh size, bc_flag is the boundary condition flay, calculation is different for different BC
//output: y = du/dx + dv/dy as an array, the size of y is (Nx+1)*Ny

void Build_Poisson_RHS_STG_GRID(const Node_Data& u, const Node_Data& v, double dx, double dy, double* y, int bc_flag);

//************************************************************************************************************************

// decay rate of [Ca2+], or the precipitation rate of CaCO3
//vol_coef: volume change from [Ca2+] to CaCO3, kp: precipitation rate coef, S: CaCO3 saturation state
double g_c(double vol_coef, double kp, double S);

// growth rate of the network
double g_b(double phi_b, double c);

// consumputation rate of the substrate due to biofilm growth
double g_s(double phi_b, double c);

// nonlinear part of the chemical potential due to phi_c (CaCO3)
double chem_c(double phi_c, double phi_b, double phi_s);

// nonlinear part of the chemical potential due to phi_b (biofilm)
double chem_b(double phi_c, double phi_b, double phi_s);

// nonlinear part of the chemical potential due to phi_s (Solvent)
double chem_s(double phi_c, double phi_b, double phi_s);

//urea hydrolysis rate, return the value at a single position
double hydrolysis(double kur, double phi_b, double Ur);

//CaCO3 precipitation rate, return value at a single position
double precipitate(double kp, double S);

//Calculate saturation state of CaCO3
Node_Data Saturation_State(const Node_Data& Ca, const Node_Data& CO3, const Node_Data& Ion_Stren, int Ion_act_flag);

//Build the growth function for different species, kind of species is indicated by k
double Growth(int k, double phi_b, double Ur, double S, double c);

/* // -(urea reaction rate) + (Ca2+ reaction rate), rate at which dissolved inorganic Carbon increases, */
/* // it is also half of the rate at which [OH-] increases */
/* Node_Data DIC_rate_N_x_N_y(const Node_Data& Ur, const Node_Data& phi_b, const Node_Data& S, double kur, double kp); */

//************************************************************************************************************************

//build the matrix(A) and rhs(F) for the transportation equation of phi_c (CaCO3),
//No-flux BC for phic and Lapacian phic at x = 0 and y-direction
//Out-flow BC for phic and Lapacian phic at x = L 
//Arguments: phib_n: bioiflm volume fraction at previous time step n
//           phic_n: CaCO3 volume fraction at previous time step n

//           un, vn: veloticy at time step n
//           S_n: CaCO3 saturation state at time step n
//This is the modified Cahn-Hilliard equation(MCH) ,
//namely, the velocity is proportional to the chemical potential gradient
void Build_Phi_c_F_x_N_y_MCH(PhiMat& A, valarray<double>& F, const valarray<double>& Rph_np1, 
                           const valarray<double>& Rph_cor, const Node_Data& un, const Node_Data& vn, const  Node_Data& phic_n,
		       const Node_Data& phib_n, const Node_Data& phis_n, const Node_Data& S_n, double dx, double dy, double dt, double alpha);

//************************************************************************************************************************

//build the matrix(A) and rhs(F) for the transportation equation of phi_b (biofilm),
//No-flux BC for phi and Lapacian phi_b in x-direction and y-direction 
//Arguments: phib_n: bioiflm volume fraction at previous time step n
//           phic_n: CaCO3 volume fraction at previous time step n
//           phis_n: solvent volume fraction at previous time step n
//           un, vn: veloticy at time step n+1
//           c_n: substrate concentration extropolated at time n
//This is the modified Cahn-Hilliard equation(MCH) ,
//namely, the velocity is proportional to the chemical potential gradient
void Build_Phi_b_F_x_N_y_MCH(PhiMat& A, valarray<double>& F, const valarray<double>& Rph_np1, const valarray<double>&  Rph_cor, const Node_Data& un, const Node_Data& vn,  const  Node_Data& phib_n,
		       const Node_Data& phic_n, const Node_Data& phis_n, const Node_Data& c_n, double dx, double dy, double dt, double alpha);

/* //\************************************************************************************************************************ */

/* //Crank-Nicolson SCHEME in TIME */
/* //function to bulid the matrix(A) and rhs(F) for the trnasportation equation of the substrate (c) */
/* //No-Flux BC x and y = 0, Dirichelet BC at y = 1 */
/* //arguments: phis_n, phis_np1 : volume fraction of solvent at time step  n and n + 1 */
/* //           phib_n: volume fraction of network at time step n */
/* //           c_n: substrate concentration at time step n  */
/* //           un, vn : velocity at time step n+1 */
/* //           dx, dy, dt */
/* void Build_c_F_x_N_y(c_Mat& A, valarray<double>& F, const valarray<double>& Lc, const valarray<double>& Rc, const Node_Data& phib_n, const Node_Data& phis_np1_sin, const Node_Data& phis_n_sin,  */
/*                        const Node_Data& c_n, const Node_Data& un, const Node_Data& vn, double dx, double dy, double dt, */
/*                        const Control_Parameter& CP); */

/* //\************************************************************************************************************************ */

/* //Crank-Nicolson SCHEME in TIME */
/* //function to bulid the matrix(A) and rhs(F) for the reactioin, advection, diffusion equation of [Urea] */
/* //No-Flux BC in x and y */
/* //arguments: phis_n, phis_np1 : volume fraction of solvent at time step  n and n + 1 */
/* //           Ur_n: Urea concentration at time step n */
/* //           phib_n: biofilm volume fraction at time step n  */
/* //           un, vn : velocity at time step n */
/* //           dx, dy, dt */
/* void Build_Ur_F_x_N_y(c_Mat& A, valarray<double>& F, const valarray<double>& Lc, const valarray<double>& Rc, const Node_Data& phib_n, const Node_Data& phis_np1, const Node_Data& phis_n,  */
/*                        const Node_Data& Ur_n, const Node_Data& un, const Node_Data& vn, double kur, double dx, double dy, double dt); */

/* //\************************************************************************************************************************ */

/* //Crank-Nicolson SCHEME in TIME */
/* //function to bulid the matrix(A) and rhs(F) for the reactioin, advection, diffusion equation of [Ca2+] */
/* //No-Flux BC in x and y */
/* //arguments: phis_n, phis_np1 : volume fraction of solvent at time step  n and n + 1 */
/* //           Ca_n: [Ca2+] concentration at time step n */
/* //           S_n: CaCO3 saturatoin state at time step n  */
/* //           un, vn : velocity at time step n */
/* //           dx, dy, dt */
/* void Build_Ca_F_x_N_y(c_Mat& A, valarray<double>& F, const valarray<double>& Lc, const valarray<double>& Rc, const Node_Data& S_n, const Node_Data& phis_np1, const Node_Data& phis_n,  */
/*                        const Node_Data& Ca_n, const Node_Data& un, const Node_Data& vn, double kp, double dx, double dy, double dt); */
                      
/* //\************************************************************************************************************************ */

/* //Crank-Nicolson SCHEME in TIME */
/* //function to bulid the matrix(A) and rhs(F) for the reactioin, advection, diffusion equation of [OH-] */
/* //No-Flux BC in x and y */
/* //arguments: phis_n, phis_np1 : volume fraction of solvent at time step  n and n + 1 */
/* //           OH_n: [OH-] concentration at time step n */
/* //           OH_rate_n: growth rate of [OH-] at time step n, = 2*DIC_rate  */
/* //           un, vn : velocity at time step n */
/* //           dx, dy, dt */
/* void Build_OH_F_x_N_y(c_Mat& A, valarray<double>& F, const valarray<double>& Lc, const valarray<double>& Rc, const Node_Data& OH_rate_n, const Node_Data& phis_np1, const Node_Data& phis_n,  */
/*                        const Node_Data& OH_n, const Node_Data& un, const Node_Data& vn, double dx, double dy, double dt); */
                       

/* //\************************************************************************************************************************ */

/* //Crank-Nicolson SCHEME in TIME */
/* //function to bulid the matrix(A) and rhs(F) for the reactioin, advection, diffusion equation of [CT] (Dissolved Inorganic Carbon) */
/* //No-Flux BC in x and y */
/* //arguments: phis_n, phis_np1 : volume fraction of solvent at time step  n and n + 1 */
/* //           CT_n: [CT] concentration at time step n */
/* //           CT_rate_n: growth rate of [DIC] at time step n  */
/* //           un, vn : velocity at time step n */
/* //           dx, dy, dt */
/* void Build_CT_F_x_N_y(c_Mat& A, valarray<double>& F, const valarray<double>& Lc, const valarray<double>& Rc, const Node_Data& CT_rate_n, const Node_Data& phis_np1, const Node_Data& phis_n,  */
/*                        const Node_Data& CT_n, const Node_Data& un, const Node_Data& vn, double dx, double dy, double dt); */
//************************************************************************************************************************

//function to build the linear system describing the slow process including slow chemical reactions, convection, diffusion due to concentration gradient
//and electric potential gradient, In-Out flow BC in x and No-Flux BC in y
//output: coefficient matrix (A), rhs (F)
//arguments: k: species index, takes value 1, 2, 3, 4, 6, 7, 8, 9
//           D: variable diffusion coefficient (depending on phi_c), Z: electrical charge of the species
//           C_n: concentration at time step n, Lc, Rc: known values at left and right boundaries due to In-Out flow BC
//           phi_b (biofilm), Ur (Urea), S (Calcite saturation), used in calculating growth rate
//           phis_n, phis_np1 : volume fraction of solvent at time step  n and n + 1
//           un, vn : velocity at time step n
//           EP: electric potential
//           a: unit constant for electric field; dx, dy, dt: grid size in space and time

void Build_Slow_Process(c_Mat& A, valarray<double>& F, int k, const Node_Data& D_n, const Node_Data& D_np1, double Z, const Node_Data& C_n, 
                        const valarray<double>& Lc, const valarray<double>& Rc, const Node_Data& phi_b, const Node_Data& Ur, const Node_Data& S, 
                        const Node_Data& c, const Node_Data& phis_np1_sin, const Node_Data& phis_n_sin, const Node_Data& un, const Node_Data& vn, 
                        const Node_Data& EP, double a, double dx, double dy, double dt, const Control_Parameter& CP);

//************************************************************************************************************************
//same as above, but treat the convection term by upwind scheme

void Build_Slow_Process_Upwind(c_Mat& A, valarray<double>& F, int k, const Node_Data& D_n, const Node_Data& D_np1, double Z, const Node_Data& C_n, 
                        const valarray<double>& Lc, const valarray<double>& Rc, const Node_Data& phi_b, const Node_Data& Ur, const Node_Data& S, 
                        const Node_Data& c, const Node_Data& phis_np1_sin, const Node_Data& phis_n_sin, const Node_Data& un, const Node_Data& vn, 
			const Node_Data& EP, double a, double dx, double dy, double dt, const Control_Parameter& CP);


//************************************************************************************************************************
//same as above, but use physical boundary condtion (balancing flux)

void Build_Slow_Process_FluxBC(c_Mat& A, valarray<double>& F, int k, const Node_Data& D_n, const Node_Data& D_np1, double Z, const Node_Data& C_n,  const valarray<double>& Rc, 
                               const Node_Data& phi_b, const Node_Data& Ur, const Node_Data& S, const Node_Data& c, const Node_Data& phis_np1_sin, 
                               const Node_Data& phis_n_sin, const Node_Data& un, const Node_Data& vn, const Node_Data& EP, 
                               double a, double dx, double dy, double dt, double C_Ini, const Control_Parameter& CP);
//************************************************************************************************************************

//function to bulid the matrix for solving the heat (diffusion) equation with theta-method
//here alph = theta, eps_dif is the diffusion constant

void Build_Diffusion(c_NoFlux_x_NoFlux_y_mat& A, valarray<double>& F, const Node_Data& U_n,
                     double dx, double dy, double dt, double alpha, double eps_dif);

//************************************************************************************************************************

//function to calculate the diffusion coefficient affected by calcite concentration, there is a critical calcite concentration phic_tol,
//below this, diffusion coef D is unaffected, above this, D drops to delta*D quickly
double Diffusion_coef_phic(double D, double phic, double phic_tol, double delta);


//************************************************************************************************************************
 
//function to calculate the diffusion coefficient affected by calcite concentration,
//return the pointwise diffusion coefficient in a Node_Data variable
void Diffusion_Node(Node_Data& DN, const Node_Data& phic, double D, double phic_tol, double delta);


//************************************************************************************************************************
//function to calculate the effective diffusion coefficient by treating calcite as suspended impermeable spheres 
// D/D_0 = (1 - phi_c)/(1 + .5*phi_c)
double Eff_Diff_coef_phic_sis(double D, double phic, double delta);

//************************************************************************************************************************
 
//function to calculate the diffusion coefficient affected by calcite concentration by treating calcite as suspended impermeable spheres ,
//return the pointwise diffusion coefficient in a Node_Data variable
void Eff_Diff_Node_sis(Node_Data& DN, const Node_Data& phic, double D, double delta);


//************************************************************************************************************************
//function to build the matrix and rhs of discretized the electric potential equation (variable coefficient
//Poisson equation), with no-flux BC in x and y direction
//Input: charged species (except [H+] and [OH-]) concentration and their diffusion coefficients (depending on 
//calcite concentration phic), phis, dx, dy, dimensional parameter (a)
void Build_EP_N_x_N_y(c_Mat& A, valarray<double>& F, const Node_Data_Ext_IO_x_N_y& NH4, const Node_Data_Ext_IO_x_N_y& HCO3, const Node_Data_Ext_IO_x_N_y& CO3, 
                      const Node_Data_Ext_IO_x_N_y& Ca, const Node_Data_Ext_IO_x_N_y& Cl, const Node_Data_Ext_IO_x_N_y& Dv_nh4, const Node_Data_Ext_IO_x_N_y& Dv_hco3, 
                      const Node_Data_Ext_IO_x_N_y& Dv_co3, const Node_Data_Ext_IO_x_N_y& Dv_ca, const Node_Data_Ext_IO_x_N_y& Dv_cl, const Node_Data_Ext_IO_x_N_y& phis,
                      double dx, double dy, double a);

//************************************************************************************************************************
//function to build the matrix and rhs of discretized the electric potential equation (variable coefficient
//Poisson equation), with no-flux BC in x and y direction only in the interior domain, excluding i = 0, i = Nx + 1
//Input: charged species (except [H+] and [OH-]) concentration and their diffusion coefficients (depending on 
//calcite concentration phic), phis, dx, dy, dimensional parameter (a)
//source is the integral of the RHS of the equation, check if it's zero
void Build_EP_Interior(c_Mat& A, valarray<double>& F, const Node_Data& NH4, const Node_Data& HCO3, const Node_Data& CO3, const Node_Data& Ca, const Node_Data& Cl, const Node_Data& H, const Node_Data& OH,
                       const Node_Data& Dv_nh4, const Node_Data& Dv_hco3, const Node_Data& Dv_co3, const Node_Data& Dv_ca, const Node_Data& Dv_cl, const Node_Data& Dv_h, const Node_Data& Dv_oh,
                       const Node_Data& phis, double dx, double dy, double dt, double a, double& source, int Inter_file_NO);

//************************************************************************************************************************

//Nonlinear function for [H+] (F([H+] = 0, derived from charge neutral and chemical reaction equilibrium, 
//need to be solved by Newton's method)
double NLF_H(double H, double NT, double K3, double Kw, double CT, double K1, double K2, double Ca, double Cl, double cha_sum);

//Derivative of the previous nonlinear function (NLF_H), used in Newton's method, same for charge neutrality and conservation of [H+]
double Deri_NLF_H(double H, double NT, double K3, double Kw, double CT, double K1, double K2);

//Function to find the root of NLF_H, parameter H holds initial guess at input, and the root as output
//Return 0 if successful, 1 not successful
int Newton_Root(double& H, double NT, double K3, double Kw, double CT, double K1, double K2, double Ca, double Cl, double cha_sum,  int& N_iter);

//function to calculate the concentration of species participating fast reactions
void Cal_Fast_Species(double H, double NT, double K3, double Kw, double CT, double K1, double K2, 
                      double& NH3, double& NH4, double& H2CO3, double& HCO3, double& CO3, double& OH);

//function to calculate the concentration of species participating fast reactions at all grid points
void Cal_Fast_Species_Node(Node_Data& H, const Node_Data& NT, double K3, double Kw, const Node_Data& CT, double K1, double K2, const Node_Data& Ca,
                           const Node_Data& Cl, const Node_Data& Cha_Sum, Node_Data& NH3, Node_Data& NH4, Node_Data& H2CO3, Node_Data& HCO3, Node_Data& CO3, Node_Data& OH, 
                           int Ion_act_flag);

//************************************************************************************************************************
//function to calculate the sum of charges at each point, used as RHS of the fast reaction equation (Charge conservation)
void Cal_Charge_Sum(Node_Data& Cha_Sum, const Node_Data& H, const Node_Data& OH, const Node_Data& Ca, const Node_Data& Cl, 
                                       const Node_Data& NH4, const Node_Data& HCO3, const Node_Data& CO3,  const Node_Data& CaHCO3);

//************************************************************************************************************************

//************************************************************************************************************************
//function to calculate the ionic strength at each point, used for calculating the activity coefficients
void Cal_Ionic_Strength_Node(Node_Data& Ion_Stren, const Node_Data& H, const Node_Data& OH, const Node_Data& Ca, const Node_Data& Cl, 
                                       const Node_Data& NH4, const Node_Data& HCO3, const Node_Data& CO3,  const Node_Data& CaHCO3);

//************************************************************************************************************************
//function to calculate the adjusted equilibrium constants K1_bar, K2_bar, K3_bar, Kw_bar, 
//through ionic strength and activitiy coefficients, only used when ionic strength and activity coefficients are used
void Cal_Equi_Const(double K1, double K2, double K3, double Kw, double K4, double K5, double H, double OH, double Ca, double Cl, double NH4,
		    double HCO3, double CO3, double CaHCO3, double& K3_bar, double& Kw_bar, double& Kb_h2co3, double& Kb_hco3,
                    double& Kb_ca_4, double& Kb_ca_524, double& Kb_ca_52);

//************************************************************************************************************************
#endif
