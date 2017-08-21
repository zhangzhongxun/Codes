#include "Build.h"

#define DEBUG_Build

//sign function
double sign(double x) {

if (x > 0)
  return 1.0;
else if (x < 0)
  return -1.0;
else
  return 0.0;
}

//************************************************************************************************************************
//RHS of the NSE equation, No-flux BC in both x and y direction
void Build_NSE_RHS_REG_GRID(double* Fx, double* Fy, const Node_Data& rho_n, const Node_Data& mu, const valarray<double>& In_u, 
                            const valarray<double>& Out_u, const valarray<double>& Out_v,const Node_Data& un, const Node_Data& vn, 
                            const Node_Data& s, const Node_Data& phi_b, const Node_Data& phi_c, const double& dx, 
                            const double& dy, const double& dt, const Control_Parameter& CP){


  int i, j, Nx, Ny, curr;

  Nx = un.get_Nx();
  Ny = un.get_Ny();

  Node_Data phi_s(Nx, Ny);

  for(i = 0; i <= Nx + 1; i++)
    for(j = 0; j <= Ny; j++)
      phi_s[i][j] = 1.0 - phi_b[i][j] - phi_c[i][j];

  double dx2, dy2, tmp1, tmp2, tmp3, dE_dphb, dE_dphc, sum_In, sum_Out, Scale;
 
  sum_In = In_u.sum();
  sum_Out = Out_u.sum();

  //If no out flow, use scale as 1, ensure incompressibility, what comes in goes out
  Scale = (sum_Out > 1e-6) ? (sum_In/sum_Out) : 1.0;;
 
  dx2 = dx*dx;
  dy2 = dy*dy;

  //calculate the average of mu in x and y direction

  Node_Data_x_Average mu_x(mu);
  Node_Data_y_Average mu_y(mu);

  //calculate the RHS at the interior nodes
  for(j = 1; j <= Ny - 1; j++)
    for(i = 1; i <= Nx; i++)
      {
	curr = j*(Nx+2)+i;

	//calculate the variational derivative of the mixed energy at node (i,j)
	dE_dphb = -Gamma_1*((phi_b[i+1][j] - 2.0*phi_b[i][j] + phi_b[i-1][j])/dx2 + 
			    (phi_b[i][j+1] - 2.0*phi_b[i][j] + phi_b[i][j-1])/dy2) + 
                	     chem_b(phi_c[i][j], phi_b[i][j], phi_s[i][j]);


	dE_dphc = -Gamma_0*((phi_c[i+1][j] - 2.0*phi_c[i][j] + phi_c[i-1][j])/dx2 + 
			    (phi_c[i][j+1] - 2.0*phi_c[i][j] + phi_c[i][j-1])/dy2) + 
                	     chem_c(phi_c[i][j], phi_b[i][j], phi_s[i][j]);
        

	tmp1 = 2.0*(mu_x[i][j]*(un[i+1][j] - un[i][j]) - mu_x[i-1][j]*(un[i][j] - un[i-1][j]))/dx2 + 
           	   (mu_y[i][j]*(un[i][j+1] - un[i][j]) - mu_y[i][j-1]*(un[i][j] - un[i][j-1]))/dy2 + 
	     + .25*(mu_y[i][j]*(vn[i+1][j+1] + vn[i+1][j] - vn[i-1][j+1] - vn[i-1][j]) - 
	 	    mu_y[i][j-1]*(vn[i+1][j] + vn[i+1][j-1] - vn[i-1][j] - vn[i-1][j-1]))/dx/dy;


        Fx[curr] =  (- rho_n[i][j]*un[i][j]/dt + rho_n[i][j]*(un[i][j]*.5*(un[i+1][j] - un[i-1][j])/dx + vn[i][j]*.5*(un[i][j+1] - un[i][j-1])/dy)
                     - (dE_dphb*(phi_b[i+1][j] - phi_b[i-1][j]) + dE_dphc*(phi_c[i+1][j] - phi_c[i-1][j]))*0.5/dx - tmp1)/eta_ave 
	             + .5*(s[i+1][j] - s[i-1][j])/dx 
	             + ((un[i+1][j] - 2.0*un[i][j] + un[i-1][j])/dx2 + (un[i][j+1] - 2.0*un[i][j] + un[i][j-1])/dy2);


        tmp1 = 2.0*(mu_y[i][j]*(vn[i][j+1] - vn[i][j]) - mu_y[i][j-1]*(vn[i][j] - vn[i][j-1]))/dy2 + 
                   (mu_x[i][j]*(vn[i+1][j] - vn[i][j]) - 
                    mu_x[i-1][j]*(vn[i][j] - vn[i-1][j]))/dx2 + 
	     + .25*(mu_x[i][j]*(un[i+1][j+1] + un[i][j+1] - un[i+1][j-1] - un[i][j-1]) - 
		    mu_x[i-1][j]*(un[i][j+1] + un[i-1][j+1] - 
                                        un[i][j-1] - un[i-1][j-1]))/dx/dy;
	Fy[curr] = (- rho_n[i][j]*vn[i][j]/dt + rho_n[i][j]*(un[i][j]*.5*(vn[i+1][j] - vn[i-1][j])/dx + vn[i][j]*.5*(vn[i][j+1] - vn[i][j-1])/dy) 
                    - (dE_dphb*(phi_b[i][j+1] - phi_b[i][j-1]) + dE_dphc*(phi_c[i][j+1] - phi_c[i][j-1]))*0.5/dy - tmp1)/eta_ave
                    + .5*(s[i][j+1] - s[i][j-1])/dy
        	    + ((vn[i+1][j] - 2.0*vn[i][j] + vn[i-1][j])/dx2 + (vn[i][j+1] - 2.0*vn[i][j] + vn[i][j-1])/dy2); 
       }

  //now treat the boundary condition
  //in y-direction, always Dirichelet BC at y = 0 and y = 1
  for(i = 0; i <= Nx+1; i++)
    {
      j = 0;
      curr = j*(Nx+2) + i;

      Fx[curr] = un[i][j];
      Fy[curr] = vn[i][j];

      j = Ny;
      curr = j*(Nx+2) + i;

      Fx[curr] = un[i][j];
      Fy[curr] = vn[i][j];

    }

  //in x-direction, in-flow and out-flow boundary condition, use drift BC in out flow, it's Dirichlet BC
  for(j = 1; j <= Ny-1; j++)
    {
      //at x = 0

      Fx[j*(Nx+2)] = In_u[j];
      Fy[j*(Nx+2)] = 0.0;
    
      //at x = 1
      Fx[j*(Nx+2)+Nx+1] = Scale*Out_u[j];
      Fy[j*(Nx+2)+Nx+1] = Out_v[j];
    }

}

//************************************************************************************************************************
void Build_NSE_RHS_STG_GRID(double* Fx, double* Fy, const Node_Data& rho_n, const Node_Data& mu, const valarray<double>& In_u, 
                            const valarray<double>& Out_u, const valarray<double>& Out_v, const Node_Data& un, const Node_Data& vn, 
			    double* s, const Node_Data& phi_b,  const Node_Data& phi_c, const double& dx, const double& dy, 
                            const double& dt, const Control_Parameter& CP){


  int i, j, Nx, Ny, curr;

  Nx = un.get_Nx();
  Ny = un.get_Ny();

  Node_Data phi_s(Nx, Ny);

  for(i = 0; i <= Nx + 1; i++)
    for(j = 0; j <= Ny; j++)
      phi_s[i][j] = 1.0 - phi_b[i][j] - phi_c[i][j];

  valarray<double> px(0.0,(Nx+2)*(Ny-1));
  valarray<double> py(0.0,(Nx+2)*(Ny-1));

  double dx2, dy2, tmp1, tmp2, tmp3, dE_dphb, dE_dphc,  sum_In, sum_Out, Scale, y;

  sum_In = In_u.sum();
  sum_Out = Out_u.sum();

  //If no out flow, use scale as 1, ensure incompressibility, what comes in goes out
  Scale = (sum_Out > 1e-6) ? (sum_In/sum_Out) : 1.0;
 
  dx2 = dx*dx;
  dy2 = dy*dy;

  //get the gradient of the pressure, No-flux BC
  STG_pressure_gradient(Nx, Ny, dx, dy, s, px, py, CP.BC_flag); 

  //calculate the average of mu in x and y direction

  Node_Data_x_Average mu_x(mu);
  Node_Data_y_Average mu_y(mu);

  //calculate the RHS at the interior nodes
  for(j = 1; j <= Ny - 1; j++)
    for(i = 1; i <= Nx; i++)
      {
	curr = j*(Nx+2)+i;

	//calculate the variational derivative of the mixed energy at node (i,j)
	dE_dphb = -Gamma_1*((phi_b[i+1][j] - 2.0*phi_b[i][j] + phi_b[i-1][j])/dx2 + 
			    (phi_b[i][j+1] - 2.0*phi_b[i][j] + phi_b[i][j-1])/dy2) + 
                	     chem_b(phi_c[i][j], phi_b[i][j], phi_s[i][j]);


	dE_dphc = -Gamma_0*((phi_c[i+1][j] - 2.0*phi_c[i][j] + phi_c[i-1][j])/dx2 + 
			    (phi_c[i][j+1] - 2.0*phi_c[i][j] + phi_c[i][j-1])/dy2) + 
                	     chem_c(phi_c[i][j], phi_b[i][j], phi_c[i][j]);
        

	tmp1 = 2.0*(mu_x[i][j]*(un[i+1][j] - un[i][j]) - mu_x[i-1][j]*(un[i][j] - un[i-1][j]))/dx2 + 
           	   (mu_y[i][j]*(un[i][j+1] - un[i][j]) - mu_y[i][j-1]*(un[i][j] - un[i][j-1]))/dy2 + 
	     + .25*(mu_y[i][j]*(vn[i+1][j+1] + vn[i+1][j] - vn[i-1][j+1] - vn[i-1][j]) - 
	 	    mu_y[i][j-1]*(vn[i+1][j] + vn[i+1][j-1] - vn[i-1][j] - vn[i-1][j-1]))/dx/dy;


        Fx[curr] =  (- rho_n[i][j]*un[i][j]/dt + rho_n[i][j]*(un[i][j]*.5*(un[i+1][j] - un[i-1][j])/dx + vn[i][j]*.5*(un[i][j+1] - un[i][j-1])/dy)
                     - (dE_dphb*(phi_b[i+1][j] - phi_b[i-1][j]) + dE_dphc*(phi_c[i+1][j] - phi_c[i-1][j]))*0.5/dx - tmp1)/eta_ave 
	             + px[(j-1)*(Nx+2)+i]
	             + ((un[i+1][j] - 2.0*un[i][j] + un[i-1][j])/dx2 + (un[i][j+1] - 2.0*un[i][j] + un[i][j-1])/dy2);


        tmp1 = 2.0*(mu_y[i][j]*(vn[i][j+1] - vn[i][j]) - mu_y[i][j-1]*(vn[i][j] - vn[i][j-1]))/dy2 + 
                   (mu_x[i][j]*(vn[i+1][j] - vn[i][j]) - 
                    mu_x[i-1][j]*(vn[i][j] - vn[i-1][j]))/dx2 + 
	     + .25*(mu_x[i][j]*(un[i+1][j+1] + un[i][j+1] - un[i+1][j-1] - un[i][j-1]) - 
		    mu_x[i-1][j]*(un[i][j+1] + un[i-1][j+1] - 
                                        un[i][j-1] - un[i-1][j-1]))/dx/dy;
	Fy[curr] = (- rho_n[i][j]*vn[i][j]/dt + rho_n[i][j]*(un[i][j]*.5*(vn[i+1][j] - vn[i-1][j])/dx + vn[i][j]*.5*(vn[i][j+1] - vn[i][j-1])/dy) 
                    - (dE_dphb*(phi_b[i][j+1] - phi_b[i][j-1]) + dE_dphc*(phi_c[i][j+1] - phi_c[i][j-1]))*0.5/dy - tmp1)/eta_ave
                    + py[(j-1)*(Nx+2)+i]
        	    + ((vn[i+1][j] - 2.0*vn[i][j] + vn[i-1][j])/dx2 + (vn[i][j+1] - 2.0*vn[i][j] + vn[i][j-1])/dy2); 
       }

  //now treat the boundary condition
  //in y-direction, always Dirichelet BC at y = 0 and y = 1
  for(i = 0; i <= Nx+1; i++)
    {
      j = 0;
      curr = j*(Nx+2) + i;

      Fx[curr] = un[i][j];
      Fy[curr] = vn[i][j];

      j = Ny;
      curr = j*(Nx+2) + i;

      Fx[curr] = un[i][j];
      Fy[curr] = vn[i][j];

    }

  //in x-direction, in-flow and out-flow boundary condition, use drift BC in out flow, it's Dirichlet BC
  for(j = 1; j <= Ny-1; j++)
    {
      //at x = 0

      Fx[j*(Nx+2)] = In_u[j];
      Fy[j*(Nx+2)] = 0.0;
    
      //at x = 1
      Fx[j*(Nx+2)+Nx+1] = Scale*Out_u[j];
      Fy[j*(Nx+2)+Nx+1] = Out_v[j];
    }

}

//************************************************************************************************************************

//function to compute the RHS of the Pressure equation (Poisson equation with Neumann BCs), i.e., div(U)
//input:  u, v as the Node_Data, 
//output: du/dx + dv/dy as an object of MultipleRHS
//note here u, v and div(U) are both of type Node_Data, may need to modify such that div(U) is only defined
//for 1 <= j <= Ny-1

void Build_Poisson_RHS_REG_GRID(const Node_Data& u, const Node_Data& v, double dx, double dy, MultipleRHS& y, int bc_flag){
 
  if(y.get_dim() != (u.get_Nx() + 2) || y.get_col_num() != (u.get_Ny() + 1))
    {
      cout << endl << "Dimension Incompatible in Build_Poisson_RHS(), " 
           << endl << "Error in Build_2D_new.cpp !!" << endl;
      exit(-1);
    }

  if(dx != dy)
    {
      cout << endl << "dx != dy in in Build_Poisson_RHS(), they must be equal in order to apply the "
           << "Fast Poisson Solver. ERROR in Build_2D_new.cpp !!" << endl;
      exit(-1);
    }
 
  int i, j, Nx, Ny;

  Nx = u.get_Nx();
  Ny = u.get_Ny();

  //now compute the inside nodes, same for either BC
  for(j = 1; j <= Ny-1; j++)
    for(i = 1; i <= Nx; i++)
      y[j][i] =  (.5*(u[i+1][j] - u[i-1][j])/dx + .5*(v[i][j+1] - v[i][j-1])/dy);

  //top and bottom boundary except the corner, same for either BC

  j = 0; //the bottom boundary
  for(i = 1; i <= Nx; i++)
    y[j][i] = (.5*(u[i+1][j] - u[i-1][j])/dx + .5*(-3.0*v[i][0] + 4.0*v[i][1] - v[i][2])/dy);  // 2nd order

  j = Ny;//the top boundary
  for(i = 1; i <= Nx; i++)
    y[j][i] = (.5*(u[i+1][j] - u[i-1][j])/dx + .5*(3.0*v[i][Ny] - 4.0*v[i][Ny-1] + v[i][Ny-2])/dy);  // 2nd order

  //The following are at i = 0 and i = Nx+1, calcuation varies according to BC
  if(bc_flag == 0) //no-flux BC, 
    {
      i = 0; // at the left boundary
      {
	for(j = 1; j <= Ny-1; j++)
	  y[j][i] =  (0.5*(-3.0*u[0][j] + 4.0*u[1][j] - u[2][j])/dx + .5*(v[i][j+1] - v[i][j-1])/dy); //2nd order

	j = 0; // lower corner
	y[j][i] =  (0.5*(-3.0*u[0][j] + 4.0*u[1][j] - u[2][j])/dx + .5*(-3.0*v[i][0] + 4.0*v[i][1] - v[i][2])/dy);  

	j = Ny; // upper corner
	y[j][i] =  (0.5*(-3.0*u[0][j] + 4.0*u[1][j] - u[2][j])/dx + .5*(3.0*v[i][Ny] - 4.0*v[i][Ny-1] + v[i][Ny-2])/dy); 
      }

      i = Nx+1; //at the right boundary
      {
	for(j = 1; j <= Ny-1; j++)
	  y[j][i] =  (.5*(3.0*u[i][j] - 4.0*u[i-1][j] + u[i-2][j])/dx + .5*(v[i][j+1] - v[i][j-1])/dy); // 2nd order

	j = 0; // lower corner
	y[j][i] =  (0.5*(3.0*u[i][j] - 4.0*u[i-1][j] + u[i-2][j])/dx + .5*(-3.0*v[i][0] + 4.0*v[i][1] - v[i][2])/dy);  

	j = Ny; // upper corner
	y[j][i] =  (0.5*(3.0*u[i][j] - 4.0*u[i-1][j] + u[i-2][j])/dx + .5*(3.0*v[i][Ny] - 4.0*v[i][Ny-1] + v[i][Ny-2])/dy); 
      }
    }
  else if(bc_flag == 1) //periodic BC in x-direction, 
    {
      i = 0; // at the left boundary
      {
	for(j = 1; j <= Ny-1; j++)
	  y[j][i] =  (0.5*(u[1][j] - u[Nx][j])/dx + .5*(v[i][j+1] - v[i][j-1])/dy); //2nd order

	j = 0; // lower corner
	y[j][i] =  (0.5*(u[1][j] - u[Nx][j])/dx + .5*(-3.0*v[i][0] + 4.0*v[i][1] - v[i][2])/dy);  

	j = Ny; // upper corner
	y[j][i] =  (0.5*(u[1][j] - u[Nx][j])/dx + .5*(3.0*v[i][Ny] - 4.0*v[i][Ny-1] + v[i][Ny-2])/dy); 
      }

      i = Nx+1; //at the right boundary
      {
        y[j][i] = y[j][0];
      }
    }
  else
    {
      cout << endl << " bc_flag in Build_Poisson_RHS() must be 0 or 1, ERROR in BuildLinearSystem.cpp !!! " << endl << endl;
      exit(-1);
    }

}

//************************************************************************************************************************

void Build_Poisson_RHS_STG_GRID(const Node_Data& u, const Node_Data& v, double dx, double dy, double* y, int bc_flag){

  //have to make sure the dimension of y and u, v are compatible, i.e., y.size() = (Nx + 1)*Ny

   int i, j, Nx, Ny;

  Nx = u.get_Nx();
  Ny = u.get_Ny();

  //now compute the inside nodes, same for either BC
  for(j = 0; j <= Ny-1; j++)
    for(i = 0; i <= Nx; i++)
      y[j*(Nx+1)+i] =   0.5*((u[i+1][j] - u[i][j]) + (u[i+1][j+1] - u[i][j+1]))/dx 
	             +  0.5*((v[i][j+1] - v[i][j]) + (v[i+1][j+1] - v[i+1][j]))/dy;  
}

//************************************************************************************************************************

// decay rate of [Ca2+], or the precipitation rate of CaCO3
//vol_coef: volume change from [Ca2+] to CaCO3, kp: precipitation rate coef, S: CaCO3 saturation state
double g_c(double vol_coef, double kp, double S)
{
  return -vol_coef*precipitate(kp,S);
}

//function model the growth rate of the network
double g_b(double phi_b, double c)
{
  return sub_epsilon*sub_mu*phi_b*c/(sub_K_c + c);
}

// consumputation rate of the substrate due to biofilm growth
double g_s(double phi_b, double c)
{
  return -phi_b*sub_A*c/(sub_K_c + c);
}

// nonlinear part of the chemical potential due to phi_c (CaCO3)
double chem_c(double phi_c, double phi_b, double phi_s)
{
  if(chemc_flag == 0)
    return Ca_phase*2.0*phi_c*(1.0 - phi_c)*(1.0 - 2.0*phi_c) + 2.0*Vol_Penal*(phi_c + phi_b + phi_s - 1.0);
  else if(chemc_flag == 1)
    return Ca_phase*(1.0 - 2.0*phi_c) + (1e-1)*Gamma_2*((log(phi_c + PHI_reg) + 1.0) - log(1.0 - phi_c + PHI_reg) - 1) 
           + 2.0*Vol_Penal*(phi_c + phi_b + phi_s - 1.0);
  else 
    {
      cout << endl << " chemc_flag = " << chemc_flag << ", it must be 0 or 1, error at line " << __LINE__ << " of file " << __FILE__ << endl;
      exit(-1);
    }
}

// nonlinear part of the chemical potential due to phi_b (biofilm)
// This part is tricky and may cause problem !!!
double chem_b(double phi_c, double phi_b, double phi_s)
{ 		    
  //  return Gamma_2*(NP_inv*(log(phi_b + PHI_reg) + 1.0) + Kai_CH*phi_s) + 2.0*Vol_Penal*(phi_c + phi_b + phi_s - 1.0);
  return Gamma_2*(NP_inv*(log(phi_b + PHI_reg) + 1.0) - log(1.0 - phi_b + PHI_reg) - 1 + Kai_CH*(1.0 - 2*phi_b)) + 2.0*Vol_Penal*(phi_c + phi_b + phi_s - 1.0);
}

// nonlinear part of the chemical potential due to phi_s (Solvent)
double chem_s(double phi_c, double phi_b, double phi_s)
{
  return Gamma_2*(log(phi_s + PHI_reg) + 1 + Kai_CH*phi_b) + 2.0*Vol_Penal*(phi_c + phi_b + phi_s - 1.0);
}

//urea hydrolysis rate, return the value at a single position
double hydrolysis(double kur, double phi_b, double Ur)
{
  return (-kur*phi_b*Ur);
}

//CaCO3 precipitation rate, return value at a single position
double precipitate(double kp, double S)
{
 //  return  (S < S_crit) ? 0.0 : (-kp*(S - 1.0)*(S - 1.0));
  return .5*(tanh(S - S_crit) + 1.0)*(-kp*(S - 1.0)*(S - 1.0));
} 

//************************************************************************************************************************
//Calculate saturation state of CaCO3
Node_Data Saturation_State(const Node_Data& Ca, const Node_Data& CO3){

  int i, j, Nx, Ny;

  Nx = Ca.get_Nx();
  Ny = Ca.get_Ny();

  Node_Data S(Nx, Ny);

  for(i = 0; i <= Nx + 1; i++)
    for(j = 0; j <= Ny; j++)
      {
        //note the unit of [Ca2+] and [CT] are mM (milli-mole), need to multipley by 1e-3
	S[i][j] = Ca[i][j]*CO3[i][j]/K_so;
      }

  return S;
}

//************************************************************************************************************************
// auxiliary function to creat a consumption rate on the right half domain
// double Ind_Consump(double x, double y) {

//   if(x >= 2.0)
//     return 1.0;
//   else 
//     return 0.0;

// }


// double Growth(int k, double phi_b, double Ur, double S, double c, double x, double y){

//   double result;

//   if(k == 1)  // [Urea]
//     result =  hydrolysis(k_ur, phi_b, Ur);
//   else if(k == 2) // [NH3]
//     result = -2.0*hydrolysis(k_ur, phi_b, Ur);
//   else if(k == 3) // [H2CO3]
//     result = -hydrolysis(k_ur, phi_b, Ur);
//   else if(k == 4 || k == 6 || k == 9) // [NH4+], [HCO3-], [Cl-]
//     result = 0.0;
//   else if(k == 7 || k == 8) // [Ca2+], [CO3 2-]
//     result = precipitate(k_p, S);
//   else if(k == 11)
//     result = g_s(phi_b, c);
//   else
//     {
//       cout << endl << "Species index k not in the correct range, error in function Growth() at line " << __LINE__ << " of file " << __FILE__ << endl;
//       exit(-1);
//     }

//   return result;
// }

double Growth(int k, double phi_b, double Ur, double S, double c){

  double result;

  if(k == 1)  // [Urea]
    result =  hydrolysis(k_ur, phi_b, Ur);
  else if(k == 2) // [NH3]
    result = -2.0*hydrolysis(k_ur, phi_b, Ur);
  else if(k == 3) // [H2CO3]
    result = -hydrolysis(k_ur, phi_b, Ur);
  else if(k == 4 || k == 6 || k == 9 || k == 5 || k == 10) // [NH4+], [HCO3-], [Cl-], [OH-], [H+]
    result = 0.0;
  else if(k == 7 || k == 8) // [Ca2+], [CO3 2-]
    result = precipitate(k_p, S);
  else if(k == 11)
    result = g_s(phi_b, c);
  else
    {
      cout << endl << "Species index k not in the correct range, error in function Growth() at line " << __LINE__ << " of file " << __FILE__ << endl;
      exit(-1);
    }

  return result;
}

//************************************************************************************************************************
void Build_Phi_b_F_x_N_y_MCH(PhiMat& A, valarray<double>& F, const valarray<double>& Rph_np1, const valarray<double>& Rph_cor,
                             const Node_Data& un, const Node_Data& vn,  const  Node_Data& phib_n, const Node_Data& phic_n, 
                             const Node_Data& phis_n, const Node_Data& c_n, double dx, double dy, double dt, double alpha){

  int i, j, Nx, Ny, center, curr, LX;

  Nx = A.get_Nx();
  Ny = A.get_Ny();
  LX = A.get_LX();

  if(A.get_dim() != F.size() || Rph_np1.size() != Ny + 1 || Rph_cor.size() != Ny + 1)
    {
      cout << endl << "Dimension incompatible in Build_Phi_Flow_x_N_y_MCH(), ERROR at line " << __LINE__ 
           << " of file " << __FILE__ << " !! "  << endl;
      exit(-1);
    }

  char debfile[80];


  double Idx2, Idy2, Idx4, Idy4, tmp;

  Idx2 = 1.0/(dx*dx);
  Idx4 = Idx2*Idx2;
  Idy2 = 1.0/(dy*dy);
  Idy4 = Idy2*Idy2;

  Node_Data_Ext_N_x_N_y phb_E(phib_n);
  Node_Data_Ext_N_x_N_y phc_E(phic_n);
  Node_Data_Ext_N_x_N_y phs_E(phis_n);

  Node_Data mob(Nx, Ny); //used to calculate the mobility parameter

  //a_11
  for(i = 0; i <= Nx + 1; i++)
    for(j = 0; j <= Ny; j++)
      {
        tmp = Lambda_c*phic_n[i][j]*(1.0 - phic_n[i][j]);
        mob[i][j] = (tmp >= 0) ? tmp : 0.0;
      }

  Node_Data_x_Ave_Ext_N a_11_x(mob);
  Node_Data_y_Ave_Ext_N a_11_y(mob);

//   //*************************************DEBUG*********************************************
//   sprintf(debfile,"%sa_11_x",SaveDir);
//   Write_Node_Data_x_Ave_Ext_N(debfile, a_11_x);
//   sprintf(debfile,"%sa_11_y",SaveDir);
//   Write_Node_Data_y_Ave_Ext_N(debfile, a_11_y);
//   //*************************************DEBUG*********************************************

  //a_22
  for(i = 0; i <= Nx + 1; i++)
    for(j = 0; j <= Ny; j++)
      {
        tmp = Lambda_b*phib_n[i][j]*(1.0 - phib_n[i][j]);
        mob[i][j] = (tmp >= 0) ? tmp : 0.0;
      }

  Node_Data_x_Ave_Ext_N a_22_x(mob);
  Node_Data_y_Ave_Ext_N a_22_y(mob);

//   //*************************************DEBUG*********************************************
//   sprintf(debfile,"%sa_22_x",SaveDir);
//   Write_Node_Data_x_Ave_Ext_N(debfile, a_22_x);
//   sprintf(debfile,"%sa_22_y",SaveDir);
//   Write_Node_Data_y_Ave_Ext_N(debfile, a_22_y);
//   //*************************************DEBUG*********************************************


  //a_33
  for(i = 0; i <= Nx + 1; i++)
    for(j = 0; j <= Ny; j++)
      {
        tmp = Lambda_s*phis_n[i][j]*(1.0 - phis_n[i][j]);
        mob[i][j] = (tmp >= 0) ? tmp : 0.0;
      }

  Node_Data_x_Ave_Ext_N a_33_x(mob);
  Node_Data_y_Ave_Ext_N a_33_y(mob);

//   //*************************************DEBUG*********************************************
//   sprintf(debfile,"%sa_33_x",SaveDir);
//   Write_Node_Data_x_Ave_Ext_N(debfile, a_33_x);
//   sprintf(debfile,"%sa_33_y",SaveDir);
//   Write_Node_Data_y_Ave_Ext_N(debfile, a_33_y);
//   //*************************************DEBUG*********************************************

  //a_12
  Node_Data_x_Ave_Ext_N a_12_x(Nx, Ny);
  Node_Data_y_Ave_Ext_N a_12_y(Nx, Ny);

  for(i = 0; i <= Nx + 2; i++)
    for(j = 0; j <= Ny; j++)
      a_12_x[i][j] = .5*(a_33_x[i][j] - a_11_x[i][j] - a_22_x[i][j]);

  for(i = 0; i <= Nx + 1; i++)
    for(j = 0; j <= Ny + 1; j++)
      a_12_y[i][j] = .5*(a_33_y[i][j] - a_11_y[i][j] - a_22_y[i][j]);

//   //*************************************DEBUG*********************************************
//   sprintf(debfile,"%sa_12_x",SaveDir);
//   Write_Node_Data_x_Ave_Ext_N(debfile, a_12_x);
//   sprintf(debfile,"%sa_12_y",SaveDir);
//   Write_Node_Data_y_Ave_Ext_N(debfile, a_12_y);
//   //*************************************DEBUG*********************************************

  //a_23
  Node_Data_x_Ave_Ext_N a_23_x(Nx, Ny);
  Node_Data_y_Ave_Ext_N a_23_y(Nx, Ny);

  for(i = 0; i <= Nx + 2; i++)
    for(j = 0; j <= Ny; j++)
      a_23_x[i][j] = .5*(a_11_x[i][j] - a_22_x[i][j] - a_33_x[i][j]);

  for(i = 0; i <= Nx + 1; i++)
    for(j = 0; j <= Ny + 1; j++)
      a_23_y[i][j] = .5*(a_11_y[i][j] - a_22_y[i][j] - a_33_y[i][j]);

//   //*************************************DEBUG*********************************************
//   sprintf(debfile,"%sa_23_x",SaveDir);
//   Write_Node_Data_x_Ave_Ext_N(debfile, a_23_x);
//   sprintf(debfile,"%sa_23_y",SaveDir);
//   Write_Node_Data_y_Ave_Ext_N(debfile, a_23_y);
//   //*************************************DEBUG*********************************************

  // vector of length 13 to put the coefficients from the \Delta\Delta phi
  valarray<double> DLAP(13);  //double Laplacian

  // vector of length 5 to put the coefficients fromt the div(phi*grad(phi)) term
  valarray<double> LAP(5);


  //***********************DEBUG***************************
  double bW[7] = {1.0, 1.0, 1.0, 1.0, .0, .0, .0};
//   double bW[7] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  //***********************DEBUG***************************

 
  //Calculation of A[i][j] and  the RHS, it is same for any (i,j) since we use extended data structure containing
  //nodes outside the boundary 
  for(j = 0; j <= Ny; j++)
    {
      center = j*13 + 6; 
      for(i = 0; i <= Nx; i++)
	{
	  DLAP[0]  = Idy4*a_22_y[i][j];
	  DLAP[1]  = Idx2*Idy2*(a_22_x[i][j] + a_22_y[i][j]);
	  DLAP[2]  = - (Idx2*Idy2*(a_22_x[i+1][j] + a_22_x[i][j]) + Idy4 *a_22_y[i][j+1] + Idy2*(3.0*Idy2 + 2.0*Idx2)*a_22_y[i][j]);
	  DLAP[3]  = Idx2*Idy2*(a_22_x[i+1][j] + a_22_y[i][j]);
	  DLAP[4]  = Idx4*a_22_x[i][j];
	  DLAP[5]  = -(Idx4*a_22_x[i+1][j] + Idx2*(3.0*Idx2 + 2.0*Idy2)*a_22_x[i][j] + Idx2*Idy2*(a_22_y[i][j+1] + a_22_y[i][j]));
	  DLAP[6]  = (Idx2*(3.0*Idx2 + 2.0*Idy2)*(a_22_x[i+1][j] + a_22_x[i][j]) + 
                      Idy2*(3.0*Idy2 + 2.0*Idx2)*(a_22_y[i][j+1] + a_22_y[i][j]));
	  DLAP[7]  = -(Idx2*(3.0*Idx2 + 2.0*Idy2)*a_22_x[i+1][j] + Idx4*a_22_x[i][j] + Idx2*Idy2*(a_22_y[i][j+1] + a_22_y[i][j]));
	  DLAP[8]  = Idx4*a_22_x[i+1][j];
	  DLAP[9]  = Idx2*Idy2*(a_22_x[i][j] + a_22_y[i][j+1]);
	  DLAP[10] = -(Idx2*Idy2*(a_22_x[i+1][j] + a_22_x[i][j]) + Idy4*a_22_y[i][j] + Idy2*(2.0*Idx2 + 3.0*Idy2)*a_22_y[i][j+1]);
	  DLAP[11] = Idx2*Idy2*(a_22_x[i+1][j] + a_22_y[i][j+1]);
	  DLAP[12] = Idy4*a_22_y[i][j+1];

	  A[center][i]   = 1.0/dt + alpha*Gamma_1*DLAP[6];
	  A[center-1][i] = -.5*un[i][j]/dx + alpha*Gamma_1* DLAP[5];
	  A[center-2][i] = alpha*Gamma_1*DLAP[4];
	  A[center-3][i] = alpha*Gamma_1*DLAP[3];
	  A[center-4][i] = -.5*vn[i][j]/dy + alpha*Gamma_1*DLAP[2];
	  A[center-5][i] = alpha*Gamma_1*DLAP[1];
	  A[center-6][i] = alpha*Gamma_1*DLAP[0];
	  A[center+1][i] = .5*un[i][j]/dx + alpha*Gamma_1* DLAP[7];
	  A[center+2][i] = alpha*Gamma_1*DLAP[8];
	  A[center+3][i] = alpha*Gamma_1*DLAP[9];
	  A[center+4][i] = .5*vn[i][j]/dy + alpha*Gamma_1*DLAP[10];
	  A[center+5][i] = alpha*Gamma_1*DLAP[11];
	  A[center+6][i] = alpha*Gamma_1*DLAP[12];

	  curr = j*LX + i;
      
	  F[curr] = bW[0]*phib_n[i][j]/dt + bW[1]*g_b(phib_n[i][j], c_n[i][j]);

	  tmp = (alpha - 1.0)*Gamma_1*(DLAP[0]*phb_E.get(i,j-2) +
                           DLAP[1]*phb_E.get(i-1,j-1) + DLAP[2]*phb_E.get(i,j-1)    + DLAP[3]*phb_E.get(i+1,j-1) +
			   DLAP[4]*phb_E.get(i-2,j)   + DLAP[5]*phb_E.get(i-1,j)    + DLAP[6]*phb_E.get(i,j)  + 
			   DLAP[7]*phb_E.get(i+1,j)   + DLAP[8]*phb_E.get(i+2,j)    + DLAP[9]*phb_E.get(i-1,j+1) +
			   DLAP[10]*phb_E.get(i,j+1)  + DLAP[11]*phb_E.get(i+1,j+1) + DLAP[12]*phb_E.get(i,j+2) ); 

	  F[curr] += bW[2]*tmp;  //the double Laplacian term due to phi_b

          LAP[0] = Idy2*a_22_y[i][j];
          LAP[1] = Idx2*a_22_x[i][j];
	  LAP[2] = -(Idx2*(a_22_x[i+1][j] + a_22_x[i][j]) + Idy2*(a_22_y[i][j+1] + a_22_y[i][j]));
          LAP[3] = Idx2*a_22_x[i+1][j];
          LAP[4] = Idy2*a_22_y[i][j+1];
 
	  tmp = LAP[0]*chem_b(phc_E.get(i,j-1), phb_E.get(i,j-1), phs_E.get(i,j-1)) + LAP[1]*chem_b(phc_E.get(i-1,j), phb_E.get(i-1,j), phs_E.get(i-1,j)) + 
                LAP[2]*chem_b(phc_E.get(i,j), phb_E.get(i,j), phs_E.get(i,j)) +
                LAP[3]*chem_b(phc_E.get(i+1,j), phb_E.get(i+1,j), phs_E.get(i+1,j)) + LAP[4]*chem_b(phc_E.get(i,j+1), phb_E.get(i,j+1), phs_E.get(i,j+1));

	  F[curr] += bW[3]*tmp;  //the Laplacian term due to phi_b

          LAP[0] = Idy2*a_23_y[i][j];
          LAP[1] = Idx2*a_23_x[i][j];
	  LAP[2] = -(Idx2*(a_23_x[i+1][j] + a_23_x[i][j]) + Idy2*(a_23_y[i][j+1] + a_23_y[i][j]));
          LAP[3] = Idx2*a_23_x[i+1][j];
          LAP[4] = Idy2*a_23_y[i][j+1];
 
	  tmp = LAP[0]*chem_s(phc_E.get(i,j-1), phb_E.get(i,j-1), phs_E.get(i,j-1)) + LAP[1]*chem_s(phc_E.get(i-1,j), phb_E.get(i-1,j), phs_E.get(i-1,j)) + 
                LAP[2]*chem_s(phc_E.get(i,j), phb_E.get(i,j), phs_E.get(i,j)) +
                LAP[3]*chem_s(phc_E.get(i+1,j), phb_E.get(i+1,j), phs_E.get(i+1,j)) + LAP[4]*chem_s(phc_E.get(i,j+1), phb_E.get(i,j+1), phs_E.get(i,j+1));

	  F[curr] += bW[4]*tmp;  //the Laplacian term due to phi_s

          LAP[0] = Idy2*a_12_y[i][j];
          LAP[1] = Idx2*a_12_x[i][j];
	  LAP[2] = -(Idx2*(a_12_x[i+1][j] + a_12_x[i][j]) + Idy2*(a_12_y[i][j+1] + a_12_y[i][j]));
          LAP[3] = Idx2*a_12_x[i+1][j];
          LAP[4] = Idy2*a_12_y[i][j+1];
 
	  tmp = LAP[0]*chem_c(phc_E.get(i,j-1), phb_E.get(i,j-1), phs_E.get(i,j-1)) + LAP[1]*chem_c(phc_E.get(i-1,j), phb_E.get(i-1,j), phs_E.get(i-1,j)) + 
                LAP[2]*chem_c(phc_E.get(i,j), phb_E.get(i,j), phs_E.get(i,j)) +
                LAP[3]*chem_c(phc_E.get(i+1,j), phb_E.get(i+1,j), phs_E.get(i+1,j)) + LAP[4]*chem_c(phc_E.get(i,j+1), phb_E.get(i,j+1), phs_E.get(i,j+1));

	  F[curr] += bW[5]*tmp;  //the Laplacian term due to phi_c

	  DLAP[0]  = Idy4*a_12_y[i][j];
	  DLAP[1]  = Idx2*Idy2*(a_12_x[i][j] + a_12_y[i][j]);
	  DLAP[2]  = - (Idx2*Idy2*(a_12_x[i+1][j] + a_12_x[i][j]) + Idy4 *a_12_y[i][j+1] + Idy2*(3.0*Idy2 + 2.0*Idx2)*a_12_y[i][j]);
	  DLAP[3]  = Idx2*Idy2*(a_12_x[i+1][j] + a_12_y[i][j]);
	  DLAP[4]  = Idx4*a_12_x[i][j];
	  DLAP[5]  = -(Idx4*a_12_x[i+1][j] + Idx2*(3.0*Idx2 + 2.0*Idy2)*a_12_x[i][j] + Idx2*Idy2*(a_12_y[i][j+1] + a_12_y[i][j]));
	  DLAP[6]  = (Idx2*(3.0*Idx2 + 2.0*Idy2)*(a_12_x[i+1][j] + a_12_x[i][j]) + 
                      Idy2*(3.0*Idy2 + 2.0*Idx2)*(a_12_y[i][j+1] + a_12_y[i][j]));
	  DLAP[7]  = -(Idx2*(3.0*Idx2 + 2.0*Idy2)*a_12_x[i+1][j] + Idx4*a_12_x[i][j] + Idx2*Idy2*(a_12_y[i][j+1] + a_12_y[i][j]));
	  DLAP[8]  = Idx4*a_12_x[i+1][j];
	  DLAP[9]  = Idx2*Idy2*(a_12_x[i][j] + a_12_y[i][j+1]);
	  DLAP[10] = -(Idx2*Idy2*(a_12_x[i+1][j] + a_12_x[i][j]) + Idy4*a_12_y[i][j] + Idy2*(2.0*Idx2 + 3.0*Idy2)*a_12_y[i][j+1]);
	  DLAP[11] = Idx2*Idy2*(a_12_x[i+1][j] + a_12_y[i][j+1]);
	  DLAP[12] = Idy4*a_12_y[i][j+1];

	  tmp = - Gamma_0*(DLAP[0]*phc_E.get(i,j-2) +
                           DLAP[1]*phc_E.get(i-1,j-1) + DLAP[2]*phc_E.get(i,j-1)    + DLAP[3]*phc_E.get(i+1,j-1) +
			   DLAP[4]*phc_E.get(i-2,j)   + DLAP[5]*phc_E.get(i-1,j)    + DLAP[6]*phc_E.get(i,j)  + 
			   DLAP[7]*phc_E.get(i+1,j)   + DLAP[8]*phc_E.get(i+2,j)    + DLAP[9]*phc_E.get(i-1,j+1) +
			   DLAP[10]*phc_E.get(i,j+1)  + DLAP[11]*phc_E.get(i+1,j+1) + DLAP[12]*phc_E.get(i,j+2) ); 

	  F[curr] += bW[6]*tmp; //the double Lapacian term due to phi_c
	}
    }

  //now enforce the Boundary conditions, no flux at x = 0, y = 0, y = 1; Dirichlet at x = L

  j = 0;
  center = j*13 + 6; 

  i = 0;

  A[center+6][i] += A[center-6][i];
  A[center+4][i] += A[center-4][i];      
  A[center+5][i] += A[center-3][i] +  A[center+3][i] + A[center-5][i];    
  A[center+1][i] += A[center-1][i];
  A[center+2][i] += A[center-2][i];   

  i = 1;

  A[center+6][i] += A[center-6][i];
  A[center+3][i] += A[center-5][i];      
  A[center+4][i] += A[center-4][i];      
  A[center+5][i] += A[center-3][i]; 
  A[center][i]   += A[center-2][i];   

  for(i = 2; i <= Nx-2; i++)
    {
      A[center+6][i] += A[center-6][i];
      A[center+3][i] += A[center-5][i];      
      A[center+4][i] += A[center-4][i];      
      A[center+5][i] += A[center-3][i]; 
    } 

  i = Nx-1;

  A[center+6][i] += A[center-6][i];
  A[center+3][i] += A[center-5][i];      
  A[center+4][i] += A[center-4][i];      
  A[center+5][i] += A[center-3][i]; 

  curr = j*LX + i;

  F[curr] -= A[center+2][i]*Rph_np1[j];

  i = Nx;

  A[center+6][i] += A[center-6][i];
  A[center+4][i] += A[center-4][i];      
  A[center+3][i] += A[center-5][i];
  A[center][i] -= A[center+2][i];
 
  curr = j*LX + i;

  F[curr] -= A[center+2][i]*Rph_cor[j] + A[center+1][i]*Rph_np1[j] + 
             A[center-3][i]*Rph_np1[j+1] + A[center+5][i]*Rph_np1[j+1];

  j = 1; 
  center = j*13 + 6; 

  i = 0;

  A[center+1][i] += A[center-1][i];
  A[center+2][i] += A[center-2][i];      
  A[center+5][i] += A[center+3][i];      
  A[center-3][i] += A[center-5][i];      
  A[center][i]   += A[center-6][i];

  i = 1;

  A[center][i] += A[center-6][i] + A[center-2][i];

  for(i = 2; i <= Nx-2; i++)
    {
      A[center][i] += A[center-6][i];
    } 

  i = Nx - 1;

  A[center][i] += A[center-6][i]; 

  curr = j*LX + i;
  F[curr] -= A[center+2][i]*Rph_np1[j];

  i = Nx;

  A[center][i] += A[center-6][i] - A[center+2][i];

  curr = j*LX + i;
  F[curr] -= A[center+1][i]*Rph_np1[j] + A[center+2][i]*Rph_cor[j] + A[center-3][i]*Rph_np1[j-1] + A[center+5][i]*Rph_np1[j+1];


  for(j = 2; j <= Ny-2; j++)
    {
      center = j*13 + 6; 

      i = 0;

      A[center+1][i] += A[center-1][i];
      A[center+2][i] += A[center-2][i];      
      A[center+5][i] += A[center+3][i];      
      A[center-3][i] += A[center-5][i];      

      i = 1;

      A[center][i] += A[center-2][i];

      i = Nx - 1;
      
      curr = j*LX + i;
      F[curr] -= A[center+2][i]*Rph_np1[j];

      i = Nx;

      A[center][i] -= A[center+2][i];

      curr = j*LX + i;
      F[curr] -= A[center+1][i]*Rph_np1[j] + A[center+2][i]*Rph_cor[j] + A[center-3][i]*Rph_np1[j-1] + A[center+5][i]*Rph_np1[j+1];
    }

  j = Ny-1 ;
  center = j*13 + 6; 

  i = 0;

  A[center+1][i] += A[center-1][i];
  A[center+2][i] += A[center-2][i];      
  A[center+5][i] += A[center+3][i];      
  A[center-3][i] += A[center-5][i];      
  A[center][i]   += A[center+6][i];

  i = 1;
  A[center][i] += A[center+6][i] + A[center-2][i];

  for(i = 2; i <= Nx-2; i++)
    {
      A[center][i] += A[center+6][i];
    }

  i = Nx - 1;

  A[center][i] += A[center+6][i];

  curr = j*LX + i;
  F[curr] -= A[center+2][i]*Rph_np1[j];

  i = Nx;

  A[center][i] += A[center+6][i] - A[center+2][i];

  curr = j*LX + i;
  F[curr] -= A[center+1][i]*Rph_np1[j] + A[center+2][i]*Rph_cor[j] + A[center+5][i]*Rph_np1[j+1] + A[center-3][i]*Rph_np1[j-1];

  j = Ny;
  center = j*13 + 6; 

  i = 0;

  A[center-6][i] += A[center+6][i];
  A[center-4][i] += A[center+4][i];      
  A[center-3][i] += A[center+5][i] +  A[center+3][i] + A[center-5][i];    
  A[center+1][i] += A[center-1][i];
  A[center+2][i] += A[center-2][i];   

  i = 1;

  A[center-6][i] += A[center+6][i];
  A[center-5][i] += A[center+3][i];      
  A[center-4][i] += A[center+4][i];      
  A[center-3][i] += A[center+5][i];
  A[center][i]   += A[center-2][i];   

  for(i = 2; i <= Nx-2; i++)
    {
      A[center-6][i] += A[center+6][i];
      A[center-5][i] += A[center+3][i];      
      A[center-4][i] += A[center+4][i];      
      A[center-3][i] += A[center+5][i]; 
    } 

  i = Nx - 1;

  A[center-6][i] += A[center+6][i];
  A[center-5][i] += A[center+3][i];      
  A[center-4][i] += A[center+4][i];      
  A[center-3][i] += A[center+5][i];

  curr = j*LX + i;
  F[curr] -= A[center+2][i]*Rph_np1[j];
 

  i = Nx;

  A[center-6][i] += A[center+6][i];
  A[center-4][i] += A[center+4][i];      
  A[center-5][i] += A[center+3][i];
  A[center][i] -= A[center+2][i];

  curr = j*LX + i;
  F[curr] -= A[center+1][i]*Rph_np1[j] + A[center+2][i]*Rph_cor[j] + A[center+5][i]*Rph_np1[j-1] + A[center-3][i]*Rph_np1[j-1];
 
}

//************************************************************************************************************************
void Build_Phi_c_F_x_N_y_MCH(PhiMat& A, valarray<double>& F, const valarray<double>& Rph_np1, const valarray<double>& Rph_cor, 
                             const Node_Data& un, const Node_Data& vn, const  Node_Data& phic_n, const Node_Data& phib_n, const Node_Data& phis_n, 
                             const Node_Data& S_n, double dx, double dy, double dt, double alpha){

  int i, j, Nx, Ny, center, curr, LX;

  Nx = A.get_Nx();
  Ny = A.get_Ny();
  LX = A.get_LX();

  if(A.get_dim() != F.size() || Rph_np1.size() != Ny + 1 || Rph_cor.size() != Ny + 1)
    {
      cout << endl << "Dimension incompatible in Build_Phi_Flow_x_N_y_MCH(), ERROR at line " << __LINE__ 
           << " of file " << __FILE__ << " !! "  << endl;
      exit(-1);
    }

  char debfile[80];

  double Idx2, Idy2, Idx4, Idy4, tmp;

  Idx2 = 1.0/(dx*dx);
  Idx4 = Idx2*Idx2;
  Idy2 = 1.0/(dy*dy);
  Idy4 = Idy2*Idy2;

  Node_Data_Ext_N_x_N_y phb_E(phib_n);
  Node_Data_Ext_N_x_N_y phc_E(phic_n);
  Node_Data_Ext_N_x_N_y phs_E(phis_n);

  Node_Data mob(Nx, Ny); //used to calculate the mobility parameter

  //a_11
  for(i = 0; i <= Nx + 1; i++)
    for(j = 0; j <= Ny; j++)
      {
        tmp = Lambda_c*phic_n[i][j]*(1.0 - phic_n[i][j]);
        mob[i][j] = (tmp >= 0) ? tmp : 0.0;
      }

  Node_Data_x_Ave_Ext_N a_11_x(mob);
  Node_Data_y_Ave_Ext_N a_11_y(mob);

  //a_22
  for(i = 0; i <= Nx + 1; i++)
    for(j = 0; j <= Ny; j++)
      {
        tmp = Lambda_b*phib_n[i][j]*(1.0 - phib_n[i][j]);
        mob[i][j] = (tmp >= 0) ? tmp : 0.0;
      }

  Node_Data_x_Ave_Ext_N a_22_x(mob);
  Node_Data_y_Ave_Ext_N a_22_y(mob);

  //a_33
  for(i = 0; i <= Nx + 1; i++)
    for(j = 0; j <= Ny; j++)
      {
        tmp = Lambda_s*phis_n[i][j]*(1.0 - phis_n[i][j]);
        mob[i][j] = (tmp >= 0) ? tmp : 0.0;
      }

  Node_Data_x_Ave_Ext_N a_33_x(mob);
  Node_Data_y_Ave_Ext_N a_33_y(mob);

  //a_12
  Node_Data_x_Ave_Ext_N a_12_x(Nx, Ny);
  Node_Data_y_Ave_Ext_N a_12_y(Nx, Ny);

  for(i = 0; i <= Nx + 2; i++)
    for(j = 0; j <= Ny; j++)
      a_12_x[i][j] = .5*(a_33_x[i][j] - a_11_x[i][j] - a_22_x[i][j]);

  for(i = 0; i <= Nx + 1; i++)
    for(j = 0; j <= Ny + 1; j++)
      a_12_y[i][j] = .5*(a_33_y[i][j] - a_11_y[i][j] - a_22_y[i][j]);

  //a_13
  Node_Data_x_Ave_Ext_N a_13_x(Nx, Ny);
  Node_Data_y_Ave_Ext_N a_13_y(Nx, Ny);

  for(i = 0; i <= Nx + 2; i++)
    for(j = 0; j <= Ny; j++)
      a_13_x[i][j] = .5*(a_22_x[i][j] - a_11_x[i][j] - a_33_x[i][j]);

  for(i = 0; i <= Nx + 1; i++)
    for(j = 0; j <= Ny + 1; j++)
      a_13_y[i][j] = .5*(a_22_y[i][j] - a_11_y[i][j] - a_33_y[i][j]);

//   //*************************************DEBUG*********************************************
//   sprintf(debfile,"%sa_13_x",SaveDir);
//   Write_Node_Data_x_Ave_Ext_N(debfile, a_13_x);
//   sprintf(debfile,"%sa_13_y",SaveDir);
//   Write_Node_Data_y_Ave_Ext_N(debfile, a_13_y);
//   //*************************************DEBUG*********************************************

  // vector of length 13 to put the coefficients from the \Delta\Delta phi
  valarray<double> DLAP(13);  //double Laplacian

  // vector of length 5 to put the coefficients fromt the div(phi*grad(phi)) term
  valarray<double> LAP(5);


  //***********************DEBUG***************************
  double cW[7] = {1.0, 1.0, 1.0, 1.0, .0, .0, .0};
 //  double cW[7] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  //***********************DEBUG***************************

  //Calculation of A[i][j] and  the RHS, it is same for any (i,j) since we use extended data structure containing
  //nodes outside the boundary 
  for(j = 0; j <= Ny; j++)
    {
      center = j*13 + 6; 
      for(i = 0; i <= Nx; i++)
	{
	  DLAP[0]  = Idy4*a_11_y[i][j];
	  DLAP[1]  = Idx2*Idy2*(a_11_x[i][j] + a_11_y[i][j]);
	  DLAP[2]  = - (Idx2*Idy2*(a_11_x[i+1][j] + a_11_x[i][j]) + Idy4 *a_11_y[i][j+1] + Idy2*(3.0*Idy2 + 2.0*Idx2)*a_11_y[i][j]);
	  DLAP[3]  = Idx2*Idy2*(a_11_x[i+1][j] + a_11_y[i][j]);
	  DLAP[4]  = Idx4*a_11_x[i][j];
	  DLAP[5]  = -(Idx4*a_11_x[i+1][j] + Idx2*(3.0*Idx2 + 2.0*Idy2)*a_11_x[i][j] + Idx2*Idy2*(a_11_y[i][j+1] + a_11_y[i][j]));
	  DLAP[6]  = (Idx2*(3.0*Idx2 + 2.0*Idy2)*(a_11_x[i+1][j] + a_11_x[i][j]) + 
                      Idy2*(3.0*Idy2 + 2.0*Idx2)*(a_11_y[i][j+1] + a_11_y[i][j]));
	  DLAP[7]  = -(Idx2*(3.0*Idx2 + 2.0*Idy2)*a_11_x[i+1][j] + Idx4*a_11_x[i][j] + Idx2*Idy2*(a_11_y[i][j+1] + a_11_y[i][j]));
	  DLAP[8]  = Idx4*a_11_x[i+1][j];
	  DLAP[9]  = Idx2*Idy2*(a_11_x[i][j] + a_11_y[i][j+1]);
	  DLAP[10] = -(Idx2*Idy2*(a_11_x[i+1][j] + a_11_x[i][j]) + Idy4*a_11_y[i][j] + Idy2*(2.0*Idx2 + 3.0*Idy2)*a_11_y[i][j+1]);
	  DLAP[11] = Idx2*Idy2*(a_11_x[i+1][j] + a_11_y[i][j+1]);
	  DLAP[12] = Idy4*a_11_y[i][j+1];

	  A[center][i]   = 1.0/dt + alpha*Gamma_0*DLAP[6];
	  A[center-1][i] = -.5*un[i][j]/dx + alpha*Gamma_0* DLAP[5];
	  A[center-2][i] = alpha*Gamma_0*DLAP[4];
	  A[center-3][i] = alpha*Gamma_0*DLAP[3];
	  A[center-4][i] = -.5*vn[i][j]/dy + alpha*Gamma_0*DLAP[2];
	  A[center-5][i] = alpha*Gamma_0*DLAP[1];
	  A[center-6][i] = alpha*Gamma_0*DLAP[0];
	  A[center+1][i] = .5*un[i][j]/dx + alpha*Gamma_0* DLAP[7];
	  A[center+2][i] = alpha*Gamma_0*DLAP[8];
	  A[center+3][i] = alpha*Gamma_0*DLAP[9];
	  A[center+4][i] = .5*vn[i][j]/dy + alpha*Gamma_0*DLAP[10];
	  A[center+5][i] = alpha*Gamma_0*DLAP[11];
	  A[center+6][i] = alpha*Gamma_0*DLAP[12];

	  curr = j*LX + i;
      
	  F[curr] = cW[0]*phic_n[i][j]/dt + cW[1]*g_c(Ca_vol_coef, k_p, S_n[i][j]);

	  tmp = (alpha - 1.0)*Gamma_0*(DLAP[0]*phc_E.get(i,j-2) +
                           DLAP[1]*phc_E.get(i-1,j-1) + DLAP[2]*phc_E.get(i,j-1)    + DLAP[3]*phc_E.get(i+1,j-1) +
			   DLAP[4]*phc_E.get(i-2,j)   + DLAP[5]*phc_E.get(i-1,j)    + DLAP[6]*phc_E.get(i,j)  + 
			   DLAP[7]*phc_E.get(i+1,j)   + DLAP[8]*phc_E.get(i+2,j)    + DLAP[9]*phc_E.get(i-1,j+1) +
			   DLAP[10]*phc_E.get(i,j+1)  + DLAP[11]*phc_E.get(i+1,j+1) + DLAP[12]*phc_E.get(i,j+2) ); 

	  F[curr] += cW[2]*tmp;  //the double Laplacian term due to phi_c

          LAP[0] = Idy2*a_11_y[i][j];
          LAP[1] = Idx2*a_11_x[i][j];
	  LAP[2] = -(Idx2*(a_11_x[i+1][j] + a_11_x[i][j]) + Idy2*(a_11_y[i][j+1] + a_11_y[i][j]));
          LAP[3] = Idx2*a_11_x[i+1][j];
          LAP[4] = Idy2*a_11_y[i][j+1];

	  tmp = LAP[0]*chem_c(phc_E.get(i,j-1), phb_E.get(i,j-1), phs_E.get(i,j-1)) + LAP[1]*chem_c(phc_E.get(i-1,j), phb_E.get(i-1,j), phs_E.get(i-1,j)) + 
                LAP[2]*chem_c(phc_E.get(i,j), phb_E.get(i,j), phs_E.get(i,j)) +
                LAP[3]*chem_c(phc_E.get(i+1,j), phb_E.get(i+1,j), phs_E.get(i+1,j)) + LAP[4]*chem_c(phc_E.get(i,j+1), phb_E.get(i,j+1), phs_E.get(i,j+1));

	  F[curr] += cW[3]*tmp;  //the Laplacian term due to phi_c

          LAP[0] = Idy2*a_13_y[i][j];
          LAP[1] = Idx2*a_13_x[i][j];
	  LAP[2] = -(Idx2*(a_13_x[i+1][j] + a_13_x[i][j]) + Idy2*(a_13_y[i][j+1] + a_13_y[i][j]));
          LAP[3] = Idx2*a_13_x[i+1][j];
          LAP[4] = Idy2*a_13_y[i][j+1];
 
	  tmp = LAP[0]*chem_s(phc_E.get(i,j-1), phb_E.get(i,j-1), phs_E.get(i,j-1)) + LAP[1]*chem_s(phc_E.get(i-1,j), phb_E.get(i-1,j), phs_E.get(i-1,j)) + 
                LAP[2]*chem_s(phc_E.get(i,j), phb_E.get(i,j), phs_E.get(i,j)) +
                LAP[3]*chem_s(phc_E.get(i+1,j), phb_E.get(i+1,j), phs_E.get(i+1,j)) + LAP[4]*chem_s(phc_E.get(i,j+1), phb_E.get(i,j+1), phs_E.get(i,j+1)); 

	  F[curr] += cW[4]*tmp;  //the Laplacian term due to phi_s

          LAP[0] = Idy2*a_12_y[i][j];
          LAP[1] = Idx2*a_12_x[i][j];
	  LAP[2] = -(Idx2*(a_12_x[i+1][j] + a_12_x[i][j]) + Idy2*(a_12_y[i][j+1] + a_12_y[i][j]));
          LAP[3] = Idx2*a_12_x[i+1][j];
          LAP[4] = Idy2*a_12_y[i][j+1];
 
	  tmp = LAP[0]*chem_b(phc_E.get(i,j-1), phb_E.get(i,j-1), phs_E.get(i,j-1)) + LAP[1]*chem_b(phc_E.get(i-1,j), phb_E.get(i-1,j), phs_E.get(i-1,j)) + 
                LAP[2]*chem_b(phc_E.get(i,j), phb_E.get(i,j), phs_E.get(i,j)) +
                LAP[3]*chem_b(phc_E.get(i+1,j), phb_E.get(i+1,j), phs_E.get(i+1,j)) + LAP[4]*chem_b(phc_E.get(i,j+1), phb_E.get(i,j+1), phs_E.get(i,j+1)); 

	  F[curr] += cW[5]*tmp;  //the Laplacian term due to phi_b

	  DLAP[0]  = Idy4*a_12_y[i][j];
	  DLAP[1]  = Idx2*Idy2*(a_12_x[i][j] + a_12_y[i][j]);
	  DLAP[2]  = - (Idx2*Idy2*(a_12_x[i+1][j] + a_12_x[i][j]) + Idy4 *a_12_y[i][j+1] + Idy2*(3.0*Idy2 + 2.0*Idx2)*a_12_y[i][j]);
	  DLAP[3]  = Idx2*Idy2*(a_12_x[i+1][j] + a_12_y[i][j]);
	  DLAP[4]  = Idx4*a_12_x[i][j];
	  DLAP[5]  = -(Idx4*a_12_x[i+1][j] + Idx2*(3.0*Idx2 + 2.0*Idy2)*a_12_x[i][j] + Idx2*Idy2*(a_12_y[i][j+1] + a_12_y[i][j]));
	  DLAP[6]  = (Idx2*(3.0*Idx2 + 2.0*Idy2)*(a_12_x[i+1][j] + a_12_x[i][j]) + 
                      Idy2*(3.0*Idy2 + 2.0*Idx2)*(a_12_y[i][j+1] + a_12_y[i][j]));
	  DLAP[7]  = -(Idx2*(3.0*Idx2 + 2.0*Idy2)*a_12_x[i+1][j] + Idx4*a_12_x[i][j] + Idx2*Idy2*(a_12_y[i][j+1] + a_12_y[i][j]));
	  DLAP[8]  = Idx4*a_12_x[i+1][j];
	  DLAP[9]  = Idx2*Idy2*(a_12_x[i][j] + a_12_y[i][j+1]);
	  DLAP[10] = -(Idx2*Idy2*(a_12_x[i+1][j] + a_12_x[i][j]) + Idy4*a_12_y[i][j] + Idy2*(2.0*Idx2 + 3.0*Idy2)*a_12_y[i][j+1]);
	  DLAP[11] = Idx2*Idy2*(a_12_x[i+1][j] + a_12_y[i][j+1]);
	  DLAP[12] = Idy4*a_12_y[i][j+1];

	  tmp = - Gamma_1*(DLAP[0]*phb_E.get(i,j-2) +
                           DLAP[1]*phb_E.get(i-1,j-1) + DLAP[2]*phb_E.get(i,j-1)    + DLAP[3]*phb_E.get(i+1,j-1) +
			   DLAP[4]*phb_E.get(i-2,j)   + DLAP[5]*phb_E.get(i-1,j)    + DLAP[6]*phb_E.get(i,j)  + 
			   DLAP[7]*phb_E.get(i+1,j)   + DLAP[8]*phb_E.get(i+2,j)    + DLAP[9]*phb_E.get(i-1,j+1) +
			   DLAP[10]*phb_E.get(i,j+1)  + DLAP[11]*phb_E.get(i+1,j+1) + DLAP[12]*phb_E.get(i,j+2) ); 

	  F[curr] += cW[6]*tmp; //the double Lapacian term due to phi_b
	}
    }

  //now enforce the Boundary conditions, no flux at x = 0, y = 0, y = 1; Dirichlet at x = L

  j = 0;
  center = j*13 + 6; 

  i = 0;

  A[center+6][i] += A[center-6][i];
  A[center+4][i] += A[center-4][i];      
  A[center+5][i] += A[center-3][i] +  A[center+3][i] + A[center-5][i];    
  A[center+1][i] += A[center-1][i];
  A[center+2][i] += A[center-2][i];   

  i = 1;

  A[center+6][i] += A[center-6][i];
  A[center+3][i] += A[center-5][i];      
  A[center+4][i] += A[center-4][i];      
  A[center+5][i] += A[center-3][i]; 
  A[center][i]   += A[center-2][i];   

  for(i = 2; i <= Nx-2; i++)
    {
      A[center+6][i] += A[center-6][i];
      A[center+3][i] += A[center-5][i];      
      A[center+4][i] += A[center-4][i];      
      A[center+5][i] += A[center-3][i]; 
    } 

  i = Nx-1;

  A[center+6][i] += A[center-6][i];
  A[center+3][i] += A[center-5][i];      
  A[center+4][i] += A[center-4][i];      
  A[center+5][i] += A[center-3][i]; 

  curr = j*LX + i;

  F[curr] -= A[center+2][i]*Rph_np1[j];

  i = Nx;

  A[center+6][i] += A[center-6][i];
  A[center+4][i] += A[center-4][i];      
  A[center+3][i] += A[center-5][i];
  A[center][i] -= A[center+2][i];
 
  curr = j*LX + i;

  F[curr] -= A[center+2][i]*Rph_cor[j] + A[center+1][i]*Rph_np1[j] + 
             A[center-3][i]*Rph_np1[j+1] + A[center+5][i]*Rph_np1[j+1];

  j = 1; 
  center = j*13 + 6; 

  i = 0;

  A[center+1][i] += A[center-1][i];
  A[center+2][i] += A[center-2][i];      
  A[center+5][i] += A[center+3][i];      
  A[center-3][i] += A[center-5][i];      
  A[center][i]   += A[center-6][i];

  i = 1;

  A[center][i] += A[center-6][i] + A[center-2][i];

  for(i = 2; i <= Nx-2; i++)
    {
      A[center][i] += A[center-6][i];
    } 

  i = Nx - 1;

  A[center][i] += A[center-6][i]; 

  curr = j*LX + i;
  F[curr] -= A[center+2][i]*Rph_np1[j];

  i = Nx;

  A[center][i] += A[center-6][i] - A[center+2][i];

  curr = j*LX + i;
  F[curr] -= A[center+1][i]*Rph_np1[j] + A[center+2][i]*Rph_cor[j] + A[center-3][i]*Rph_np1[j-1] + A[center+5][i]*Rph_np1[j+1];


  for(j = 2; j <= Ny-2; j++)
    {
      center = j*13 + 6; 

      i = 0;

      A[center+1][i] += A[center-1][i];
      A[center+2][i] += A[center-2][i];      
      A[center+5][i] += A[center+3][i];      
      A[center-3][i] += A[center-5][i];      

      i = 1;

      A[center][i] += A[center-2][i];

      i = Nx - 1;
      
      curr = j*LX + i;
      F[curr] -= A[center+2][i]*Rph_np1[j];

      i = Nx;

      A[center][i] -= A[center+2][i];

      curr = j*LX + i;
      F[curr] -= A[center+1][i]*Rph_np1[j] + A[center+2][i]*Rph_cor[j] + A[center-3][i]*Rph_np1[j-1] + A[center+5][i]*Rph_np1[j+1];
    }

  j = Ny-1 ;
  center = j*13 + 6; 

  i = 0;

  A[center+1][i] += A[center-1][i];
  A[center+2][i] += A[center-2][i];      
  A[center+5][i] += A[center+3][i];      
  A[center-3][i] += A[center-5][i];      
  A[center][i]   += A[center+6][i];

  i = 1;
  A[center][i] += A[center+6][i] + A[center-2][i];

  for(i = 2; i <= Nx-2; i++)
    {
      A[center][i] += A[center+6][i];
    }

  i = Nx - 1;

  A[center][i] += A[center+6][i];

  curr = j*LX + i;
  F[curr] -= A[center+2][i]*Rph_np1[j];

  i = Nx;

  A[center][i] += A[center+6][i] - A[center+2][i];

  curr = j*LX + i;
  F[curr] -= A[center+1][i]*Rph_np1[j] + A[center+2][i]*Rph_cor[j] + A[center+5][i]*Rph_np1[j+1] + A[center-3][i]*Rph_np1[j-1];

  j = Ny;
  center = j*13 + 6; 

  i = 0;

  A[center-6][i] += A[center+6][i];
  A[center-4][i] += A[center+4][i];      
  A[center-3][i] += A[center+5][i] +  A[center+3][i] + A[center-5][i];    
  A[center+1][i] += A[center-1][i];
  A[center+2][i] += A[center-2][i];   

  i = 1;

  A[center-6][i] += A[center+6][i];
  A[center-5][i] += A[center+3][i];      
  A[center-4][i] += A[center+4][i];      
  A[center-3][i] += A[center+5][i];
  A[center][i]   += A[center-2][i];   

  for(i = 2; i <= Nx-2; i++)
    {
      A[center-6][i] += A[center+6][i];
      A[center-5][i] += A[center+3][i];      
      A[center-4][i] += A[center+4][i];      
      A[center-3][i] += A[center+5][i]; 
    } 

  i = Nx - 1;

  A[center-6][i] += A[center+6][i];
  A[center-5][i] += A[center+3][i];      
  A[center-4][i] += A[center+4][i];      
  A[center-3][i] += A[center+5][i];

  curr = j*LX + i;
  F[curr] -= A[center+2][i]*Rph_np1[j];
 

  i = Nx;

  A[center-6][i] += A[center+6][i];
  A[center-4][i] += A[center+4][i];      
  A[center-5][i] += A[center+3][i];
  A[center][i] -= A[center+2][i];

  curr = j*LX + i;
  F[curr] -= A[center+1][i]*Rph_np1[j] + A[center+2][i]*Rph_cor[j] + A[center+5][i]*Rph_np1[j-1] + A[center-3][i]*Rph_np1[j-1];
}


/*

//************************************************************************************************************************
void Build_c_F_x_N_y(c_Mat& A, valarray<double>& F, const valarray<double>& Lc, const valarray<double>& Rc, const Node_Data& phib_n, const Node_Data& phis_np1_sin, const Node_Data& phis_n_sin, 
                       const Node_Data& c_n, const Node_Data& un, const Node_Data& vn, double dx, double dy, double dt,
		     const Control_Parameter& CP){

  int i, j, Nx, Ny, center, curr, LX;

  Nx = A.get_Nx();
  Ny = A.get_Ny();
  LX = A.get_LX();

  if(A.get_dim() != F.size() || Lc.size() != Ny+1 || Rc.size() != Ny + 1)
    {
      cout << endl << "Dimension incompatible in Build_c_Flow_x_N_y(), ERROR at line " << __LINE__ << " of file " 
           << __FILE__ << " !! " << endl;

      exit(-1);
    }

  double dx2, dy2, tmp, alpha;

  dx2 = dx*dx;
  dy2 = dy*dy;

  alpha = 0.5;

  Node_Data phis_np1(phis_np1_sin);
  Node_Data phis_n(phis_n_sin);

  Confine_Node_Data(phis_np1, PHIS_reg, 1.0);
  Confine_Node_Data(phis_n, PHIS_reg, 1.0);

  Node_Data_Ext_N_x_N_y phis_np1_e(phis_np1);
  Node_Data_Ext_N_x_N_y phis_n_e(phis_n);
  Node_Data_Ext_N_x_N_y c_n_e(c_n);

  Node_Data_x_Ave_Ext_N phis_np1_x(phis_np1);
  Node_Data_x_Ave_Ext_N phis_n_x(phis_n); 

  Node_Data_y_Ave_Ext_N phis_np1_y(phis_np1);
  Node_Data_y_Ave_Ext_N phis_n_y(phis_n); 
 

  //build the matrix and rhs
  
  for(j = 0; j <= Ny; j++)
    {
      center = j*5 + 2;

      for(i = 1; i <= Nx; i++)
	{
	  A[center][i-1] = phis_np1[i][j]/dt +
	                alpha*Ds*(phis_np1_x[i+1][j] + phis_np1_x[i][j])/dx2 +
	                alpha*Ds*(phis_np1_y[i][j+1] + phis_np1_y[i][j])/dy2;

          A[center-1][i-1] =  -(.5)*un[i][j]*phis_np1_e.get(i-1,j)/dx - alpha*Ds*phis_np1_x[i][j]/dx2;
	  A[center+1][i-1] =   (.5)*un[i][j]*phis_np1_e.get(i+1,j)/dx - alpha*Ds*phis_np1_x[i+1][j]/dx2;
          A[center-2][i-1] =  -(.5)*vn[i][j]*phis_np1_e.get(i,j-1)/dy - alpha*Ds*phis_np1_y[i][j]/dy2;
          A[center+2][i-1] =   (.5)*vn[i][j]*phis_np1_e.get(i,j+1)/dy - alpha*Ds*phis_np1_y[i][j+1]/dy2;

          //calculate the RHS
          curr = j*LX + i - 1;

	  F[curr] = phis_n[i][j]*c_n[i][j]/dt - g_s(phib_n[i][j],c_n[i][j]);
	
	  //add the explict part from the diffusion term from Crank-Nicolson scheme
	  tmp = ((1.0 - alpha)*Ds/dx2)*(phis_n_x[i+1][j]*(c_n_e.get(i+1,j) - c_n_e.get(i,j)) - 
			       phis_n_x[i][j]*(c_n_e.get(i,j) - c_n_e.get(i-1,j))) +
	        ((1.0 - alpha)*Ds/dy2)*(phis_n_y[i][j+1]*(c_n_e.get(i,j+1) - c_n_e.get(i,j)) - 
                             phis_n_y[i][j]*(c_n_e.get(i,j) - c_n_e.get(i,j-1)));
						
	  F[curr] += tmp;
	}				     
					
    }

  //Enfore the no-flux BC at j = 0 and j = Ny,
  //add the coef of j = -1, to coef of j = 1, j = Ny + 1 to j = Ny - 1

  j = 0;
  center = j*5 + 2;

  for(i = 1; i <= Nx; i++)
    {
      A[center+2][i-1] += A[center-2][i-1];
      A[center-2][i-1] = 0.0;
    }

  j = Ny;
  center = j*5 + 2;

  for(i = 1; i <= Nx; i++)
    {
      A[center-2][i-1] += A[center+2][i-1];
      A[center+2][i-1] = 0.0;
    }

  // enforce  Dirichelet BC at x = 0 and x = L by substrating the known terms

  for(j = 0; j <= Ny; j++)
    {
      center = j*5 + 2;

      i = 1;
 
      curr = j*LX + i - 1;

      F[curr] -= A[center-1][i-1]*Lc[j]; //(CP.Ini_c + CP.Ini_ce*sin(CP.Ini_c_omega*t_current)); 

      i = Nx;

      curr = j*LX + i - 1;

      F[curr] -= A[center+1][i-1]*Rc[j];
    }

}

//************************************************************************************************************************
void Build_Ur_F_x_N_y(c_Mat& A, valarray<double>& F, const valarray<double>& Lc, const valarray<double>& Rc, const Node_Data& phib_n, const Node_Data& phis_np1_sin, const Node_Data& phis_n_sin, 
		      const Node_Data& Ur_n, const Node_Data& un, const Node_Data& vn, double kur, double dx, double dy, double dt){

  int i, j, Nx, Ny, center, curr, LX;

  Nx = A.get_Nx();
  Ny = A.get_Ny();
  LX = A.get_LX();

  if(A.get_dim() != F.size() || Lc.size() != Ny+1 || Rc.size() != Ny + 1)
    {
      cout << endl << "Dimension incompatible in Build_c_Flow_x_N_y(), ERROR at line " << __LINE__ << " of file " 
           << __FILE__ << " !! " << endl;

      exit(-1);
    }
 
  double Idx2, Idy2, tmp;

  Idx2 = 1.0/dx/dx;
  Idy2 = 1.0/dy/dy;

  Node_Data phis_np1(phis_np1_sin);
  Node_Data phis_n(phis_n_sin);

  Confine_Node_Data(phis_np1, PHIS_reg, 1.0);
  Confine_Node_Data(phis_n, PHIS_reg, 1.0);


  Node_Data_Ext_N_x_N_y Ur_n_E(Ur_n);
  Node_Data_Ext_N_x_N_y phs_np1_E(phis_np1);
          
  Node_Data_x_Ave_Ext_N phs_n_x(phis_n);
  Node_Data_y_Ave_Ext_N phs_n_y(phis_n);
    
  Node_Data_x_Ave_Ext_N phs_np1_x(phis_np1);
  Node_Data_y_Ave_Ext_N phs_np1_y(phis_np1);
    
  for(j = 0; j <= Ny; j++)
      {
          center = j*5 + 2;
          for(i = 1; i <= Nx ; i++)
              {
                  A[center][i-1] = phis_np1[i][j]/dt +
                                .5*D_ur*((phs_np1_x[i][j] + phs_np1_x[i+1][j])*Idx2 +
                                         (phs_np1_y[i][j] + phs_np1_y[i][j+1])*Idy2);

                  A[center-1][i-1] =  -.5*un[i][j]*phs_np1_E.get(i-1,j)/dx - .5*D_ur*phs_np1_x[i][j]*Idx2;
                  A[center+1][i-1] =   .5*un[i][j]*phs_np1_E.get(i+1,j)/dx - .5*D_ur*phs_np1_x[i+1][j]*Idx2;
                  A[center-2][i-1] =  -.5*vn[i][j]*phs_np1_E.get(i,j-1)/dy - .5*D_ur*phs_np1_y[i][j]*Idy2;
                  A[center+2][i-1] =   .5*vn[i][j]*phs_np1_E.get(i,j+1)/dy - .5*D_ur*phs_np1_y[i][j+1]*Idy2;

                  //calculate the RHS
                  curr = j*LX + i - 1;

                  F[curr] = phis_n[i][j]*Ur_n[i][j]/dt + phis_n[i][j]*hydrolysis(kur, phib_n[i][j], Ur_n[i][j]);
 
                  //add the explict part from the diffusion term from Crank-Nicolson scheme
                  tmp = .5*D_ur*(Idx2*(phs_n_x[i+1][j]*(Ur_n_E.get(i+1,j) - Ur_n_E.get(i,j)) - 
                                       phs_n_x[i][j]*(Ur_n_E.get(i,j) - Ur_n_E.get(i-1,j))) +
                                 Idy2*(phs_n_y[i][j+1]*(Ur_n_E.get(i,j+1) - Ur_n_E.get(i,j)) - 
                                       phs_n_y[i][j]*(Ur_n_E.get(i,j) - Ur_n_E.get(i,j-1))));

                  F[curr] += tmp;
              }

      }

  //Enfore the no-flux BC at j = 0 and j = Ny,
  //add the coef of j = -1, to coef of j = 1, j = Ny + 1 to j = Ny - 1

  j = 0;
  center = j*5 + 2;

  for(i = 1; i <= Nx; i++)
    {
      A[center+2][i-1] += A[center-2][i-1];
      A[center-2][i-1] = 0.0;
    }

  j = Ny;
  center = j*5 + 2;

  for(i = 1; i <= Nx; i++)
    {
      A[center-2][i-1] += A[center+2][i-1];
      A[center+2][i-1] = 0.0;
    }

  // enforce  Dirichelet BC at x = 0 and x = L by substrating the known terms

  for(j = 0; j <= Ny; j++)
    {
      center = j*5 + 2;

      i = 1;
 
      curr = j*LX + i - 1;

      F[curr] -= A[center-1][i-1]*Lc[j]; //(CP.Ini_c + CP.Ini_ce*sin(CP.Ini_c_omega*t_current)); 

      i = Nx;

      curr = j*LX + i - 1;

      F[curr] -= A[center+1][i-1]*Rc[j];
    }

}

//************************************************************************************************************************
void Build_Ca_F_x_N_y(c_Mat& A, valarray<double>& F, const valarray<double>& Lc, const valarray<double>& Rc, const Node_Data& S_n, const Node_Data& phis_np1_sin, const Node_Data& phis_n_sin, 
		      const Node_Data& Ca_n, const Node_Data& un, const Node_Data& vn, double kp, double dx, double dy, double dt){

  int i, j, Nx, Ny, center, curr, LX;

  Nx = A.get_Nx();
  Ny = A.get_Ny();
  LX = A.get_LX();

  if(A.get_dim() != F.size() || Lc.size() != Ny+1 || Rc.size() != Ny + 1)
    {
      cout << endl << "Dimension incompatible in Build_c_Flow_x_N_y(), ERROR at line " << __LINE__ << " of file " 
           << __FILE__ << " !! " << endl;

      exit(-1);
    }

  double Idx2, Idy2, tmp;

  Idx2 = 1.0/dx/dx;
  Idy2 = 1.0/dy/dy;

  Node_Data phis_np1(phis_np1_sin);
  Node_Data phis_n(phis_n_sin);

  Confine_Node_Data(phis_np1, PHIS_reg, 1.0);
  Confine_Node_Data(phis_n, PHIS_reg, 1.0);

  Node_Data_Ext_N_x_N_y Ca_n_E(Ca_n);
  Node_Data_Ext_N_x_N_y phs_np1_E(phis_np1);
          
  Node_Data_x_Ave_Ext_N phs_n_x(phis_n);
  Node_Data_y_Ave_Ext_N phs_n_y(phis_n);
    
  Node_Data_x_Ave_Ext_N phs_np1_x(phis_np1);
  Node_Data_y_Ave_Ext_N phs_np1_y(phis_np1);
    
  for(j = 0; j <= Ny; j++)
      {
          center = j*5 + 2;
          for(i = 1; i <= Nx; i++)
              {
                  A[center][i-1] = phis_np1[i][j]/dt +
                                .5*D_ca*((phs_np1_x[i][j] + phs_np1_x[i+1][j])*Idx2 +
                                         (phs_np1_y[i][j] + phs_np1_y[i][j+1])*Idy2);

                  A[center-1][i-1] =  -.5*un[i][j]*phs_np1_E.get(i-1,j)/dx - .5*D_ca*phs_np1_x[i][j]*Idx2;
                  A[center+1][i-1] =   .5*un[i][j]*phs_np1_E.get(i+1,j)/dx - .5*D_ca*phs_np1_x[i+1][j]*Idx2;
                  A[center-2][i-1] =  -.5*vn[i][j]*phs_np1_E.get(i,j-1)/dy - .5*D_ca*phs_np1_y[i][j]*Idy2;
                  A[center+2][i-1] =   .5*vn[i][j]*phs_np1_E.get(i,j+1)/dy - .5*D_ca*phs_np1_y[i][j+1]*Idy2;

                  //calculate the RHS
                  curr = j*LX + i - 1;

                  F[curr] = phis_n[i][j]*Ca_n[i][j]/dt + phis_n[i][j]*precipitate(kp, S_n[i][j]);
 
                  //add the explict part from the diffusion term from Crank-Nicolson scheme
                  tmp = .5*D_ca*(Idx2*(phs_n_x[i+1][j]*(Ca_n_E.get(i+1,j) - Ca_n_E.get(i,j)) - 
                                       phs_n_x[i][j]*(Ca_n_E.get(i,j) - Ca_n_E.get(i-1,j))) +
                                 Idy2*(phs_n_y[i][j+1]*(Ca_n_E.get(i,j+1) - Ca_n_E.get(i,j)) - 
                                       phs_n_y[i][j]*(Ca_n_E.get(i,j) - Ca_n_E.get(i,j-1)))); 

                  F[curr] += tmp;
              }

      }


  //Enfore the no-flux BC at j = 0 and j = Ny,
  //add the coef of j = -1, to coef of j = 1, j = Ny + 1 to j = Ny - 1

  j = 0;
  center = j*5 + 2;

  for(i = 1; i <= Nx; i++)
    {
      A[center+2][i-1] += A[center-2][i-1];
      A[center-2][i-1] = 0.0;
    }

  j = Ny;
  center = j*5 + 2;

  for(i = 1; i <= Nx; i++)
    {
      A[center-2][i-1] += A[center+2][i-1];
      A[center+2][i-1] = 0.0;
    }

  // enforce  Dirichelet BC at x = 0 and x = L by substrating the known terms

  for(j = 0; j <= Ny; j++)
    {
      center = j*5 + 2;

      i = 1;
 
      curr = j*LX + i - 1;

      F[curr] -= A[center-1][i-1]*Lc[j]; //(CP.Ini_c + CP.Ini_ce*sin(CP.Ini_c_omega*t_current)); 

      i = Nx;

      curr = j*LX + i - 1;

      F[curr] -= A[center+1][i-1]*Rc[j];
    }

}

//************************************************************************************************************************
void Build_CT_F_x_N_y(c_Mat& A, valarray<double>& F, const valarray<double>& Lc, const valarray<double>& Rc, const Node_Data& CT_rate_n, const Node_Data& phis_np1_sin, const Node_Data& phis_n_sin, 
		      const Node_Data& CT_n, const Node_Data& un, const Node_Data& vn, double dx, double dy, double dt){

  int i, j, Nx, Ny, center, curr, LX;

  Nx = A.get_Nx();
  Ny = A.get_Ny();
  LX = A.get_LX();

  if(A.get_dim() != F.size() || Lc.size() != Ny+1 || Rc.size() != Ny + 1)
    {
      cout << endl << "Dimension incompatible in Build_c_Flow_x_N_y(), ERROR at line " << __LINE__ << " of file " 
           << __FILE__ << " !! " << endl;

      exit(-1);
    }

  double Idx2, Idy2, tmp;

  Idx2 = 1.0/dx/dx;
  Idy2 = 1.0/dy/dy;

  Node_Data phis_np1(phis_np1_sin);
  Node_Data phis_n(phis_n_sin);

  Confine_Node_Data(phis_np1, PHIS_reg, 1.0);
  Confine_Node_Data(phis_n, PHIS_reg, 1.0);


  Node_Data_Ext_N_x_N_y CT_n_E(CT_n);
  Node_Data_Ext_N_x_N_y phs_np1_E(phis_np1);
          
  Node_Data_x_Ave_Ext_N phs_n_x(phis_n);
  Node_Data_y_Ave_Ext_N phs_n_y(phis_n);
    
  Node_Data_x_Ave_Ext_N phs_np1_x(phis_np1);
  Node_Data_y_Ave_Ext_N phs_np1_y(phis_np1);
    
  for(j = 0; j <= Ny; j++)
      {
          center = j*5 + 2;
          for(i = 1; i <= Nx; i++)
              {
                  A[center][i-1] = phis_np1[i][j]/dt +
                                .5*D_ct*((phs_np1_x[i][j] + phs_np1_x[i+1][j])*Idx2 +
                                         (phs_np1_y[i][j] + phs_np1_y[i][j+1])*Idy2);

                  A[center-1][i-1] =  -.5*un[i][j]*phs_np1_E.get(i-1,j)/dx - .5*D_ct*phs_np1_x[i][j]*Idx2;
                  A[center+1][i-1] =   .5*un[i][j]*phs_np1_E.get(i+1,j)/dx - .5*D_ct*phs_np1_x[i+1][j]*Idx2;
                  A[center-2][i-1] =  -.5*vn[i][j]*phs_np1_E.get(i,j-1)/dy - .5*D_ct*phs_np1_y[i][j]*Idy2;
                  A[center+2][i-1] =   .5*vn[i][j]*phs_np1_E.get(i,j+1)/dy - .5*D_ct*phs_np1_y[i][j+1]*Idy2;

                  //calculate the RHS
                  curr = j*LX + i - 1;

                  F[curr] = phis_n[i][j]*CT_n[i][j]/dt + phis_n[i][j]*CT_rate_n[i][j];
 
                  //add the explict part from the diffusion term from Crank-Nicolson scheme
                  tmp = .5*D_ct*(Idx2*(phs_n_x[i+1][j]*(CT_n_E.get(i+1,j) - CT_n_E.get(i,j)) - 
                                       phs_n_x[i][j]*(CT_n_E.get(i,j) - CT_n_E.get(i-1,j))) +
                                 Idy2*(phs_n_y[i][j+1]*(CT_n_E.get(i,j+1) - CT_n_E.get(i,j)) - 
                                       phs_n_y[i][j]*(CT_n_E.get(i,j) - CT_n_E.get(i,j-1))));
                  
                  F[curr] += tmp;
              }

      }


  //Enfore the no-flux BC at j = 0 and j = Ny,
  //add the coef of j = -1, to coef of j = 1, j = Ny + 1 to j = Ny - 1

  j = 0;
  center = j*5 + 2;

  for(i = 1; i <= Nx; i++)
    {
      A[center+2][i-1] += A[center-2][i-1];
      A[center-2][i-1] = 0.0;
    }

  j = Ny;
  center = j*5 + 2;

  for(i = 1; i <= Nx; i++)
    {
      A[center-2][i-1] += A[center+2][i-1];
      A[center+2][i-1] = 0.0;
    }

  // enforce  Dirichelet BC at x = 0 and x = L by substrating the known terms

  for(j = 0; j <= Ny; j++)
    {
      center = j*5 + 2;

      i = 1;
 
      curr = j*LX + i - 1;

      F[curr] -= A[center-1][i-1]*Lc[j]; //(CP.Ini_c + CP.Ini_ce*sin(CP.Ini_c_omega*t_current)); 

      i = Nx;

      curr = j*LX + i - 1;

      F[curr] -= A[center+1][i-1]*Rc[j];
    }

}

//************************************************************************************************************************
void Build_OH_F_x_N_y(c_Mat& A, valarray<double>& F, const valarray<double>& Lc, const valarray<double>& Rc, const Node_Data& OH_rate_n, const Node_Data& phis_np1_sin, const Node_Data& phis_n_sin, 
		      const Node_Data& OH_n, const Node_Data& un, const Node_Data& vn, double dx, double dy, double dt){

  int i, j, Nx, Ny, center, curr, LX;

  Nx = A.get_Nx();
  Ny = A.get_Ny();
  LX = A.get_LX();

  if(A.get_dim() != F.size() || Lc.size() != Ny+1 || Rc.size() != Ny + 1)
    {
      cout << endl << "Dimension incompatible in Build_c_Flow_x_N_y(), ERROR at line " << __LINE__ << " of file " 
           << __FILE__ << " !! " << endl;

      exit(-1);
    }

 
  double Idx2, Idy2, tmp;

  Idx2 = 1.0/dx/dx;
  Idy2 = 1.0/dy/dy;

  Node_Data phis_np1(phis_np1_sin);
  Node_Data phis_n(phis_n_sin);

  Confine_Node_Data(phis_np1, PHIS_reg, 1.0);
  Confine_Node_Data(phis_n, PHIS_reg, 1.0);

  Node_Data_Ext_N_x_N_y OH_n_E(OH_n);
  Node_Data_Ext_N_x_N_y phs_np1_E(phis_np1);
          
  Node_Data_x_Ave_Ext_N phs_n_x(phis_n);
  Node_Data_y_Ave_Ext_N phs_n_y(phis_n);
    
  Node_Data_x_Ave_Ext_N phs_np1_x(phis_np1);
  Node_Data_y_Ave_Ext_N phs_np1_y(phis_np1);
    
  for(j = 0; j <= Ny; j++)
      {
          center = j*5 + 2;
          for(i = 1; i <= Nx; i++)
              {
                  A[center][i-1] = phis_np1[i][j]/dt +
                                .5*D_oh*((phs_np1_x[i][j] + phs_np1_x[i+1][j])*Idx2 +
                                         (phs_np1_y[i][j] + phs_np1_y[i][j+1])*Idy2);

                  A[center-1][i-1] =  -.5*un[i][j]*phs_np1_E.get(i-1,j)/dx - .5*D_oh*phs_np1_x[i][j]*Idx2;
                  A[center+1][i-1] =   .5*un[i][j]*phs_np1_E.get(i+1,j)/dx - .5*D_oh*phs_np1_x[i+1][j]*Idx2;
                  A[center-2][i-1] =  -.5*vn[i][j]*phs_np1_E.get(i,j-1)/dy - .5*D_oh*phs_np1_y[i][j]*Idy2;
                  A[center+2][i-1] =   .5*vn[i][j]*phs_np1_E.get(i,j+1)/dy - .5*D_oh*phs_np1_y[i][j+1]*Idy2;

                  //calculate the RHS
                  curr = j*LX + i - 1;

                  F[curr] = phis_n[i][j]*OH_n[i][j]/dt + phis_n[i][j]*OH_rate_n[i][j];
 
                  //add the explict part from the diffusion term from Crank-Nicolson scheme
                  tmp = .5*D_oh*(Idx2*(phs_n_x[i+1][j]*(OH_n_E.get(i+1,j) - OH_n_E.get(i,j)) - 
                                       phs_n_x[i][j]*(OH_n_E.get(i,j) - OH_n_E.get(i-1,j))) +
                                 Idy2*(phs_n_y[i][j+1]*(OH_n_E.get(i,j+1) - OH_n_E.get(i,j)) - 
                                       phs_n_y[i][j]*(OH_n_E.get(i,j) - OH_n_E.get(i,j-1))));

                  F[curr] += tmp;
              }

      }

  //Enfore the no-flux BC at j = 0 and j = Ny,
  //add the coef of j = -1, to coef of j = 1, j = Ny + 1 to j = Ny - 1

  j = 0;
  center = j*5 + 2;

  for(i = 1; i <= Nx; i++)
    {
      A[center+2][i-1] += A[center-2][i-1];
      A[center-2][i-1] = 0.0;
    }

  j = Ny;
  center = j*5 + 2;

  for(i = 1; i <= Nx; i++)
    {
      A[center-2][i-1] += A[center+2][i-1];
      A[center+2][i-1] = 0.0;
    }

  // enforce  Dirichelet BC at x = 0 and x = L by substrating the known terms

  for(j = 0; j <= Ny; j++)
    {
      center = j*5 + 2;

      i = 1;
 
      curr = j*LX + i - 1;

      F[curr] -= A[center-1][i-1]*Lc[j]; //(CP.Ini_c + CP.Ini_ce*sin(CP.Ini_c_omega*t_current)); 

      i = Nx;

      curr = j*LX + i - 1;

      F[curr] -= A[center+1][i-1]*Rc[j];
    }
 
}

*/ 

                     
//************************************************************************************************************************

void Build_Slow_Process(c_Mat& A, valarray<double>& F, int k, const Node_Data& D_n, const Node_Data& D_np1, double Z, const Node_Data& C_n, 
                        const valarray<double>& Lc, const valarray<double>& Rc, const Node_Data& phi_b, const Node_Data& Ur, const Node_Data& S, 
                        const Node_Data& c, const Node_Data& phis_np1_sin, const Node_Data& phis_n_sin, const Node_Data& un, const Node_Data& vn, 
                        const Node_Data& EP, double a, double dx, double dy, double dt, const Control_Parameter& CP){


  int i, j, Nx, Ny, center, curr, LX;

  Nx = A.get_Nx();
  Ny = A.get_Ny();
  LX = A.get_LX();

  if(A.get_dim() != F.size() || Lc.size() != Ny+1 || Rc.size() != Ny + 1)
    {
      cout << endl << "Dimension incompatible in Build_Slow_Process(), ERROR at line " << __LINE__ << " of file " 
           << __FILE__ << " !! " << endl;

      exit(-1);
    }

 
  double Idx2, Idy2, tmp1, tmp2, tmp3, tmp4, EP_ind, alpha;

  alpha = CP.alph;

  //********************DEBUG PURPOSE*****************************************************
  double SW[5]; //switcher to turn various terms on an off
  SW[0] = 1.0;  //time derivative term
  SW[1] = 1.0;  //Diffusion term
  SW[2] = 1.0;  //Electric potential term
  SW[3] = 1.0;  //Convection term
  SW[4] = 1.0;  //Growth term
  //********************DEBUG PURPOSE*****************************************************

  Idx2 = 1.0/dx/dx;
  Idy2 = 1.0/dy/dy;

  tmp1 = 0.0;
  tmp2 = 0.0;
  tmp3 = 0.0;
  tmp4 = 0.0;

  if(CP.EP_flag == 1)
    EP_ind = 1.0;
  else 
    EP_ind = 0.0;

  Node_Data phis_np1(phis_np1_sin);
  Node_Data phis_n(phis_n_sin);

  Confine_Node_Data(phis_np1, PHIS_reg, 1.0);
  Confine_Node_Data(phis_n, PHIS_reg, 1.0);

  Node_Data_Ext_N_x_N_y C_n_E(C_n);
  Node_Data_Ext_N_x_N_y phs_np1_E(phis_np1);
  Node_Data_Ext_N_x_N_y EP_E(EP);
          
  Node_Data_x_Ave_Ext_N phs_n_x(phis_n);
  Node_Data_y_Ave_Ext_N phs_n_y(phis_n);
    
  Node_Data_x_Ave_Ext_N phs_np1_x(phis_np1);
  Node_Data_y_Ave_Ext_N phs_np1_y(phis_np1);

  Node_Data_x_Ave_Ext_N Dn_x(D_n);
  Node_Data_y_Ave_Ext_N Dn_y(D_n);

  Node_Data_x_Ave_Ext_N Dnp1_x(D_np1);
  Node_Data_y_Ave_Ext_N Dnp1_y(D_np1);

  for(j = 0; j <= Ny; j++)
    {
      center = j*5 + 2;
      for(i = 1; i <= Nx; i++)
	{
          //calculate the terms due to electric potential gradient
	  if(CP.EP_flag == 1)
	    {
	      tmp1 = -a*Z*Idx2*.5*(EP_E.get(i+1,j) - EP_E.get(i,j));  //(i+1,j)
	      tmp2 =  a*Z*Idx2*.5*(EP_E.get(i,j) - EP_E.get(i-1,j));  //(i-1,j)
	      tmp3 = -a*Z*Idy2*.5*(EP_E.get(i,j+1) - EP_E.get(i,j));  //(i,j+1)      
	      tmp4 =  a*Z*Idy2*.5*(EP_E.get(i,j) - EP_E.get(i,j-1));  //(i,j-1)      
	    }

	  A[center][i-1] = SW[0]*phis_np1[i][j]/dt +
	                   SW[1]*alpha*((Dnp1_x[i][j]*phs_np1_x[i][j] + Dnp1_x[i+1][j]*phs_np1_x[i+1][j])*Idx2  +
		                     (Dnp1_y[i][j]*phs_np1_y[i][j] + Dnp1_y[i][j+1]*phs_np1_y[i][j+1])*Idy2) +
	                   SW[2]*EP_ind*(tmp1 + tmp2 + tmp3 + tmp4);

	  A[center-1][i-1] =  SW[3]*-.5*un[i][j]*phs_np1_E.get(i-1,j)/dx - SW[1]*alpha*Dnp1_x[i][j]*phs_np1_x[i][j]*Idx2 + SW[2]*EP_ind*tmp2;
	  A[center+1][i-1] =  SW[3]*.5*un[i][j]*phs_np1_E.get(i+1,j)/dx - SW[1]*alpha*Dnp1_x[i+1][j]*phs_np1_x[i+1][j]*Idx2 + SW[2]*EP_ind*tmp1;
	  A[center-2][i-1] =  SW[3]*-.5*vn[i][j]*phs_np1_E.get(i,j-1)/dy - SW[1]*alpha*Dnp1_y[i][j]*phs_np1_y[i][j]*Idy2 + SW[2]*EP_ind*tmp4;
	  A[center+2][i-1] =  SW[3]*.5*vn[i][j]*phs_np1_E.get(i,j+1)/dy - SW[1]*alpha*Dnp1_y[i][j+1]*phs_np1_y[i][j+1]*Idy2 + SW[2]*EP_ind*tmp3;

	  //calculate the RHS
	  curr = j*LX + i - 1;

	  F[curr] = SW[0]*phis_n[i][j]*C_n[i][j]/dt + SW[4]*Growth(k, phi_b[i][j], Ur[i][j], S[i][j], c[i][j]); 
 
	  //add the explict part from the diffusion term from Crank-Nicolson scheme
	  tmp1 = SW[1]*(1.0 - alpha)*(Idx2*(Dn_x[i+1][j]*phs_n_x[i+1][j]*(C_n_E.get(i+1,j) - C_n_E.get(i,j)) - 
			                    Dn_x[i][j]*phs_n_x[i][j]*(C_n_E.get(i,j) - C_n_E.get(i-1,j))) +
		                      Idy2*(Dn_y[i][j+1]*phs_n_y[i][j+1]*(C_n_E.get(i,j+1) - C_n_E.get(i,j)) - 
			                    Dn_y[i][j]*phs_n_y[i][j]*(C_n_E.get(i,j) - C_n_E.get(i,j-1))));

	  F[curr] += tmp1;
	}
    }

  //Enfore the no-flux BC at j = 0 and j = Ny,
  //add the coef of j = -1 to coef of j = 1, j = Ny + 1 to j = Ny - 1

  j = 0;
  center = j*5 + 2;

  for(i = 1; i <= Nx; i++)
    {
      A[center+2][i-1] += A[center-2][i-1];
      A[center-2][i-1] = 0.0;
    }

  j = Ny;
  center = j*5 + 2;

  for(i = 1; i <= Nx; i++)
    {
      A[center-2][i-1] += A[center+2][i-1];
      A[center+2][i-1] = 0.0;
    }

  // enforce  Dirichelet BC at x = 0 and x = L by substrating the known terms

  for(j = 0; j <= Ny; j++)
    {
      center = j*5 + 2;

      i = 1;
 
      curr = j*LX + i - 1;

      F[curr] -= A[center-1][i-1]*Lc[j]; //(CP.Ini_c + CP.Ini_ce*sin(CP.Ini_c_omega*t_current)); 

      i = Nx;

      curr = j*LX + i - 1;

      F[curr] -= A[center+1][i-1]*Rc[j];
    }
}

                     
//************************************************************************************************************************

void Build_Slow_Process_Upwind(c_Mat& A, valarray<double>& F, int k, const Node_Data& D_n, const Node_Data& D_np1, double Z, const Node_Data& C_n, 
                        const valarray<double>& Lc, const valarray<double>& Rc, const Node_Data& phi_b, const Node_Data& Ur, const Node_Data& S, 
                        const Node_Data& c, const Node_Data& phis_np1_sin, const Node_Data& phis_n_sin, const Node_Data& un, const Node_Data& vn, 
                        const Node_Data& EP, double a, double dx, double dy, double dt, const Control_Parameter& CP){


  int i, j, Nx, Ny, center, curr, LX;

  Nx = A.get_Nx();
  Ny = A.get_Ny();
  LX = A.get_LX();

  if(A.get_dim() != F.size() || Lc.size() != Ny+1 || Rc.size() != Ny + 1)
    {
      cout << endl << "Dimension incompatible in Build_Slow_Process_Upwind(), ERROR at line " << __LINE__ << " of file " 
           << __FILE__ << " !! " << endl;

      exit(-1);
    }

 
  double Idx2, Idy2, tmp0, tmp1, tmp2, tmp3, tmp4, EP_ind, alpha, sg_u, sg_v;

  alpha = CP.alph;

  //********************DEBUG PURPOSE*****************************************************
  double SW[5]; //switcher to turn various terms on an off
  SW[0] = 1.0;  //time derivative term
  SW[1] = 1.0;  //Diffusion term
  SW[2] = 1.0;  //Electric potential term
  SW[3] = 1.0;  //Convection term
  SW[4] = 1.0;  //Growth term
  //********************DEBUG PURPOSE*****************************************************

  Idx2 = 1.0/dx/dx;
  Idy2 = 1.0/dy/dy;

  tmp1 = 0.0;
  tmp2 = 0.0;
  tmp3 = 0.0;
  tmp4 = 0.0;

  if(CP.EP_flag == 1)
    EP_ind = 1.0;
  else 
    EP_ind = 0.0;

  Node_Data phis_np1(phis_np1_sin);
  Node_Data phis_n(phis_n_sin);

  Confine_Node_Data(phis_np1, PHIS_reg, 1.0);
  Confine_Node_Data(phis_n, PHIS_reg, 1.0);

  Node_Data_Ext_N_x_N_y C_n_E(C_n);
  Node_Data_Ext_N_x_N_y phs_np1_E(phis_np1);
  Node_Data_Ext_N_x_N_y EP_E(EP);
          
  Node_Data_x_Ave_Ext_N phs_n_x(phis_n);
  Node_Data_y_Ave_Ext_N phs_n_y(phis_n);
    
  Node_Data_x_Ave_Ext_N phs_np1_x(phis_np1);
  Node_Data_y_Ave_Ext_N phs_np1_y(phis_np1);

  Node_Data_x_Ave_Ext_N Dn_x(D_n);
  Node_Data_y_Ave_Ext_N Dn_y(D_n);

  Node_Data_x_Ave_Ext_N Dnp1_x(D_np1);
  Node_Data_y_Ave_Ext_N Dnp1_y(D_np1);

  for(j = 0; j <= Ny; j++)
    {
      center = j*5 + 2;
      for(i = 1; i <= Nx; i++)
	{
          //calculate the terms due to electric potential gradient
	  if(CP.EP_flag == 1)
	    {
	      tmp1 = -a*Z*Idx2*.5*(EP_E.get(i+1,j) - EP_E.get(i,j));  //(i+1,j)
	      tmp2 =  a*Z*Idx2*.5*(EP_E.get(i,j) - EP_E.get(i-1,j));  //(i-1,j)
	      tmp3 = -a*Z*Idy2*.5*(EP_E.get(i,j+1) - EP_E.get(i,j));  //(i,j+1)      
	      tmp4 =  a*Z*Idy2*.5*(EP_E.get(i,j) - EP_E.get(i,j-1));  //(i,j-1)    

              tmp0 = (tmp1 + tmp2 + tmp3 + tmp4)*D_np1[i][j];

              tmp1 = tmp1*D_np1[i+1][j];
	      tmp2 = tmp2*D_np1[i-1][j];
	      tmp3 = tmp3*D_np1[i][j+1];
	      tmp4 = tmp4*D_np1[i][j-1];  
	    }

          //sign of the velocity at (i,j) node
	  sg_u = sign(un[i][j]);
	  sg_v = sign(vn[i][j]);

	  A[center][i-1] = SW[0]*phis_np1[i][j]/dt +
	                   SW[1]*alpha*((Dnp1_x[i][j]*phs_np1_x[i][j] + Dnp1_x[i+1][j]*phs_np1_x[i+1][j])*Idx2  +
		                     (Dnp1_y[i][j]*phs_np1_y[i][j] + Dnp1_y[i][j+1]*phs_np1_y[i][j+1])*Idy2) +
	                   SW[2]*EP_ind*tmp0 + SW[3]*phs_np1_E.get(i,j)*(un[i][j]*sg_u/dx + vn[i][j]*sg_v/dy);

	  A[center-1][i-1] =  SW[3]*(-.5)*un[i][j]*(1.0 + sg_u)*phs_np1_E.get(i-1,j)/dx - SW[1]*alpha*Dnp1_x[i][j]*phs_np1_x[i][j]*Idx2 + SW[2]*EP_ind*tmp2;
	  A[center+1][i-1] =  SW[3]*.5*un[i][j]*(1.0 - sg_u)*phs_np1_E.get(i+1,j)/dx - SW[1]*alpha*Dnp1_x[i+1][j]*phs_np1_x[i+1][j]*Idx2 + SW[2]*EP_ind*tmp1;
	  A[center-2][i-1] =  SW[3]*(-.5)*vn[i][j]*(1.0 + sg_v)*phs_np1_E.get(i,j-1)/dy - SW[1]*alpha*Dnp1_y[i][j]*phs_np1_y[i][j]*Idy2 + SW[2]*EP_ind*tmp4;
	  A[center+2][i-1] =  SW[3]*.5*vn[i][j]*(1.0 - sg_v)*phs_np1_E.get(i,j+1)/dy - SW[1]*alpha*Dnp1_y[i][j+1]*phs_np1_y[i][j+1]*Idy2 + SW[2]*EP_ind*tmp3;

	  //calculate the RHS
	  curr = j*LX + i - 1;

	  F[curr] = SW[0]*phis_n[i][j]*C_n[i][j]/dt + SW[4]*Growth(k, phi_b[i][j], Ur[i][j], S[i][j], c[i][j]); 
 
	  //add the explict part from the diffusion term from Crank-Nicolson scheme
	  tmp1 = SW[1]*(1.0 - alpha)*(Idx2*(Dn_x[i+1][j]*phs_n_x[i+1][j]*(C_n_E.get(i+1,j) - C_n_E.get(i,j)) - 
			                    Dn_x[i][j]*phs_n_x[i][j]*(C_n_E.get(i,j) - C_n_E.get(i-1,j))) +
		                      Idy2*(Dn_y[i][j+1]*phs_n_y[i][j+1]*(C_n_E.get(i,j+1) - C_n_E.get(i,j)) - 
			                    Dn_y[i][j]*phs_n_y[i][j]*(C_n_E.get(i,j) - C_n_E.get(i,j-1))));

	  F[curr] += tmp1;
	}
    }

  //Enfore the no-flux BC at j = 0 and j = Ny,
  //add the coef of j = -1 to coef of j = 1, j = Ny + 1 to j = Ny - 1

  j = 0;
  center = j*5 + 2;

  for(i = 1; i <= Nx; i++)
    {
      A[center+2][i-1] += A[center-2][i-1];
      A[center-2][i-1] = 0.0;
    }

  j = Ny;
  center = j*5 + 2;

  for(i = 1; i <= Nx; i++)
    {
      A[center-2][i-1] += A[center+2][i-1];
      A[center+2][i-1] = 0.0;
    }

  if(CP.velo_flag == 1) //shear flow, enforce  Dirichelet BC at x = 0 and x = L by substrating the known terms
    {
      for(j = 0; j <= Ny; j++)
	{
	  center = j*5 + 2;

	  i = 1;
 
	  curr = j*LX + i - 1;

	  F[curr] -= A[center-1][i-1]*Lc[j]; //(CP.Ini_c + CP.Ini_ce*sin(CP.Ini_c_omega*t_current)); 

	  i = Nx;

	  curr = j*LX + i - 1;

	  F[curr] -= A[center+1][i-1]*Rc[j];
	}
    }
  else if(CP.velo_flag == 0) //no flow, enforce no-flux BC at x = 0 and x = L
    {
      for(j = 0; j <= Ny; j++)
	{
	  center = j*5 + 2;

	  i = 1;

          A[center+1][i-1] += A[center-1][i-1];
 
	  i = Nx;

          A[center-1][i-1] += A[center+1][i-1];
	}
    }
}


//************************************************************************************************************************
void Build_Slow_Process_FluxBC(c_Mat& A, valarray<double>& F, int k, const Node_Data& D_n, const Node_Data& D_np1, double Z, const Node_Data& C_n,  const valarray<double>& Rc, 
                               const Node_Data& phi_b, const Node_Data& Ur, const Node_Data& S, const Node_Data& c, const Node_Data& phis_np1_sin, 
                               const Node_Data& phis_n_sin, const Node_Data& un, const Node_Data& vn, const Node_Data& EP, 
                               double a, double dx, double dy, double dt, double C_Ini, const Control_Parameter& CP){

  int i, j, Nx, Ny, center, curr, LX;

  Nx = A.get_Nx();
  Ny = A.get_Ny();
  LX = A.get_LX();

  if(A.get_dim() != F.size())
    {
      cout << endl << "Dimension incompatible in Build_Slow_Process_Upwind(), ERROR at line " << __LINE__ << " of file " 
           << __FILE__ << " !! " << endl;

      exit(-1);
    }

 
  double Idx2, Idy2, tmp0, tmp1, tmp2, tmp3, tmp4, EP_ind, alpha, sg_u, sg_v;

  alpha = CP.alph;

  //********************DEBUG PURPOSE*****************************************************
  double SW[5]; //switcher to turn various terms on an off
  SW[0] = 1.0;  //time derivative term
  SW[1] = 1.0;  //Diffusion term
  SW[2] = 1.0;  //Electric potential term
  SW[3] = 1.0;  //Convection term
  SW[4] = 1.0;  //Growth term
  //********************DEBUG PURPOSE*****************************************************

  Idx2 = 1.0/dx/dx;
  Idy2 = 1.0/dy/dy;

  tmp1 = 0.0;
  tmp2 = 0.0;
  tmp3 = 0.0;
  tmp4 = 0.0;

  if(CP.EP_flag == 1)
    EP_ind = 1.0;
  else 
    EP_ind = 0.0;

  Node_Data phis_np1(phis_np1_sin);
  Node_Data phis_n(phis_n_sin);

  Confine_Node_Data(phis_np1, PHIS_reg, 1.0);
  Confine_Node_Data(phis_n, PHIS_reg, 1.0);

  Node_Data_Ext_N_x_N_y C_n_E(C_n);
  Node_Data_Ext_N_x_N_y phs_np1_E(phis_np1);
  Node_Data_Ext_N_x_N_y EP_E(EP);
  Node_Data_Ext_N_x_N_y D_np1_E(D_np1);          

  Node_Data_x_Ave_Ext_N phs_n_x(phis_n);
  Node_Data_y_Ave_Ext_N phs_n_y(phis_n);
    
  Node_Data_x_Ave_Ext_N phs_np1_x(phis_np1);
  Node_Data_y_Ave_Ext_N phs_np1_y(phis_np1);

  Node_Data_x_Ave_Ext_N Dn_x(D_n);
  Node_Data_y_Ave_Ext_N Dn_y(D_n);

  Node_Data_x_Ave_Ext_N Dnp1_x(D_np1);
  Node_Data_y_Ave_Ext_N Dnp1_y(D_np1);


  for(j = 0; j <= Ny; j++)
    {
      center = j*5 + 2;
      for(i = 0; i <= Nx+1; i++)
	{

          //calculate the terms due to electric potential gradient
	  if(CP.EP_flag == 1)
	    {
	      tmp1 = -a*Z*Idx2*.5*(EP_E.get(i+1,j) - EP_E.get(i,j));  //(i+1,j)
	      tmp2 =  a*Z*Idx2*.5*(EP_E.get(i,j) - EP_E.get(i-1,j));  //(i-1,j)
	      tmp3 = -a*Z*Idy2*.5*(EP_E.get(i,j+1) - EP_E.get(i,j));  //(i,j+1)      
	      tmp4 =  a*Z*Idy2*.5*(EP_E.get(i,j) - EP_E.get(i,j-1));  //(i,j-1)    

              tmp0 = (tmp1 + tmp2 + tmp3 + tmp4)*D_np1_E.get(i,j);

              tmp1 = tmp1*D_np1_E.get(i+1,j); 
	      tmp2 = tmp2*D_np1_E.get(i-1,j); 
	      tmp3 = tmp3*D_np1_E.get(i,j+1); 
	      tmp4 = tmp4*D_np1_E.get(i,j-1); 
	    }


          //sign of the velocity at (i,j) node
	  sg_u = sign(un[i][j]);
	  sg_v = sign(vn[i][j]);

	  A[center][i] = SW[0]*phis_np1[i][j]/dt +
	                   SW[1]*alpha*((Dnp1_x[i][j]*phs_np1_x[i][j] + Dnp1_x[i+1][j]*phs_np1_x[i+1][j])*Idx2  +
		                     (Dnp1_y[i][j]*phs_np1_y[i][j] + Dnp1_y[i][j+1]*phs_np1_y[i][j+1])*Idy2) +
	                   SW[2]*EP_ind*tmp0 + SW[3]*phs_np1_E.get(i,j)*(un[i][j]*sg_u/dx + vn[i][j]*sg_v/dy);

	  A[center-1][i] =  SW[3]*(-.5)*un[i][j]*(1.0 + sg_u)*phs_np1_E.get(i-1,j)/dx - SW[1]*alpha*Dnp1_x[i][j]*phs_np1_x[i][j]*Idx2 + SW[2]*EP_ind*tmp2;
	  A[center+1][i] =  SW[3]*.5*un[i][j]*(1.0 - sg_u)*phs_np1_E.get(i+1,j)/dx - SW[1]*alpha*Dnp1_x[i+1][j]*phs_np1_x[i+1][j]*Idx2 + SW[2]*EP_ind*tmp1;
	  A[center-2][i] =  SW[3]*(-.5)*vn[i][j]*(1.0 + sg_v)*phs_np1_E.get(i,j-1)/dy - SW[1]*alpha*Dnp1_y[i][j]*phs_np1_y[i][j]*Idy2 + SW[2]*EP_ind*tmp4;
	  A[center+2][i] =  SW[3]*.5*vn[i][j]*(1.0 - sg_v)*phs_np1_E.get(i,j+1)/dy - SW[1]*alpha*Dnp1_y[i][j+1]*phs_np1_y[i][j+1]*Idy2 + SW[2]*EP_ind*tmp3;

	  //calculate the RHS
	  curr = j*LX + i;

	  F[curr] = SW[0]*phis_n[i][j]*C_n[i][j]/dt + SW[4]*Growth(k, phi_b[i][j], Ur[i][j], S[i][j], c[i][j]); 
 
	  //add the explict part from the diffusion term from Crank-Nicolson scheme
	  tmp1 = SW[1]*(1.0 - alpha)*(Idx2*(Dn_x[i+1][j]*phs_n_x[i+1][j]*(C_n_E.get(i+1,j) - C_n_E.get(i,j)) - 
			                    Dn_x[i][j]*phs_n_x[i][j]*(C_n_E.get(i,j) - C_n_E.get(i-1,j))) +
		                      Idy2*(Dn_y[i][j+1]*phs_n_y[i][j+1]*(C_n_E.get(i,j+1) - C_n_E.get(i,j)) - 
			                    Dn_y[i][j]*phs_n_y[i][j]*(C_n_E.get(i,j) - C_n_E.get(i,j-1))));

	  F[curr] += tmp1;
	}
    }



  //Enfore the no-flux BC at j = 0 and j = Ny,
  //add the coef of j = -1 to coef of j = 1, j = Ny + 1 to j = Ny - 1

  j = 0;
  center = j*5 + 2;

  for(i = 0; i <= Nx+1; i++)
    {
      A[center+2][i] += A[center-2][i];
      A[center-2][i] = 0.0;
    }

  j = Ny;
  center = j*5 + 2;

  for(i = 0; i <= Nx+1; i++)
    {
      A[center-2][i] += A[center+2][i];
      A[center+2][i] = 0.0;
    }

  //Boundary condition in x-direction is enforced to balance flux
  for(j = 0; j <= Ny; j++)
    {
     center = j*5 + 2;

     i = 0;   //Left boundary, balance flux -D*c_x + vx*c = vx_0*c_0

     if(CP.velo_flag == 1) //shear flow, out-flow BC  
       {
	  curr = j*LX + i;

	  F[curr] -= A[center-1][i]*C_Ini;
       }
     else if(CP.velo_flag == 0) //No external flow, no-flux BC
       {
         A[center+1][i] += A[center-1][i];
       }

     // central difference for gradient
//      A[center+1][i] += A[center-1][i];
//      A[center][i] -= A[center-1][i]*2*dx*un[0][j]/D_np1[0][j];
//      F[curr] -= A[center-1][i]*2*dx*un[0][j]*C_Ini/D_np1[0][j];

     // one-side difference for gradient
//      A[center][i] += A[center-1][i]*(1.0 - dx*un[0][j]/D_np1[0][j]);
//      F[curr] -= A[center-1][i]*dx*un[0][j]*C_Ini/D_np1[0][j];

   

     i = Nx+1; //Right boundary
     if(CP.velo_flag == 1) //shear flow, out-flow BC  
       {
	  curr = j*LX + i;

	  F[curr] -= A[center+1][i]*Rc[j];
       }
     else if(CP.velo_flag == 0) //No external flow, no-flux BC
       {
         A[center-1][i] += A[center+1][i];
       }
//      A[center][i] += A[center+1][i]*(Dnp1_x[Nx+2][j] + Dnp1_x[Nx+1][j] + dx*un[Nx+1][j])/Dnp1_x[Nx+2][j];
//      A[center-1][i] -= A[center+1][i]*(Dnp1_x[Nx+1][j] + dx*un[Nx][j])/Dnp1_x[Nx+2][j];

    }
}

//************************************************************************************************************************
void Build_Diffusion(c_NoFlux_x_NoFlux_y_mat& A, valarray<double>& F, const Node_Data& U_n, 
                     double dx, double dy, double dt, double alpha, double eps_dif){

  if(A.get_dim() != F.size() )
    {
      cout << endl << "Dimension incompatible in Build_Diffusion(), ERROR in Build_2D_new.cpp !! " 
           << endl;
      exit(-1);
    }

  int Nx, Ny, i, j, center, curr;

  double mu_x, mu_y;


  Nx = A.get_Nx();
  Ny = A.get_Ny();

  mu_x = dt/dx/dx;
  mu_y = dt/dy/dy;

  Node_Data_Ext_N_x_N_y U(U_n);

     
  for(j = 0; j <= Ny; j++)
    {
      center = j*5 + 2;

      for(i = 0; i <= Nx + 1; i++)
	{
          //calculate the entries of the matrix
	  A[center-2][i] = - eps_dif*alpha*mu_y;
          A[center-1][i] = - eps_dif*alpha*mu_x;
          A[center][i]   = 1.0 + 2.0*eps_dif*alpha*(mu_x + mu_y);
          A[center+1][i] = - eps_dif*alpha*mu_x;
	  A[center+2][i] = - eps_dif*alpha*mu_y;

	  //calculate the entry of the RHS vector
          curr = j*(Nx+2) + i;
          F[curr] =  eps_dif*(1 - alpha)*mu_y*U.get(i,j-1) + eps_dif*(1 - alpha)*mu_x*U.get(i-1,j) +
                    (1 - 2.0*eps_dif*(1 - alpha)*(mu_x + mu_y))*U.get(i,j) + 
	            eps_dif*(1 - alpha)*mu_x*U.get(i+1,j) + eps_dif*(1 - alpha)*mu_y*U.get(i,j+1);
	}
    }

  //now reinforce the no-flux boundary condition
      
  j = 0;
  center = j*5 + 2;

  for(i = 0; i <= Nx + 1; i++)
    {
      A[center+2][i] += A[center-2][i];
      A[center-2][i] = 0.0;
    }

  j = Ny;
  center = j*5 + 2;

  for(i = 0; i <= Nx + 1; i++)
    {
      A[center-2][i] += A[center+2][i];
      A[center+2][i] = 0.0;
    }

  i = 0;
  for(j = 0; j <= Ny; j++)
    {
      center = j*5 + 2;
      A[center+1][i] += A[center-1][i];
      A[center-1][i] = 0.0;
    }

  i = Nx + 1;
  for(j = 0; j <= Ny; j++)
    {
      center = j*5 + 2;
      A[center-1][i] += A[center+1][i];
      A[center+1][i] = 0.0;
    }

}


//************************************************************************************************************************

double Diffusion_coef_phic(double D, double phic, double phic_tol, double delta){

  if(phic_tol <= 0 || phic_tol >= 1)
    {
      cout << endl << "phic_tol is not between 0 and 1, error at line " << __LINE__ << " of file " << __FILE__ << endl;
      exit(-1);
    }

  double result;

  result = D + ((1.0 - delta)*D/2)*(-tanh(5.0*(phic - phic_tol)/(1.0 - phic_tol)) - 1.0);

  return result;
}


//************************************************************************************************************************
 
void Diffusion_Node(Node_Data& DN, const Node_Data& phic, double D, double phic_tol, double delta){

  int i, j, Nx, Ny;

  Nx = phic.get_Nx();
  Ny = phic.get_Ny();

  for(i = 0; i <= Nx+1; i++)
    for(j = 0; j <= Ny; j++)
      DN[i][j] = Diffusion_coef_phic(D, phic[i][j], phic_tol, delta);
}

//************************************************************************************************************************

double Eff_Diff_coef_phic_sis(double D, double phic, double delta){

  return D*((1.0 - phic)/(1 + 0.5*phic) + delta)/(1.0 + delta);

}

//************************************************************************************************************************
void Eff_Diff_Node_sis(Node_Data& DN, const Node_Data& phic, double D, double delta){

  int i, j, Nx, Ny;

  Nx = phic.get_Nx();
  Ny = phic.get_Ny();

  for(i = 0; i <= Nx+1; i++)
    for(j = 0; j <= Ny; j++)
      DN[i][j] = Eff_Diff_coef_phic_sis(D, phic[i][j], delta);
}

//************************************************************************************************************************
void Build_EP_N_x_N_y(c_Mat& A, valarray<double>& F, const Node_Data_Ext_IO_x_N_y& NH4, const Node_Data_Ext_IO_x_N_y& HCO3, const Node_Data_Ext_IO_x_N_y& CO3, 
                      const Node_Data_Ext_IO_x_N_y& Ca, const Node_Data_Ext_IO_x_N_y& Cl, const Node_Data_Ext_IO_x_N_y& Dv_nh4, const Node_Data_Ext_IO_x_N_y& Dv_hco3, 
                      const Node_Data_Ext_IO_x_N_y& Dv_co3, const Node_Data_Ext_IO_x_N_y& Dv_ca, const Node_Data_Ext_IO_x_N_y& Dv_cl, const Node_Data_Ext_IO_x_N_y& phis,
                      double dx, double dy, double a){


  int i, j, Nx, Ny, center, curr, LX;

  Nx = A.get_Nx();
  Ny = A.get_Ny();
  LX = A.get_LX();

  if(A.get_dim() != F.size())
    {
      cout << endl << "Dimension incompatible in Build_EP_N_x_N_y(), ERROR at line " << __LINE__ << " of file " 
           << __FILE__ << " !! " << endl;

      exit(-1);
    }

  double Idx2, Idy2, tmp;

  Idx2 = 1.0/dx/dx;
  Idy2 = 1.0/dy/dy;

  //calculate the average in x and y direction
  
  Node_Data_x_Ave_Ext_IO NH4_x(NH4);
  Node_Data_y_Ave_Ext_N  NH4_y(NH4);
  Node_Data_x_Ave_Ext_IO D_NH4_x(Dv_nh4);
  Node_Data_y_Ave_Ext_N  D_NH4_y(Dv_nh4);
    
  Node_Data_x_Ave_Ext_IO HCO3_x(HCO3);
  Node_Data_y_Ave_Ext_N  HCO3_y(HCO3);
  Node_Data_x_Ave_Ext_IO D_HCO3_x(Dv_hco3);
  Node_Data_y_Ave_Ext_N  D_HCO3_y(Dv_hco3);
    
  Node_Data_x_Ave_Ext_IO CO3_x(CO3);
  Node_Data_y_Ave_Ext_N  CO3_y(CO3);
  Node_Data_x_Ave_Ext_IO D_CO3_x(Dv_co3);
  Node_Data_y_Ave_Ext_N  D_CO3_y(Dv_co3);
    
  Node_Data_x_Ave_Ext_IO Ca_x(Ca);
  Node_Data_y_Ave_Ext_N  Ca_y(Ca);
  Node_Data_x_Ave_Ext_IO D_Ca_x(Dv_ca);
  Node_Data_y_Ave_Ext_N  D_Ca_y(Dv_ca);
    
  Node_Data_x_Ave_Ext_IO Cl_x(Cl);
  Node_Data_y_Ave_Ext_N  Cl_y(Cl);
  Node_Data_x_Ave_Ext_IO D_Cl_x(Dv_cl);
  Node_Data_y_Ave_Ext_N  D_Cl_y(Dv_cl);

  Node_Data_x_Ave_Ext_IO phis_x(phis);
  Node_Data_y_Ave_Ext_N  phis_y(phis);


  for(j = 0; j <= Ny; j++)
    {
      center = j*5 + 2;
      for(i = 0; i <= Nx+1 ; i++)
	{
	 
	  A[center-1][i] = Idx2*(Z_nh4*Z_nh4*NH4_x[i][j] + Z_hco3*Z_hco3*HCO3_x[i][j] + Z_co3*Z_co3*CO3_x[i][j] +
                                 Z_ca*Z_ca*Ca_x[i][j] + Z_cl*Z_cl*Cl_x[i][j]);

                              
	  A[center+1][i] = Idx2*(Z_nh4*Z_nh4*NH4_x[i+1][j] + Z_hco3*Z_hco3*HCO3_x[i+1][j] + Z_co3*Z_co3*CO3_x[i+1][j] +
                                 Z_ca*Z_ca*Ca_x[i+1][j] + Z_cl*Z_cl*Cl_x[i+1][j]);

	  A[center-2][i] = Idy2*(Z_nh4*Z_nh4*NH4_y[i][j] + Z_hco3*Z_hco3*HCO3_y[i][j] + Z_co3*Z_co3*CO3_y[i][j] +
                                 Z_ca*Z_ca*Ca_y[i][j] + Z_cl*Z_cl*Cl_y[i][j]);

	  A[center+2][i] = Idy2*(Z_nh4*Z_nh4*NH4_y[i][j+1] + Z_hco3*Z_hco3*HCO3_y[i][j+1] + Z_co3*Z_co3*CO3_y[i][j+1] +
                                 Z_ca*Z_ca*Ca_y[i][j+1] + Z_cl*Z_cl*Cl_y[i][j+1]);

          A[center][i] = -(A[center-1][i] + A[center+1][i] + A[center-2][i] + A[center+2][i]);

	  //calculate the RHS
	  curr = j*LX + i;

	  F[curr] = (-1.0/a)*(Idx2*(phis_x[i+1][j]*(Z_nh4*D_NH4_x[i+1][j]*(NH4.get(i+1,j) - NH4.get(i,j)) + 
                                                    Z_hco3*D_HCO3_x[i+1][j]*(HCO3.get(i+1,j) - HCO3.get(i,j)) + 
                                                    Z_co3*D_CO3_x[i+1][j]*(CO3.get(i+1,j) - CO3.get(i,j)) +  
                                                    Z_ca*D_Ca_x[i+1][j]*(Ca.get(i+1,j) - Ca.get(i,j)) +
                                                    Z_cl*D_Cl_x[i+1][j]*(Cl.get(i+1,j) - Cl.get(i,j)))  -
                                      phis_x[i][j]*(Z_nh4*D_NH4_x[i][j]*(NH4.get(i,j) - NH4.get(i-1,j)) + 
                                                    Z_hco3*D_HCO3_x[i][j]*(HCO3.get(i,j) - HCO3.get(i-1,j)) + 
                                                    Z_co3*D_CO3_x[i][j]*(CO3.get(i,j) - CO3.get(i-1,j)) +  
                                                    Z_ca*D_Ca_x[i][j]*(Ca.get(i,j) - Ca.get(i-1,j)) +
                                                    Z_cl*D_Cl_x[i][j]*(Cl.get(i,j) - Cl.get(i-1,j)))) +
                              Idy2*(phis_y[i][j+1]*(Z_nh4*D_NH4_y[i][j+1]*(NH4.get(i,j+1) - NH4.get(i,j)) + 
                                                    Z_hco3*D_HCO3_y[i][j+1]*(HCO3.get(i,j+1) - HCO3.get(i,j)) + 
                                                    Z_co3*D_CO3_y[i][j+1]*(CO3.get(i,j+1) - CO3.get(i,j)) +  
                                                    Z_ca*D_Ca_y[i][j+1]*(Ca.get(i,j+1) - Ca.get(i,j)) +
                                                    Z_cl*D_Cl_y[i][j+1]*(Cl.get(i,j+1) - Cl.get(i,j)))  -
                                      phis_y[i][j]*(Z_nh4*D_NH4_y[i][j]*(NH4.get(i,j) - NH4.get(i,j-1)) + 
                                                    Z_hco3*D_HCO3_y[i][j]*(HCO3.get(i,j) - HCO3.get(i,j-1)) + 
                                                    Z_co3*D_CO3_y[i][j]*(CO3.get(i,j) - CO3.get(i,j-1)) +  
                                                    Z_ca*D_Ca_y[i][j]*(Ca.get(i,j) - Ca.get(i,j-1)) +
                                                    Z_cl*D_Cl_y[i][j]*(Cl.get(i,j) - Cl.get(i,j-1)))));
                                    
	}
    }

  //enforce the no-flux boundary condition 
  for(i = 0; i <= Nx+1; i++)
    {
      j = 0;
      center = j*5 + 2;

      A[center+2][i] += A[center-2][i];

      j = Ny;
      center = j*5 + 2;

      A[center-2][i] += A[center+2][i];
    }

  for(j = 0; j <= Ny; j++)
    {
      center = j*5 + 2;

      i = 0;
      A[center+1][i] += A[center-1][i];

      i = Nx+1;
      A[center-1][i] += A[center+1][i];
    }

}

//************************************************************************************************************************
Node_Data Reflect_Node_Data_LR(const Node_Data& U){

  Node_Data V(U);

  int j, Nx, Ny;

  Nx = V.get_Nx();
  Ny = V.get_Ny();

  for(j = 0; j <= Ny; j++)
    {
      V[0][j] = V[2][j];
      V[Nx+1][j] = V[Nx-1][j];
    }

  return V;
}

//************************************************************************************************************************

void Build_EP_Interior(c_Mat& A, valarray<double>& F, const Node_Data& NH4, const Node_Data& HCO3, const Node_Data& CO3, const Node_Data& Ca, const Node_Data& Cl, const Node_Data& H, const Node_Data& OH,
                       const Node_Data& Dv_nh4, const Node_Data& Dv_hco3, const Node_Data& Dv_co3, const Node_Data& Dv_ca, const Node_Data& Dv_cl, const Node_Data& Dv_h, const Node_Data& Dv_oh,
                       const Node_Data& phis, double dx, double dy, double dt, double a, double& source, int Inter_file_NO){


  int i, j, Nx, Ny, center, curr, LX;

  Nx = A.get_Nx();
  Ny = A.get_Ny();
  LX = A.get_LX();

  if(A.get_dim() != F.size())
    {
      cout << endl << "Dimension incompatible in Build_EP_N_x_N_y(), ERROR at line " << __LINE__ << " of file " 
           << __FILE__ << " !! " << endl;

      exit(-1);
    }

  source = 0.0;

  //***********************************  DEBUG  *******************************
  Node_Data Extra_source(H);
  char debfile[80];
  //***********************************  DEBUG  *******************************

  double Idx2, Idy2, tmp;

  //********************DEBUG PURPOSE*****************************************************
  double SW[7]; //switcher to turn various terms on an off
  SW[0] = 1.0;  //NH4
  SW[1] = 1.0;  //HCO3
  SW[2] = 1.0;  //CO3
  SW[3] = 1.0;  //Ca
  SW[4] = 1.0;  //Cl
  SW[5] = 1.0;  //H+
  SW[6] = 1.0;  //OH-
  //********************DEBUG PURPOSE*****************************************************

  Idx2 = 1.0/dx/dx;
  Idy2 = 1.0/dy/dy;
 
//   //calculate the average in x and y direction
//   //adjust the input Node_Data variables so they satisfy no-flux BC at i = 1 and i = Nx
//   Node_Data NH4(NH4_in); //(Reflect_Node_Data_LR(NH4_in));
//   Node_Data Dv_nh4(Dv_nh4_in); //(Reflect_Node_Data_LR(Dv_nh4_in));

//   Node_Data HCO3(HCO3_in); //(Reflect_Node_Data_LR(HCO3_in));
//   Node_Data Dv_hco3(Dv_hco3_in); //(Reflect_Node_Data_LR(Dv_hco3_in));

//   Node_Data CO3(CO3_in); //(Reflect_Node_Data_LR(CO3_in));
//   Node_Data Dv_co3(Dv_co3_in); //(Reflect_Node_Data_LR(Dv_co3_in));

//   Node_Data Ca(Ca_in); //(Reflect_Node_Data_LR(Ca_in));
//   Node_Data Dv_ca(Dv_ca_in); //(Reflect_Node_Data_LR(Dv_ca_in));

//   Node_Data Cl(Cl_in); //(Reflect_Node_Data_LR(Cl_in));
//   Node_Data Dv_cl(Dv_cl_in); //(Reflect_Node_Data_LR(Dv_cl_in));

//   Node_Data phis(phis_in); //(Reflect_Node_Data_LR(phis_in));
  
  Node_Data_Ext_N_x_N_y NH4_E(NH4);
  Node_Data_x_Ave_Ext_N NH4_x(NH4);
  Node_Data_y_Ave_Ext_N NH4_y(NH4);
  Node_Data_x_Ave_Ext_N D_NH4_x(Dv_nh4);
  Node_Data_y_Ave_Ext_N D_NH4_y(Dv_nh4);
    
  Node_Data_Ext_N_x_N_y HCO3_E(HCO3);
  Node_Data_x_Ave_Ext_N HCO3_x(HCO3);
  Node_Data_y_Ave_Ext_N HCO3_y(HCO3);
  Node_Data_x_Ave_Ext_N D_HCO3_x(Dv_hco3);
  Node_Data_y_Ave_Ext_N D_HCO3_y(Dv_hco3);
    
  Node_Data_Ext_N_x_N_y CO3_E(CO3);
  Node_Data_x_Ave_Ext_N CO3_x(CO3);
  Node_Data_y_Ave_Ext_N CO3_y(CO3);
  Node_Data_x_Ave_Ext_N D_CO3_x(Dv_co3);
  Node_Data_y_Ave_Ext_N D_CO3_y(Dv_co3);
    
  Node_Data_Ext_N_x_N_y Ca_E(Ca);
  Node_Data_x_Ave_Ext_N Ca_x(Ca);
  Node_Data_y_Ave_Ext_N Ca_y(Ca);
  Node_Data_x_Ave_Ext_N D_Ca_x(Dv_ca);
  Node_Data_y_Ave_Ext_N D_Ca_y(Dv_ca);
    
  Node_Data_Ext_N_x_N_y Cl_E(Cl);
  Node_Data_x_Ave_Ext_N Cl_x(Cl);
  Node_Data_y_Ave_Ext_N Cl_y(Cl);
  Node_Data_x_Ave_Ext_N D_Cl_x(Dv_cl);
  Node_Data_y_Ave_Ext_N D_Cl_y(Dv_cl);
  
  Node_Data_Ext_N_x_N_y H_E(H);
  Node_Data_x_Ave_Ext_N H_x(H);
  Node_Data_y_Ave_Ext_N H_y(H);
  Node_Data_x_Ave_Ext_N D_H_x(Dv_h);
  Node_Data_y_Ave_Ext_N D_H_y(Dv_h);
  
  Node_Data_Ext_N_x_N_y OH_E(OH);
  Node_Data_x_Ave_Ext_N OH_x(OH);
  Node_Data_y_Ave_Ext_N OH_y(OH);
  Node_Data_x_Ave_Ext_N D_OH_x(Dv_oh);
  Node_Data_y_Ave_Ext_N D_OH_y(Dv_oh);


  Node_Data_x_Ave_Ext_N phis_x(phis);
  Node_Data_y_Ave_Ext_N phis_y(phis);


  for(j = 0; j <= Ny; j++)
    {
      center = j*5 + 2;
      for(i = 1; i <= Nx; i++)
	{
	  A[center-1][i-1] = Idx2*(Z_nh4*Z_nh4*D_NH4_x[i][j]*NH4_x[i][j] + Z_hco3*Z_hco3*D_HCO3_x[i][j]*HCO3_x[i][j] + 
                                   Z_co3*Z_co3*D_CO3_x[i][j]*CO3_x[i][j] + Z_ca*Z_ca*D_Ca_x[i][j]*Ca_x[i][j] + Z_cl*Z_cl*D_Cl_x[i][j]*Cl_x[i][j] + 
                                   Z_h*Z_h*D_H_x[i][j]*H_x[i][j] + Z_oh*Z_oh*D_OH_x[i][j]*OH_x[i][j]);
                              
	  A[center+1][i-1] = Idx2*(Z_nh4*Z_nh4*D_NH4_x[i+1][j]*NH4_x[i+1][j] + Z_hco3*Z_hco3*D_HCO3_x[i+1][j]*HCO3_x[i+1][j] + 
                                   Z_co3*Z_co3*D_CO3_x[i+1][j]*CO3_x[i+1][j] + Z_ca*Z_ca*D_Ca_x[i+1][j]*Ca_x[i+1][j] + Z_cl*Z_cl*D_Cl_x[i+1][j]*Cl_x[i+1][j] + 
                                   Z_h*Z_h*D_H_x[i+1][j]*H_x[i+1][j] + Z_oh*Z_oh*D_OH_x[i+1][j]*OH_x[i+1][j]);

	  A[center-2][i-1] = Idy2*(Z_nh4*Z_nh4*D_NH4_y[i][j]*NH4_y[i][j] + Z_hco3*Z_hco3*D_HCO3_y[i][j]*HCO3_y[i][j] + 
                                   Z_co3*Z_co3*D_CO3_y[i][j]*CO3_y[i][j] + Z_ca*Z_ca*D_Ca_y[i][j]*Ca_y[i][j] + Z_cl*Z_cl*D_Cl_y[i][j]*Cl_y[i][j] + 
                                   Z_h*Z_h*D_H_y[i][j]*H_y[i][j] + Z_oh*Z_oh*D_OH_y[i][j]*OH_y[i][j]);

	  A[center+2][i-1] = Idy2*(Z_nh4*Z_nh4*D_NH4_y[i][j+1]*NH4_y[i][j+1] + Z_hco3*Z_hco3*D_HCO3_y[i][j+1]*HCO3_y[i][j+1] + 
                                   Z_co3*Z_co3*D_CO3_y[i][j+1]*CO3_y[i][j+1] + Z_ca*Z_ca*D_Ca_y[i][j+1]*Ca_y[i][j+1] + Z_cl*Z_cl*D_Cl_y[i][j+1]*Cl_y[i][j+1] + 
                                   Z_h*Z_h*D_H_y[i][j+1]*H_y[i][j+1] + Z_oh*Z_oh*D_OH_y[i][j+1]*OH_y[i][j+1]);

          A[center][i-1] = -(A[center-1][i-1] + A[center+1][i-1] + A[center-2][i-1] + A[center+2][i-1]);

	  //calculate the RHS
	  curr = j*LX + i-1;

	  F[curr] = (-1.0/a)*( (Idx2*(phis_x[i+1][j]*(SW[0]*Z_nh4*D_NH4_x[i+1][j]*(NH4_E.get(i+1,j) - NH4_E.get(i,j)) + 
                                                    SW[1]*Z_hco3*D_HCO3_x[i+1][j]*(HCO3_E.get(i+1,j) - HCO3_E.get(i,j)) + 
                                                    SW[2]*Z_co3*D_CO3_x[i+1][j]*(CO3_E.get(i+1,j) - CO3_E.get(i,j)) +  
                                                    SW[3]*Z_ca*D_Ca_x[i+1][j]*(Ca_E.get(i+1,j) - Ca_E.get(i,j)) +
                                                    SW[4]*Z_cl*D_Cl_x[i+1][j]*(Cl_E.get(i+1,j) - Cl_E.get(i,j)) +
                                                    SW[5]*Z_h*D_H_x[i+1][j]*(H_E.get(i+1,j) - H_E.get(i,j)) +
                                                    SW[6]*Z_oh*D_OH_x[i+1][j]*(OH_E.get(i+1,j) - OH_E.get(i,j))) -
                                      phis_x[i][j]*(SW[0]*Z_nh4*D_NH4_x[i][j]*(NH4_E.get(i,j) - NH4_E.get(i-1,j)) + 
                                                    SW[1]*Z_hco3*D_HCO3_x[i][j]*(HCO3_E.get(i,j) - HCO3_E.get(i-1,j)) + 
                                                    SW[2]*Z_co3*D_CO3_x[i][j]*(CO3_E.get(i,j) - CO3_E.get(i-1,j)) +  
                                                    SW[3]*Z_ca*D_Ca_x[i][j]*(Ca_E.get(i,j) - Ca_E.get(i-1,j)) +
                                                    SW[4]*Z_cl*D_Cl_x[i][j]*(Cl_E.get(i,j) - Cl_E.get(i-1,j)) +
                                                    SW[5]*Z_h*D_H_x[i][j]*(H_E.get(i,j) - H_E.get(i-1,j)) +
                                                    SW[6]*Z_oh*D_OH_x[i][j]*(OH_E.get(i,j) - OH_E.get(i-1,j)))) +
                                Idy2*(phis_y[i][j+1]*(SW[0]*Z_nh4*D_NH4_y[i][j+1]*(NH4_E.get(i,j+1) - NH4_E.get(i,j)) + 
                                                    SW[1]*Z_hco3*D_HCO3_y[i][j+1]*(HCO3_E.get(i,j+1) - HCO3_E.get(i,j)) + 
                                                    SW[2]*Z_co3*D_CO3_y[i][j+1]*(CO3_E.get(i,j+1) - CO3_E.get(i,j)) +  
                                                    SW[3]*Z_ca*D_Ca_y[i][j+1]*(Ca_E.get(i,j+1) - Ca_E.get(i,j)) +
                                                    SW[4]*Z_cl*D_Cl_y[i][j+1]*(Cl_E.get(i,j+1) - Cl_E.get(i,j)) +
                                                    SW[5]*Z_h*D_H_y[i][j+1]*(H_E.get(i,j+1) - H_E.get(i,j)) +
                                                    SW[6]*Z_oh*D_OH_y[i][j+1]*(OH_E.get(i,j+1) - OH_E.get(i,j)))  -
                                       phis_y[i][j]*(SW[0]*Z_nh4*D_NH4_y[i][j]*(NH4_E.get(i,j) - NH4_E.get(i,j-1)) + 
                                                    SW[1]*Z_hco3*D_HCO3_y[i][j]*(HCO3_E.get(i,j) - HCO3_E.get(i,j-1)) + 
                                                    SW[2]*Z_co3*D_CO3_y[i][j]*(CO3_E.get(i,j) - CO3_E.get(i,j-1)) +  
                                                    SW[3]*Z_ca*D_Ca_y[i][j]*(Ca_E.get(i,j) - Ca_E.get(i,j-1)) +
                                                    SW[4]*Z_cl*D_Cl_y[i][j]*(Cl_E.get(i,j) - Cl_E.get(i,j-1)) +
                                                    SW[5]*Z_h*D_H_y[i][j]*(H_E.get(i,j) - H_E.get(i,j-1)) +
						     SW[6]*Z_oh*D_OH_y[i][j]*(OH_E.get(i,j) - OH_E.get(i,j-1)))))//);			      
                                                 +  (phis[i][j]/dt)*(Z_nh4*NH4[i][j] + Z_hco3*HCO3[i][j] + Z_co3*CO3[i][j] +
                                                                     Z_ca*Ca[i][j] + Z_cl*Cl[i][j] +  Z_h*H[i][j] + Z_oh*OH[i][j]));

#ifdef DEBUG_Build
	  Extra_source[i][j] = (phis[i][j]/dt)*(Z_nh4*NH4[i][j] + Z_hco3*HCO3[i][j] + Z_co3*CO3[i][j] +
					    Z_ca*Ca[i][j] + Z_cl*Cl[i][j] +  Z_h*H[i][j] + Z_oh*OH[i][j]);
#endif
	}
    }

  //enforce the no-flux boundary condition in y-direction
  for(i = 1; i <= Nx; i++)
    {
      j = 0;
      center = j*5 + 2;

      A[center+2][i-1] += A[center-2][i-1];
      A[center-2][i-1] = 0.0;


      j = Ny;
      center = j*5 + 2;

      A[center-2][i-1] += A[center+2][i-1];
      A[center+2][i-1] = 0.0;
    }

  /*  //Homogeneous Dirichlet BC in x direction (since same is true for Chemical species */

  //no-flux BC in x-direction
  for(j = 0; j <= Ny; j++)
    {
      center = j*5 + 2;

      i = 1;
      A[center+1][i-1] += A[center-1][i-1];
      A[center-1][i-1] = 0.0;

      i = Nx;
      A[center-1][i-1] += A[center+1][i-1];
      A[center+1][i-1] = 0.0;
    }

  //calculate the source term
  for(j = 0; j <= Ny; j++)
    {
      i = 1; // outer normal derivative is dc/dn = (c[i][j] - c[i+1][j])/dx

      source += (phis[i][j]/dx)*(Z_nh4*Dv_nh4[i][j]*(NH4[i][j] - NH4[i+1][j]) + Z_hco3*Dv_hco3[i][j]*(HCO3[i][j] - HCO3[i+1][j]) +
                                     Z_co3*Dv_co3[i][j]*(CO3[i][j] - CO3[i+1][j]) + Z_ca*Dv_ca[i][j]*(Ca[i][j] - Ca[i+1][j]) +
                                     Z_cl*Dv_cl[i][j]*(Cl[i][j] - Cl[i+1][j]) + Z_h*Dv_h[i][j]*(H[i][j] - H[i+1][j]) + 
                                     Z_oh*Dv_oh[i][j]*(OH[i][j] - OH[i+1][j]));

      i = Nx; // outer normal derivative is dc/dn = (c[i][j] - c[i-1][j])/dx

      source += (phis[i][j]/dx)*(Z_nh4*Dv_nh4[i][j]*(NH4[i][j] - NH4[i-1][j]) + Z_hco3*Dv_hco3[i][j]*(HCO3[i][j] - HCO3[i-1][j]) +
                                     Z_co3*Dv_co3[i][j]*(CO3[i][j] - CO3[i-1][j]) + Z_ca*Dv_ca[i][j]*(Ca[i][j] - Ca[i-1][j]) +
                                     Z_cl*Dv_cl[i][j]*(Cl[i][j] - Cl[i-1][j]) + Z_h*Dv_h[i][j]*(H[i][j] - H[i-1][j]) + 
                                     Z_oh*Dv_oh[i][j]*(OH[i][j] - OH[i-1][j]));
    }

  source *= dy;

#ifdef DEBUG_Build
  sprintf(debfile,"%sEP_Extr_%d",SaveDir, Inter_file_NO);
  Write_Node_Data(debfile, Extra_source); //useful data only at 1 <= i <= Nx, 0 <= j <= Ny
#endif
}

//************************************************************************************************************************

double NLF_H(double H, double NT, double K3, double Kw, double CT, double K1, double K2, double Ca, double Cl, double cha_sum){

  double result, tmp0, tmp1, tmp2;

  tmp0 = 1.0 + K3*H; 
  tmp1 = H/K1 + 1.0 + K2/H;
  tmp2 = H*H/K1/K2 + H/K2 + 1.0;

  result = NT*K3*H/tmp0 + H - Kw/H - CT/tmp1 - 2.0*CT/tmp2 + 2.0*Ca - Cl - cha_sum;

  return result;
}


double Deri_NLF_H(double H, double NT, double K3, double Kw, double CT, double K1, double K2, double Ca, double Cl, double cha_sum){

  double result, tmp0, tmp1, tmp2;
  
  tmp0 = 1.0 + K3*H; 
  tmp1 = H/K1 + 1.0 + K2/H;
  tmp2 = H*H/K1/K2 + H/K2 + 1.0;

  result = NT*K3/tmp0/tmp0 + 1.0 + Kw/H/H + CT*(1.0/K1 - K2/H/H)/tmp1/tmp1 + 2.0*CT*(2.0*H/K1/K2 + 1.0/K2)/tmp2/tmp2;

  return result;
}

int Newton_Root(double& H, double NT, double K3, double Kw, double CT, double K1, double K2, double Ca, double Cl, double cha_sum,
                int& N_iter){

  double x_n, x_np1, deps, tmp1, tmp2, ini_N;

  int flag, iter, n;

  x_n = H;

  ini_N = 0.0;

  deps = numeric_limits<double>::epsilon();

  if(fabs(x_n) < deps)
    x_n = deps;


  while(1)
    {
      tmp1 = NLF_H(x_n, NT, K3, Kw, CT, K1, K2, Ca, Cl, cha_sum);
      tmp2 = Deri_NLF_H(x_n, NT, K3, Kw, CT, K1, K2, Ca, Cl, cha_sum); 
 
      if(fabs(tmp1) < deps) //already find a solution
	{
	  H = x_n;
          N_iter = 0;
          //cout << endl << "Exit before iteration! " << endl;
	  return 0;
	}
      else if(fabs(tmp2) > deps) //derivative not zero
	break;
      else 
	x_n += deps;
    }

  x_np1 = x_n;

  do
    {
      if(x_np1 <= 0)
	x_n = H*pow(10,ini_N);
 
      tmp1 = NLF_H(x_n, NT, K3, Kw, CT, K1, K2, Ca, Cl, cha_sum);
      tmp2 = Deri_NLF_H(x_n, NT, K3, Kw, CT, K1, K2, Ca, Cl, cha_sum); 

      for(n = 0; n < 20; n++)
	{
	  x_np1 = x_n - tmp1/tmp2;

	  if(fabs((x_np1 - x_n)/x_n) < 1e-8) //relative error small, find a solution!
	    {
	      break;
	    }

	  x_n = x_np1;

	  tmp1 = NLF_H(x_n, NT, K3, Kw, CT, K1, K2, Ca, Cl, cha_sum);
	  tmp2 = Deri_NLF_H(x_n, NT, K3, Kw, CT, K1, K2, Ca, Cl, cha_sum); 
	}

      ini_N -= 1.0;

      if(ini_N < -10)
	break; //didn't converge even with initial guess reduced by 1e-10
       
    }while(x_np1 <= 0);

  N_iter = ini_N;
  H = x_np1;
//   tmp1 = NLF_H(x_np1, NT, K3, Kw, CT, K1, K2, Ca, Cl);
//   cout << endl << "the residual is " << tmp1 << ", x_np1 = " << x_np1 << endl ;
 
  if(n < 20)
    return 0;
  else 
    return 1;
}

void Cal_Fast_Species(double H, double NT, double K3, double Kw, double CT, double K1, double K2,
                      double& NH3, double& NH4, double& H2CO3, double& HCO3, double& CO3){

  double tmp0, tmp1, tmp2, tmp3;

  tmp0 = 1.0 + K3*H;
  tmp1 = 1.0 + K1/H + K1*K2/H/H;
  tmp2 = tmp1*H/K1;
  tmp3 = tmp2*H/K2;

  NH3 = NT/tmp0;
  NH4 = NT*K3*H/tmp0;

  H2CO3 = CT/tmp1;
  HCO3  = CT/tmp2;
  CO3   = CT/tmp3;
}


void Cal_Fast_Species_Node(Node_Data& H, const Node_Data& NT, double K3, double Kw, const Node_Data& CT, double K1, double K2, const Node_Data& Ca,
                           const Node_Data& Cl, const Node_Data& Cha_Sum, Node_Data& NH3, Node_Data& NH4, Node_Data& H2CO3, Node_Data& HCO3, Node_Data& CO3){

  int i, j, Nx, Ny, flag, N_iter;

  Nx = H.get_Nx();
  Ny = H.get_Ny();

  for(i = 0; i <= Nx + 1; i++)
    for(j = 0; j <= Ny; j++)
      {
        flag = Newton_Root(H[i][j], NT[i][j], K3, Kw, CT[i][j], K1, K2, Ca[i][j], Cl[i][j], Cha_Sum[i][j], N_iter);

        if(flag == 1)
	  {
	    cout << endl << "In calculation of [H+], Newton's method didn't converge after inital guess is reduced by 1e"<< N_iter << " at position i = " << i 
                 << ", j = " << j << ", faliure in Cal_Fast_Species_Node() at line " << __LINE__ << " of file " << __FILE__ << endl;
	    exit(-1);
	  }

        Cal_Fast_Species(H[i][j], NT[i][j], K3, Kw, CT[i][j], K1, K2, NH3[i][j], NH4[i][j], H2CO3[i][j], HCO3[i][j], CO3[i][j]);
      }
}


//************************************************************************************************************************
//function to calculate the sum of charges at each point, used as RHS of the fast reaction equation (Charge conservation)
void Cal_Charge_Sum(Node_Data& Cha_Sum, const Node_Data& H, const Node_Data& OH, const Node_Data& Ca, const Node_Data& Cl, 
		   const Node_Data& NH4, const Node_Data& HCO3, const Node_Data& CO3){

  Cha_Sum = Z_h*H + Z_oh*OH + Z_ca*Ca + Z_cl*Cl + Z_nh4*NH4 + Z_hco3*HCO3 + Z_co3*CO3;
       
}
//************************************************************************************************************************
//************************************************************************************************************************

//function to calculate the ionic strength at each point, used for computing activity coefficient
void Cal_Ionic_Strength(Node_Data& Ion_Stren, const Node_Data& H, const Node_Data& OH, const Node_Data& Ca, const Node_Data& Cl, 
			const Node_Data& NH4, const Node_Data& HCO3, const Node_Data& CO3){

  Ion_Stren =  0.5*(Z_h*Z_h*H + Z_oh*Z_oh*OH + Z_ca*Z_ca*Ca + Z_cl*Z_cl*Cl + Z_nh4*Z_nh4*NH4 + Z_hco3*Z_hco3*HCO3 + Z_co3*Z_co3*CO3);

}
