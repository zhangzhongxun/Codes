#include "Driver.h"

#define PI (3.14159265359)

#// define DEBUG
 
//************************************************************************************************************************

//member function of class Global_Parameter, note here the global variables are defined
//in DataStruct.cc, and they are the non-dimensionalized parameters
void Global_Parameter::To_Nondimension(const Control_Parameter& CP) {

  double rho_0, c_0, h, t_0;
  double t_rh2, t2_rh4, r_t, t2_rh2, t_h2;

  rho_0 = CP.rho_0;
  c_0 = CP.c_0;
  h = CP.h;

  U_ref = h*h*CP.Delta_P/(8.0*eta_s*CP.LL);
  t_0 = h/U_ref;
  // CP.t_0 = t_0;
  Re_s = rho_s*U_ref*h/eta_s;
  L_theta = CP.L/CP.LL;


  t_rh2  = t_0/rho_0/h/h;
  t2_rh4 = t_0*t_0/rho_0/h/h/h/h;
  r_t    = rho_0/t_0;
  t2_rh2 = t_0*t_0/rho_0/h/h;
  t_h2   = t_0/h/h;

  a_EPU = CP.a_epu*Z_electron;

  Vol_Penal = vol_penal*KB*T*t2_rh4;
  Gamma_0 = gam_0*KB*T*t2_rh4;
  Ca_phase = ca_phase*KB*T*t2_rh2;
  Gamma_1 = gam_1*KB*T*t2_rh4;
  Gamma_2 = gam_2*KB*T*t2_rh2;
  Kai_CH = kai_CH;
  NP_inv = np_inv;
  PHI_reg = phi_reg;
  PHIS_reg = phis_reg;
  S_crit = s_crit;
  K_acid_1 = kacid_1;
  K_acid_2 = kacid_2;
  K_nh = k_nh;
  K_w = k_w;
  K_so = k_so;
  Ca_vol_coef = ca_vol_coef;
  Lambda_c = lam_c*r_t;
  Lambda_b = lam_b*r_t;
  Lambda_s = lam_s*r_t;
  sub_epsilon = s_eps;
  sub_mu = s_mu*t_0;
  sub_K_c = s_K_c/c_0;
  sub_A = s_A*t_0;
  k_ur = kur*t_0;
  k_p = kp*t_0;
  Ds = ds*t_h2;
  D_ur = d_ur*t_h2;
  D_nh3 = d_nh3*t_h2;
  D_nh4 = d_nh4*t_h2;
  D_h2co3 = d_h2co3*t_h2;
  D_hco3 = d_hco3*t_h2;
  D_co3 = d_co3*t_h2;
  D_ca = d_ca*t_h2;
  D_cl = d_cl*t_h2;
  D_h = d_h*t_h2;
  D_oh = d_oh*t_h2;
  D_cahco3 = d_cahco3*t_h2;
  D_caco3 = d_caco3*t_h2;
  eta_cac = eta_c*t_rh2;
  eta_bio = eta_b*t_rh2;
  eta_sol = eta_s*t_rh2;
  rho_cac = rho_c/rho_0;
  rho_bio = rho_b/rho_0;
  rho_sol = rho_s/rho_0;
  eta_ave = .5*eta_cac+ .1*eta_bio + .4*eta_sol;  //average viscosity


  //*****************debug**********************
  cout << endl << " Gamma_0 = " <<  Gamma_0  << endl ; 
 //*****************debug**********************
}

//************************************************************************************************************************

//function to read in the global physical parameters, such Gamma_1, Gamma_2, Kai, eta_n, eta_s, rho_n, rho_s, ....
void Read_Global_Parameters(char* filename, Global_Parameter& GP, const Control_Parameter& CP){

  ifstream ifile(filename, ios::in);

  if( !ifile)
    {
      cout << " \n\t Can not open the global parameters file " << filename << " !! " << endl;
      exit(-1);
    }

  char datacomment[200];

  ifile  >> datacomment >> GP.KB         //KB
         >> datacomment >> GP.T          //Temperature
         >> datacomment >> GP.gam_0      //Gamma_0      
         >> datacomment >> GP.gam_1      //Gamma_1
         >> datacomment >> GP.gam_2      //Gamma_2
         >> datacomment >> GP.kai_CH     //Kai_CH
         >> datacomment >> GP.np_inv     //NP_in
         >> datacomment >> GP.phi_reg    //PHI_reg
         >> datacomment >> GP.phis_reg   //PHIS_reg
         >> datacomment >> GP.lam_c      //Lambda_c
         >> datacomment >> GP.lam_b      //Lambda_b
         >> datacomment >> GP.lam_s      //Lambda_s
         >> datacomment >> GP.s_eps      //epsilon
         >> datacomment >> GP.s_mu       //sub_mu
         >> datacomment >> GP.s_K_c      //sub_K_c
         >> datacomment >> GP.s_A        //sub_A
         >> datacomment >> GP.ds         //Ds
         >> datacomment >> GP.d_ur       //D_ur
         >> datacomment >> GP.d_nh3      //D_nh3     
         >> datacomment >> GP.d_nh4      //D_nh4
         >> datacomment >> GP.d_h2co3    //D_h2co3
         >> datacomment >> GP.d_hco3     //D_hco3
         >> datacomment >> GP.d_co3      //D_co3     
         >> datacomment >> GP.d_ca       //D_ca
         >> datacomment >> GP.d_cl       //D_cl
         >> datacomment >> GP.d_h        //D_h
         >> datacomment >> GP.d_oh       //D_oh
	 >>   datacomment >> GP.d_cahco3 //D_cahco3                ** need to modify input file here **
	 >>   datacomment >> GP.d_caco3 //D_caco3                ** need to modify input file here **
         >> datacomment >> GP.eta_c      //eta_cac      
         >> datacomment >> GP.eta_b      //eta_bio
         >> datacomment >> GP.eta_s      //eta_sol
         >> datacomment >> GP.rho_c      //rho_cac      
         >> datacomment >> GP.rho_b      //rho_bio
         >> datacomment >> GP.rho_s      //rho_sol
         >> datacomment >> GP.kur        //k_ur
         >> datacomment >> GP.kp         //k_p
         >> datacomment >> GP.s_crit     //S_crit
         >> datacomment >> GP.kacid_1    //K_acid_1
         >> datacomment >> GP.kacid_2    //K_acid_2    
         >> datacomment >> GP.k_nh       //K_nh   
         >> datacomment >> GP.k_w        //K_w     
         >> datacomment >> GP.k_so       //K_so
         >> datacomment >> A_Davies   
         >> datacomment >> B_Davies
         >> datacomment >> alpha_Davies
         >> datacomment >> GP.ca_phase   //Ca_phase
         >> datacomment >> GP.ca_vol_coef //Ca_vol_coef
         >> datacomment >> GP.vol_penal;  //Vol_Penal   



  ifile.close();

  //Nondimensionalize the global parameters
  GP.To_Nondimension(CP);
}

//************************************************************************************************************************

//function to display the dimensional  global physical parameters, such Gamma_1, Gamma_2, Kai, eta_n, eta_s, rho_n,
// rho_s, i.e., the raw parameters read in ....
void Global_Parameter::display(char* filename) const {

  ofstream ofile(filename,ios::out);
  if(!ofile)
    {
      cout << endl << " Can not open output file " << filename <<", ERROR in  Global_Parameter::display, "
           << " Drive_ConstRho_2time.cpp !!!!!!!!!!!" << endl << endl;
      exit(-1);
    }

  ofile << endl << endl << "***************************** The Raw Global Parameters ******************* " << endl << endl;

  ofile   << endl << " KB = " << KB            //KB
	  << endl << " T = " <<  T             //Temperature
	  << endl << " gam_0 = " <<  gam_0     // Gamma_0      
	  << endl << " gam_1 = " <<  gam_1     // Gamma_1
	  << endl << " gam_2 = " << gam_2      //Gamma_2
	  << endl << " kai_CH = " <<  kai_CH   //Kai_CH
	  << endl << " np_inv = " <<  np_inv   //NP_in
	  << endl << " phi_reg = " <<  phi_reg //PHI_reg
	  << endl << " phis_reg = " << phis_reg //PHIS_reg
	  << endl << " lam_c = " <<  lam_c     //Lambda_c
	  << endl << " lam_b = " <<  lam_b     //Lambda_b
	  << endl << " lam_s = " <<  lam_s     //Lambda_s      
	  << endl << " s_eps = " <<  s_eps     //epsilon
	  << endl << " s_mu = " <<  s_mu       //sub_mu
	  << endl << " s_K_c = " <<  s_K_c     //sub_K_c
	  << endl << " s_A = " <<  s_A         //sub_A
	  << endl << " ds = " <<  ds           //Ds
	  << endl << " d_ur = " <<  d_ur       //D_ur
	  << endl << " d_nh3 = " <<  d_nh3     //D_nh3
	  << endl << " d_nh4 = " <<  d_nh4     //D_nh4      
	  << endl << " d_h2co3 = " <<  d_h2co3 //D_h2co3
	  << endl << " d_hco3 = " <<  d_hco3   //D_hco3
	  << endl << " d_co3 = " <<  d_co3     //D_co3
	  << endl << " d_ca = " <<  d_ca       //D_ca      
	  << endl << " d_cl = " <<  d_cl       //D_cl
	  << endl << " d_h = " <<  d_h         //D_h      
	  << endl << " d_oh = " <<  d_oh       //D_oh
	  << endl << " d_cahco3  = " <<  d_cahco3         //D_cahco3     
	  << endl << " d_caco3  = " <<  d_caco3      //D_caco3
      	  << endl << " eta_c = " << eta_c      //eta_cac
	  << endl << " eta_b = " << eta_b      //eta_bio
	  << endl << " eta_s = " <<  eta_s     //eta_sol
      	  << endl << " rho_c = " <<  rho_c     //rho_cac
	  << endl << " rho_b = " <<  rho_b     //rho_bio
	  << endl << " rho_s = " <<  rho_s     //rho_sol
	  << endl << " kur = " << kur
	  << endl << " kp = "  << kp
	  << endl << " s_crit = " << s_crit
	  << endl << " kacid_1 = " << kacid_1
	  << endl << " kacid_2 = " << kacid_2
	  << endl << " k_nh = " << k_nh
	  << endl << " k_w = "  << k_w
	  << endl << " k_so = " << k_so
	  << endl << " A_Davies = " << A_Davies
	  << endl << " B_Davies = " << B_Davies
	  << endl << " alpha_Davies = " << alpha_Davies
	  << endl << " ca_phase = " << ca_phase
	  << endl << " ca_vol_coef = " << ca_vol_coef
	  << endl << " vol_penal = " << vol_penal;
      
  ofile.close();
}

 
//************************************************************************************************************************

//function to read in the control parameters such as dx, dy , Nx, k (m = 2^(k+1), Ny = m )
//total time T_final, CFL, tolerance of GMRES, ...
void Read_Control_Parameters(char* filename, Control_Parameter& CP){

  ifstream ifile(filename, ios::in);

  if( !ifile)
    {
      cout << " \n\t Can not open the control parameters file " << filename << " !! " << endl;
      exit(-1);
    }

  char datacomment[200];

  ifile  >> datacomment >> CP.dx
         >> datacomment >> CP.dy
         >> datacomment >> CP.Nx
         >> datacomment >> CP.k
         >> datacomment >> CP.T_final
         >> datacomment >> CP.CFL
         >> datacomment >> Newton_tol
         >> datacomment >> Newton_max_k
         >> datacomment >> Activity_max_k
         >> datacomment >> Activity_tol
         >> datacomment >> CP.GMRES_tol
         >> datacomment >> CP.GMRES_max_k
         >> datacomment >> CP.Init_flag
         >> datacomment >> CP.BC_flag
         >> datacomment >> CP.velo_flag
         >> datacomment >> CP.Grid_flag
         >> datacomment >> CP.Init_N      
         >> datacomment >> CP.rho_0
         >> datacomment >> CP.c_0
         >> datacomment >> CP.h
         >> datacomment >> CP.L
         >> datacomment >> CP.LL
         >> datacomment >> CP.Delta_P
         >> datacomment >> CP.F_Stride
         >> datacomment >> CP.tfac
         >> datacomment >> CP.alph
         >> datacomment >> CP.Ini_phc
         >> datacomment >> CP.Ini_phb
         >> datacomment >> CP.Ini_Ur
         >> datacomment >> CP.Ini_NH3
         >> datacomment >> CP.Ini_NH4
         >> datacomment >> CP.Ini_H2CO3
         >> datacomment >> CP.Ini_HCO3
         >> datacomment >> CP.Ini_CO3
         >> datacomment >> CP.Ini_Ca
         >> datacomment >> CP.Ini_Cl
         >> datacomment >> CP.Ini_H
         >> datacomment >> CP.Ini_OH
	 >> datacomment >> CP.Ini_CaHCO3     // ** need to modify input file **
         >> datacomment >> CP.Ini_CaCO3        // ** need to modify input file **
         >> datacomment >> CP.Ini_c
         >> datacomment >> CP.cutoff_c_time
         >> datacomment >> CP.Shear_v
         >> datacomment >> CP.D_flag
         >> datacomment >> CP.D_phic_tol
         >> datacomment >> CP.D_delta
         >> datacomment >> CP.a_epu
         >> datacomment >> CP.EP_flag
         >> datacomment >> CP.eps_dif
         >> datacomment >> CP.N_dif_steps
         >> datacomment >> CP.Difu_Init_flag
         >> datacomment >> CP.Ion_Act_flag
         >> datacomment >> SaveDir
         >> datacomment >> chemc_flag
         >> datacomment >> CP.Ave_Vis_flag;


  //*****************debug**********************
  cout << endl << " CP.a_epu = " <<  CP. a_epu  << endl ; 
 //*****************debug**********************

  ifile.close();
}

//************************************************************************************************************************

//function to display the control parameters such as dx, dy , Nx, k (m = 2^(k+1), Ny = m )
//total time T_final, CFL, tolerance of GMRES, ...
void Control_Parameter::display(char* filename) const {

  ofstream ofile(filename,ios::out);
  if(!ofile)
    {
      cout << endl << " Can not open output file " << filename <<", ERROR in  Control_Parameter::display, "
           << " Drive_ConstRho_2time.cpp !!!!!!!!!!!" << endl << endl;
      exit(-1);
    }

  ofile << endl << endl << "********************* The Raw Control Parameters ***************************" << endl << endl;

  ofile  << endl <<  " dx = " << dx
         << endl <<  " dy = " << dy
         << endl <<  " Nx = " << Nx
         << endl <<  " k  = " << k
         << endl <<  " Tfinal = " << T_final
         << endl <<  " CFL = " << CFL
         << endl <<  " Newton_tol = " <<  Newton_tol
         << endl <<  " Newton_max_k = " <<  Newton_max_k
         << endl <<  " Activity_max_k = " <<  Activity_max_k
         << endl <<  " Activity_tol = " << Activity_tol
         << endl <<  " GMRES_tol = " << GMRES_tol
         << endl <<  " GMRES_max_k = " << GMRES_max_k
         << endl <<  " Init_flag = " << Init_flag
         << endl <<  " BC_flag = " << BC_flag
         << endl <<  " velo_flag = " << velo_flag
         << endl <<  " Grid_flag = " << Grid_flag
         << endl <<  " Init_N = " << Init_N
         << endl <<  " rho_0 = " << rho_0
         << endl <<  " c_0 = " << c_0
         << endl <<  " h = " << h
         << endl <<  " U_ref = " << U_ref
         << endl <<  " L = " << L
         << endl <<  " LL = " << LL
         << endl <<  " Delta_P = " << Delta_P
         << endl <<  " F_Stride = " << F_Stride
         << endl <<  " Tfac = " << tfac
         << endl <<  " alpha = " << alph
         << endl <<  " Ini_phc = " << Ini_phc
         << endl <<  " Ini_phb = " << Ini_phb
         << endl <<  " Ini_Ur = " << Ini_Ur
         << endl <<  " Ini_NH3 = " << Ini_NH3
         << endl <<  " Ini_NH4 = " << Ini_NH4
         << endl <<  " Ini_H2CO3 = " << Ini_H2CO3
         << endl <<  " Ini_HCO3 = " << Ini_HCO3
         << endl <<  " Ini_CO3 = " << Ini_CO3 
         << endl <<  " Ini_Ca = " << Ini_Ca
         << endl <<  " Ini_Cl = " << Ini_Cl
         << endl <<  " Ini_H = " << Ini_H
         << endl <<  " Ini_OH = " << Ini_OH
         << endl <<  " Ini_c = " << Ini_c 
         << endl <<  " cutoff_c_time = " << cutoff_c_time
         << endl <<  " Shear_velocity = " << Shear_v
         << endl <<  " D_flag = " << D_flag
         << endl <<  " D_phic_tol = " << D_phic_tol
         << endl <<  " D_delta = " << D_delta
         << endl <<  " a_epu = " << a_epu
         << endl <<  " EP_flag = " << EP_flag
         << endl <<  " eps_dif = " << eps_dif
         << endl <<  " N_dif_steps = " << N_dif_steps
         << endl <<  " Diffusion const of Initial = " << Difu_Init_flag
         << endl <<  " Ion_Act_flag = " << Ion_Act_flag
         << endl <<  " Save directory = " << SaveDir
         << endl <<  " chemc_flag = " << chemc_flag
         << endl <<  " Ave_viscosity_flag = " << Ave_Vis_flag;

  ofile.close();
}


//************************************************************************************************************************

//function to save the intermediate data for later restart, so save the result at time step n and n-1,
//also need to save the time step size at n and n-1 for the use in extrapolation, and the current time
void Save_Inter_Data(char* filename, const Node_Data& un, const Node_Data& vn, const Node_Data& p, const Node_Data& s, double* s_stg, const Node_Data& phi_c, 
                     const Node_Data& phi_b, const Node_Data& Ur, const Node_Data& NH3, const Node_Data& NH4, const Node_Data& H2CO3, const Node_Data& HCO3, 
                     const Node_Data& CO3, const Node_Data& Ca, const Node_Data& Cl, const Node_Data& H, Node_Data& OH, const Node_Data& c, 
                     double T_current, int Inter_file_NO, int n){

  ofstream ofile(filename, ios::out);

  if(!ofile)
    {
        cout << endl << "Can not open output file " << filename << ", ERROR in Save_Inter_Data() , line " << __LINE__ << "of file " << __FILE__ << endl;
        exit(-1);
    }

  ofile << un.get_Nx() << endl << un.get_Ny() << endl << T_current << endl <<  Inter_file_NO << endl << n << endl;
  ofile << un << endl;
  ofile << vn << endl;
  ofile << p << endl;
  ofile << s << endl;
  ofile << phi_c << endl;
  ofile << phi_b << endl;
  ofile << Ur << endl;
  ofile << NH3 << endl;
  ofile << NH4 << endl;
  ofile << H2CO3 << endl;
  ofile << HCO3 << endl;
  ofile << CO3 << endl;
  ofile << Ca << endl;
  ofile << Cl << endl;
  ofile << H << endl;
  ofile << OH << endl;
  ofile << c << endl;


  int i, j, Nx, Ny;
  Nx = un.get_Nx();
  Ny = un.get_Ny();

  //output the pressure for half-staggered grid
  ofile << endl;
  for(j = 0; j <= Ny-1; j++)
    {
      for(i = 0; i <= Nx; i++)
	ofile << s_stg[j*(Nx+1)+i] << "  ";
      ofile << endl;
    } 

  ofile.close();
}

//************************************************************************************************************************

//function to read in the saved intermediate data
void Read_Inter_Data(char* filename,  size_t Nx, size_t Ny, Node_Data& un,  Node_Data& vn,  Node_Data& p, Node_Data& s, double* s_stg, Node_Data& phi_c,
                     Node_Data& phi_b, Node_Data& Ur, Node_Data& NH3, Node_Data& NH4, Node_Data& H2CO3, Node_Data& HCO3, Node_Data& CO3, Node_Data& Ca,
		     Node_Data& Cl, Node_Data& H, Node_Data& OH, Node_Data& c, double& T_current, int& Inter_file_NO, int& n){

  ifstream ifile(filename, ios::in);

  if(!ifile)
    {
      cout << endl << "Can not open input file " << filename << ", ERROR in Read_Inter_Data() at line " << __LINE__ << " of file " << __FILE__  << endl;
      exit(-1);
    }

  size_t nx, ny;

  ifile >> nx >> ny;
  if(nx != Nx || ny != Ny)
    {
      cout << endl << "Read in dimension incompatible with the internal dimension. ERROR in Read_Inter_Data() in " 
           << " Driver.cc !! " << endl;
      exit(-1);
    }

  ifile >> T_current >> Inter_file_NO >> n;

  ifile >> un >> vn >>  p >> s >> phi_c >> phi_b >> Ur >> NH3 >> NH4 >> H2CO3 >> HCO3 >> CO3 >> Ca >> Cl >> H >> OH >> c;

  //read  the pressure for half-staggered grid

  for(int j = 0; j <= Ny-1; j++)
      for(int i = 0; i <= Nx; i++)
	ifile >> s_stg[j*(Nx+1)+i];
 

  ifile.close();
}

//************************************************************************************************************************

//function to setup the initial conditinos, use flag to indicate the method to set the initial condition
//flag >= 0, different value corresponding to different initial data
void Set_Initial_Condition(char fname[][80], Node_Data& u, Node_Data& v, Node_Data& phi_c, Node_Data& phi_b, Node_Data& Ur, Node_Data& NH3, Node_Data& NH4,
                           Node_Data& H2CO3, Node_Data& HCO3, Node_Data& CO3, Node_Data& Ca,  Node_Data& Cl, Node_Data& H, Node_Data& OH, Node_Data& c, 
                           double dx, double dy, const Control_Parameter& CP){

  int i, j, Nx, Ny, seed, N, flag;
  double x, y, x_bar, y_bar, d, a, tmp, dNx, Ini_phc, Ini_phb, Ini_Ur, Ini_NH3, Ini_NH4, Ini_H2CO3, Ini_HCO3, Ini_CO3, 
         Ini_Ca, Ini_Cl, Ini_H, Ini_OH, Ini_c, xL, hxL, yL, hyL, dt, decay_x, decay_y;

  N = CP.Init_N;

  double ran_r[N];   //radius of random islands for initial phi_n
  double ran_x[N];   //x-coordinate of the center of the random islands
  double ran_y[N];   //y-coordinates ..................................
  double ran_phb[3*N]; //height of the random islands

  double xbL[3];  //left of the bump
  double xbR[3];  //right
  double ybB[3];  //bottom
  double ybT[3];  //top
  double ybF[3]; //Flat
  double ySign[3];

  Nx = u.get_Nx();
  Ny = u.get_Ny();

  flag = CP.Init_flag;
  
  Ini_phc = CP.Ini_phc;
  Ini_phb = CP.Ini_phb;
  Ini_Ur = CP.Ini_Ur;
  Ini_NH3 = CP.Ini_NH3;
  Ini_NH4 = CP.Ini_NH4;
  Ini_H2CO3 = CP.Ini_H2CO3;
  Ini_HCO3 = CP.Ini_HCO3;
  Ini_CO3 = CP.Ini_CO3;
  Ini_Ca = CP.Ini_Ca;
  Ini_Cl = CP.Ini_Cl;
  Ini_H = CP.Ini_H;  
  Ini_OH = CP.Ini_OH;
  Ini_c = CP.Ini_c;

  dt = dx*CP.tfac;

  
  if(flag == 100) //read all variables from file except velocity
    {
      ifstream ifile;

      ifile.open(fname[0],ios::in);
      if(!ifile)
          {
              cout << endl << " Can not open input file " << fname[0] << ", ERROR in Set_Initial_Condition() at line " << __LINE__ << " of file " << __FILE__ << endl;
              exit(-1);
          }
      ifile >> phi_c;
      ifile.close();

      ifile.open(fname[1],ios::in);
      if(!ifile)
          {
              cout << endl << " Can not open input file " << fname[1] << ", ERROR in Set_Initial_Condition() at line " << __LINE__ << " of file " << __FILE__ << endl;
              exit(-1);
          }
      ifile >> phi_b;
      ifile.close();

      
      ifile.open(fname[2],ios::in);
      if(!ifile)
          {
              cout << endl << " Can not open input file " << fname[2] << ", ERROR in Set_Initial_Condition() at line " << __LINE__ << " of file " << __FILE__ << endl;
              exit(-1);
          }
      ifile >> Ur;
      ifile.close();

      ifile.open(fname[3],ios::in);
      if(!ifile)
          {
              cout << endl << " Can not open input file " << fname[3] << ", ERROR in Set_Initial_Condition() at line " << __LINE__ << " of file " << __FILE__ << endl;
              exit(-1);
          }
      ifile >> NH3;
      ifile.close();

      ifile.open(fname[4],ios::in);
      if(!ifile)
          {
              cout << endl << " Can not open input file " << fname[4] << ", ERROR in Set_Initial_Condition() at line " << __LINE__ << " of file " << __FILE__ << endl;
              exit(-1);
          }
      ifile >> NH4;
      ifile.close();

      
      ifile.open(fname[5],ios::in);
      if(!ifile)
          {
              cout << endl << " Can not open input file " << fname[5] << ", ERROR in Set_Initial_Condition() at line " << __LINE__ << " of file " << __FILE__ << endl;
              exit(-1);
          }
      ifile >> H2CO3;
      ifile.close();

      ifile.open(fname[6],ios::in);
      if(!ifile)
          {
              cout << endl << " Can not open input file " << fname[6] << ", ERROR in Set_Initial_Condition() at line " << __LINE__ << " of file " << __FILE__ << endl;
              exit(-1);
          }
      ifile >> HCO3;
      ifile.close();

      ifile.open(fname[7],ios::in);
      if(!ifile)
          {
              cout << endl << " Can not open input file " << fname[7] << ", ERROR in Set_Initial_Condition() at line " << __LINE__ << " of file " << __FILE__ << endl;
              exit(-1);
          }
      ifile >> CO3;
      ifile.close();

      ifile.open(fname[8],ios::in);
      if(!ifile)
          {
              cout << endl << " Can not open input file " << fname[8] << ", ERROR in Set_Initial_Condition() at line " << __LINE__ << " of file " << __FILE__ << endl;
              exit(-1);
          }
      ifile >> Ca;
      ifile.close();

      ifile.open(fname[9],ios::in);
      if(!ifile)
          {
              cout << endl << " Can not open input file " << fname[9] << ", ERROR in Set_Initial_Condition() at line " << __LINE__ << " of file " << __FILE__ << endl;
              exit(-1);
          }
      ifile >> Cl;
      ifile.close();

      
      ifile.open(fname[10],ios::in);
      if(!ifile)
          {
              cout << endl << " Can not open input file " << fname[10] << ", ERROR in Set_Initial_Condition() at line " << __LINE__ << " of file " << __FILE__ << endl;
              exit(-1);
          }
      ifile >> H;
      ifile.close();

      ifile.open(fname[11],ios::in);
      if(!ifile)
          {
              cout << endl << " Can not open input file " << fname[11] << ", ERROR in Set_Initial_Condition() at line " << __LINE__ << " of file " << __FILE__ << endl;
              exit(-1);
          }
      ifile >> c;

      ifile.close();

      ifile.open(fname[12],ios::in);
      if(!ifile)
          {
              cout << endl << " Can not open input file " << fname[12] << ", ERROR in Set_Initial_Condition() at line " << __LINE__ << " of file " << __FILE__ << endl;
              exit(-1);
          }
      ifile >> OH;
      ifile.close();
    }
      
  //get the random coefficient
  xL = dx*(Nx+1);  //the x-length of domain
  hxL = .5*xL;     //half size of x-length
  yL = dy*Ny;

  seed = time(0);
  srand(seed);
//   for(i = 0; i <= 2*N/3; i++)
//     {
//       ran_r[i] = 0.03*rand()/double(RAND_MAX) + 0.02;
//       ran_x[i] = .4*xL*rand()/double(RAND_MAX) + 0.1*xL; 
//       ran_y[i] = 0.8*rand()/double(RAND_MAX) + 0.1;
//       ran_phb[i] = 0.05*rand()/double(RAND_MAX) + 0.15;
//     }

//   for(i =  2*N/3 + 1; i < N; i++)
//     {
//       ran_r[i] = 0.03*rand()/double(RAND_MAX) + 0.02;
//       ran_x[i] = .4*xL*rand()/double(RAND_MAX) + 0.5*xL; 
//       ran_y[i] = 0.8*rand()/double(RAND_MAX) + 0.1;
//       ran_phb[i] = 0.05*rand()/double(RAND_MAX) + 0.15;
//     }

  xbL[0] = xL/4.0;
  xbR[0] = 3.0*xL/8.0;
  ybF[0] = 0.0;
  ySign[0] = 1.0;
  
  xbL[1] = 4.0*xL/8.0;
  xbR[1] = 3.0*xL/4.0;
  ybF[1] = 0.0;
  ySign[1] = 1.0;

  xbL[2] = 3.0*xL/8.0;
  xbR[2] = 3.0*xL/4.0;
  ybF[2] = yL;
  ySign[2] = -1.0; 


  seed = time(0);
  srand(seed);

  for(i = 0; i < 3*N; i++)
    ran_phb[i] = (2*rand()/double(RAND_MAX) - 1); // *(1 - N/double(2*N-2) + (N-i)/double(2*N-2));

  dNx = xL/double(N); // 1.0/double(N);

  for(j = 0; j <= Ny; j++)
    for(i = 0; i <= Nx + 1; i++)
      {
	x = i*dx;
	y = j*dy;     

	decay_x = .5*(tanh(40*(1.0/8 - x)) + 1); 

        x_bar = ((x - xL/8)/(.75*xL))*(2*PI);

        tmp = ((x - xL/8) < (7.0*xL/8 - x)) ? (x - xL/8) : (7.0*xL/8 - x);

	y_bar = .3 + .15*tanh(tmp);

        for(int iN = 0; iN < N; iN++)
          y_bar += .05*ran_phb[iN]*sin((iN+1)*x_bar);
// 	  y_bar *= (1 + .2*ran_phb[iN]*sin((iN+1)*x_bar));

        tmp = ((y < yL/2)  ? y : (yL - y));

        decay_y = tanh(10*tmp);

	u[i][j] = CP.Shear_v*(4.0/yL/yL)*y*(yL-y); // *dt;
	v[i][j] = 0.0;

	if(flag == 0)  //isolated islands, random positions, random heights of biofilm
	  {
	    phi_c[i][j] = Ini_phc;
	    Ur[i][j] = Ini_Ur*decay_x; //Ini_Ur;
	    NH3[i][j] = Ini_NH3*decay_x;
	    NH4[i][j] = Ini_NH4*decay_x;
	    H2CO3[i][j] = Ini_H2CO3*decay_x;
	    HCO3[i][j] = Ini_HCO3*decay_x;
	    CO3[i][j] = Ini_CO3*decay_x;
	    Ca[i][j] = Ini_Ca*decay_x;
	    Cl[i][j] = Ini_Cl*decay_x;
	    H[i][j] = Ini_H;
            OH[i][j] = Ini_OH;
	    c[i][j] = Ini_c;
                
	    phi_b[i][j] = 0.0;

	    if(x_bar >= 0 && x_bar <= 2*PI && y <= y_bar)
	      phi_b[i][j] = Ini_phb; //*(1 - .5*y);

	  }
	else if(flag == 1) //constant in the whole domain for chemical species
	  {
	    phi_c[i][j] = Ini_phc;
	    Ur[i][j] = Ini_Ur; 
	    NH3[i][j] = Ini_NH3;
	    NH4[i][j] = Ini_NH4;
	    H2CO3[i][j] = Ini_H2CO3;
	    HCO3[i][j] = Ini_HCO3;
	    CO3[i][j] = Ini_CO3;
	    Ca[i][j] = Ini_Ca;
	    Cl[i][j] = Ini_Cl;
	    H[i][j] = Ini_H;
            OH[i][j] = Ini_OH;
	    c[i][j] = Ini_c;
                
	    phi_b[i][j] =  Ini_phb; //0.0;

// 	    if(x_bar >= 0 && x_bar <= 2*PI && y <= y_bar)
// 	      phi_b[i][j] = Ini_phb*(1 - .5*y);
	  }
	else if(flag == 2) //constant in the whole domain for chemical species
	  {
	    phi_c[i][j] = Ini_phc;
	    Ur[i][j] = Ini_Ur; 
	    NH3[i][j] = Ini_NH3;
	    NH4[i][j] = Ini_NH4;
	    H2CO3[i][j] = Ini_H2CO3;
	    HCO3[i][j] = Ini_HCO3;
	    CO3[i][j] = Ini_CO3;
	    Ca[i][j] = Ini_Ca;
	    Cl[i][j] = Ini_Cl;
	    H[i][j] = Ini_H;
            OH[i][j] = Ini_OH;
	    c[i][j] = Ini_c;
                
	    phi_b[i][j] = 0.0;  //semi circle on top and bottom boundary  as biofilm contour

	    if( ((x - xL/2)*(x - xL/2) + (y - yL)*(y - yL)) <= 0.04 || ((x - xL/2)*(x - xL/2) + (y - 0)*(y - 0)) <= 0.04 )
	      phi_b[i][j] = Ini_phb;
	  }
	else if(flag == 3) //2 bumps at bottom and 1 at the top
	  {
	    phi_c[i][j] = Ini_phc;
	    Ur[i][j] = Ini_Ur; 
	    NH3[i][j] = Ini_NH3;
	    NH4[i][j] = Ini_NH4;
	    H2CO3[i][j] = Ini_H2CO3;
	    HCO3[i][j] = Ini_HCO3;
	    CO3[i][j] = Ini_CO3;
	    Ca[i][j] = Ini_Ca;
	    Cl[i][j] = Ini_Cl;
	    H[i][j] = Ini_H;
            OH[i][j] = Ini_OH;
	    c[i][j] = Ini_c;
                
	    phi_b[i][j] = 0.0;

	    for(int ib = 0; ib < 3; ib++)
	      {
  		x_bar = ((x - xbL[ib])/(xbR[ib] - xbL[ib]))*(2*PI);

		tmp = ((x - xbL[ib]) < (xbR[ib] - x)) ? (x - xbL[ib]) : (xbR[ib] - x);

		y_bar = .1 + .15*tanh(tmp);

// 		seed = time(0);
// 		srand(seed);

		for(int iN = 0; iN < N; iN++)
		  {
// 		    ran_phb[iN] = (2*rand()/double(RAND_MAX) - 1)*(1 - N/double(2*N-2) + (N-iN)/double(2*N-2));
		    y_bar += .05*ran_phb[ib*N + iN]*(1 - N/double(2*N-2) + (N-iN)/double(2*N-2))*sin((iN+1)*x_bar);
		  }

                if(ySign[ib] > 0)
		  {
		    ybB[ib] = ybF[ib];
                    ybT[ib] = ybF[ib] + y_bar;
		  }
		else if(ySign[ib] < 0)
		  {
                    ybB[ib] = ybF[ib] - y_bar;
                    ybT[ib] = ybF[ib];
		  } 
		if(x >= xbL[ib] && x <= xbR[ib] && y >= ybB[ib] && y <= ybT[ib])
		  {
		    phi_b[i][j] = Ini_phb;
                    break;
		  }
	      }
	  }
 	else if(flag == 100) //read from file
          {
            tmp = phi_b.max();
       
            phi_b = (Ini_phb/tmp)*phi_b;
	    continue;
	  }
  	else
	  {
	    cout << endl << "The initial flag is not in valid range, ERROR in Set_Initial_Condition() at line " << __LINE__ << " of file " << __FILE__ << endl;

	    exit(-1);
	  }
      }
}

//************************************************************************************************************************

void Diffuse_Initial_Condition(Node_Data& phi, int N_steps, double dx, double dy, double dt, double alpha, double eps_dif){

  int i, Nx, Ny, GMRES_flag, GMRES_max_k, GMRES_num_restart, GMRES_num_k;

  double GMRES_tol, GMRES_res;

  GMRES_max_k = 30;
  GMRES_tol = 1e-8;

  Nx = phi.get_Nx();
  Ny = phi.get_Ny();

  c_NoFlux_x_NoFlux_y_mat A(Nx, Ny);

  valarray<double> F(A.get_dim());
  valarray<double> rhs(F.size());

  for(i = 0; i < N_steps; i++)
    {
      Build_Diffusion(A, F, phi, dx, dy, dt, alpha, eps_dif);

      //solve for new phi by GMRES
      GMRES(A, F, rhs, 1e-12, GMRES_max_k, GMRES_flag, GMRES_res, GMRES_num_restart, GMRES_num_k);
     
      //update phi after one step of diffusion
      Update_Node_Data_N_x_N_y(phi, rhs);
    }

}

//************************************************************************************************************************

//function to update a Node_Data variable with periodic BC in x-direction and, Dirichelet BC  in y-direction,
//should be used for u, v and phi: the velocity and the network concentration

void Update_Node_Data_P_x_D_y(Node_Data& u, const valarray<double>& x){

  int i, j, Nx, Ny;

  Nx = u.get_Nx();
  Ny = u.get_Ny();

  if((Nx+1)*(Ny-1) != x.size())
    {
      cout << endl << " Dimension Incompatible in Update_Node_Data_P_x_D_y(), ERROR in Driver.cc !!!!!!!!!!! " << endl;
      exit(-1);
    }

  for(j = 1; j <= Ny - 1; j++)
    {
      for(i = 0; i <= Nx; i++)
	u[i][j] = x[(j-1)*(Nx+1) + i];
   
      //set value at i = Nx + 1 equal to those at i = 0
      u[Nx+1][j] = u[0][j];
    }
}

//************************************************************************************************************************

void Update_Node_Data_D_x_D_y(Node_Data& u, const valarray<double>& x){

  int i, j, Nx, Ny;

  Nx = u.get_Nx();
  Ny = u.get_Ny();

  if(Nx*(Ny-1) != x.size())
    {
      cout << endl << " Dimension Incompatible in Update_Node_Data_D_x_D_y(), "
           << " ERROR in Driver.cc !!!!!!!!!!! " << endl;
      exit(-1);
    }

  for(j = 1; j <= Ny - 1; j++)
    {
     for(i = 1; i <= Nx; i++)
       u[i][j] = x[(j-1)*Nx + i - 1];
    }
}

//************************************************************************************************************************

void Update_Node_Data_OutFlow_x_D_y(Node_Data& u, const valarray<double>& x){

  int i, j, Nx, Ny;

  Nx = u.get_Nx();
  Ny = u.get_Ny();

  if(Nx*(Ny-1) != x.size())
    {
      cout << endl << " Dimension Incompatible in Update_Node_Data_OutFlow_x_D_y(), "
           << " ERROR in Driver.cc !!!!!!!!!!! " << endl;
      exit(-1);
    }

  for(j = 1; j <= Ny - 1; j++)
    {
     for(i = 1; i <= Nx; i++)
       u[i][j] = x[(j-1)*Nx + i - 1];

     u[0][j] = u[1][j];     //2.0*u[1][j] - u[2][j];
     u[Nx+1][j] = u[Nx][j]; //2.0*u[Nx][j] - u[Nx-1][j];
    }
}
//************************************************************************************************************************

void Update_Node_Data_P_x_N_y(Node_Data& u, const valarray<double>& x){

  int i, j, Nx, Ny;

  Nx = u.get_Nx();
  Ny = u.get_Ny();

  if((Nx+1)*(Ny+1) != x.size())
    {
      cout << endl << " Dimension Incompatible in Update_Node_Data_P_x_N_y(), ERROR in Driver.cc !!!!!!!!!!! " << endl;
      exit(-1);
    }

  for(j = 0; j <= Ny; j++)
    {
      for(i = 0; i <= Nx; i++)
	u[i][j] = x[j*(Nx+1) + i];

      //set value at i = Nx + 1 equal to those at i = 0
      u[Nx+1][j] = u[0][j];
    }
}

//************************************************************************************************************************
//function to update Electric potential, No-Flux BC in x and y, 1 <= i <= Nx
//or No-Flux in y, and homogeneous Dirichlet in x
void Update_EP(Node_Data& u, const valarray<double>& x){

//   int i, j, Nx, Ny;

//   Nx = u.get_Nx();
//   Ny = u.get_Ny();

//   if((Nx+2)*(Ny+1) != x.size())
//     {
//       cout << endl << " Dimension Incompatible in Update_EP(), ERROR at line " << __LINE__ << " of file "
//            << __FILE__ << endl;
//       exit(-1);
//     }

//   for(j = 0; j <= Ny; j++)
//     {
//       for(i = 0; i <= Nx+1; i++)
// 	u[i][j] = x[j*(Nx+2) + i];
//     }

  int i, j, Nx, Ny;

  Nx = u.get_Nx();
  Ny = u.get_Ny();

  if(Nx*(Ny+1) != x.size())
    {
      cout << endl << " Dimension Incompatible in Update_EP(), ERROR at line " << __LINE__ << " of file "
           << __FILE__ << endl;
      exit(-1);
    }

  for(j = 0; j <= Ny; j++)
    {
      for(i = 1; i <= Nx; i++)
	u[i][j] = x[j*Nx + (i-1)];

      u[0][j] = u[2][j];
      u[Nx+1][j] = u[Nx-1][j]; 
    }
}

//************************************************************************************************************************

void Update_Node_Data_N_x_N_y(Node_Data& u, const valarray<double>& x){

  int i, j, Nx, Ny;

  Nx = u.get_Nx();
  Ny = u.get_Ny();

  if((Nx+2)*(Ny+1) != x.size())
    {
      cout << endl << " Dimension Incompatible in Update_Node_Data_N_x_N_y(), ERROR in Driver.cc !!!!!!!!!!! " << endl;
      exit(-1);
    }

  for(j = 0; j <= Ny; j++)
    {
      for(i = 0; i <= Nx + 1; i++)
	u[i][j] = x[j*(Nx+2) + i];
    }
}

//************************************************************************************************************************

void Update_Node_Data_N_x_D_y(Node_Data& u, const valarray<double>& x){

  int i, j, Nx, Ny;

  Nx = u.get_Nx();
  Ny = u.get_Ny();

  if((Nx+2)*(Ny-1) != x.size())
    {
      cout << endl << " Dimension Incompatible in Update_Node_Data_N_x_D_y(), ERROR in Driver.cc !!!!!!!!!!! " << endl;
      exit(-1);
    }

  for(j = 1; j <= Ny-1; j++)
    {
      for(i = 0; i <= Nx + 1; i++)
	u[i][j] = x[(j-1)*(Nx+2) + i];
    }
}

//************************************************************************************************************************

void Update_Node_Data_N_x_D_N_y(Node_Data& u, const valarray<double>& x){

  int i, j, Nx, Ny;

  Nx = u.get_Nx();
  Ny = u.get_Ny();

  if((Nx+2)*Ny != x.size())
    {
      cout << endl << " Dimension Incompatible in Update_Node_Data_N_x_D_N_y(), ERROR in Driver.cc !!!!!!!!!!! " << endl;
      exit(-1);
    }

  for(j = 0; j <= Ny-1; j++)
    {
      for(i = 0; i <= Nx + 1; i++)
	u[i][j] = x[j*(Nx+2) + i];
    }
}

//************************************************************************************************************************

void Update_Node_Data_P_x_D_N_y(Node_Data& u, const valarray<double>& x){

  int i, j, Nx, Ny;

  Nx = u.get_Nx();
  Ny = u.get_Ny();

  if((Nx+1)*Ny != x.size())
    {
      cout << endl << " Dimension Incompatible in Update_Node_Data_P_x_D_N_y(), ERROR in Driver.cc !!!!!!!!!!! " << endl;
      exit(-1);
    }

  for(j = 0; j <= Ny-1; j++)
    {
      for(i = 0; i <= Nx; i++)
	u[i][j] = x[j*(Nx+1) + i];

      u[Nx+1][j] = u[0][j];
    }
}

//************************************************************************************************************************
void ADD_Update_Node_Data_P_x_D_y(Node_Data& u, const valarray<double>& x){

  int i, j, Nx, Ny;

  Nx = u.get_Nx();
  Ny = u.get_Ny();

  if((Nx+2)*(Ny-1) != x.size())
    {
      cout << endl << " Dimension Incompatible in Add_Update_Node_Data_P_x_D_y(), ERROR in Driver.cc !!!!!!!!!!! "
           << endl;
      exit(-1);
    }

  for(j = 1; j <= Ny - 1; j++)
    {
      for(i = 0; i <= Nx; i++)
	u[i][j] += x[(j-1)*(Nx+2) + i];  

      //set value at i = Nx + 1 equal to those at i = 0
      u[Nx+1][j] = u[0][j];
    }
}

//************************************************************************************************************************

void ADD_Update_Node_Data_D_x_D_y(Node_Data& u, const valarray<double>& x){

  int i, j, Nx, Ny;

  Nx = u.get_Nx();
  Ny = u.get_Ny();

  if((Nx+2)*(Ny-1) != x.size())
    {
      cout << endl << " Dimension Incompatible in Add_Update_Node_Data_D_x_D_y(), ERROR in Driver.cc !!!!!!!!!!! "
           << endl;
      exit(-1);
    }

  for(j = 1; j <= Ny - 1; j++)
    for(i = 1; i <= Nx; i++)
      u[i][j] += x[(j-1)*(Nx+2) + i];
}

//************************************************************************************************************************

void ADD_Update_Node_Data_OutFlow_x_D_y(Node_Data& u, const valarray<double>& x){

  int i, j, Nx, Ny;

  Nx = u.get_Nx();
  Ny = u.get_Ny();

  if((Nx+2)*(Ny-1) != x.size())
    {
      cout << endl << " Dimension Incompatible in Add_Update_Node_Data_OutFlow_x_D_y(), ERROR in Driver.cc !!!!!!!!!!! "
           << endl;
      exit(-1);
    }

  for(j = 1; j <= Ny - 1; j++)
    {
     for(i = 1; i <= Nx; i++)
       u[i][j] += x[(j-1)*(Nx+2) + i];

     u[0][j] = u[1][j];     //2.0*u[1][j] - u[2][j];
     u[Nx+1][j] = u[Nx][j]; //2.0*u[Nx][j] - u[Nx-1][j];
    }
}


//************************************************************************************************************************

void Update_Velo_By_Momentum(Node_Data& u, const valarray<double>& x, int bc_flag){

  if(bc_flag == 0 || bc_flag == 3) //periodic BC
    Update_Node_Data_P_x_D_y(u, x);
  else if(bc_flag == 1)  //no-flux in x and y-direction (i.e. Dirichelete BC for u, v in x-directino)
    Update_Node_Data_D_x_D_y(u, x);
  else if(bc_flag == 2)  //out-flow BC in x-direction
    Update_Node_Data_OutFlow_x_D_y(u, x);
  else
    {
      cout << endl << " bc_flag must be 0, 1 , 2, 3, ERROR in Update_Velo_By_Momentum() at Driver.cc() !!! " << endl ;
      exit(-1);
    } 
}

//************************************************************************************************************************

void Update_Velo_By_Pressure(Node_Data& u, const valarray<double>& x, int bc_flag){

  if(bc_flag == 0)  //no-flux in x-direction (i.e. Dirichelete BC for u, v in x-directino)
    ADD_Update_Node_Data_D_x_D_y(u, x);

  else
    {
      cout << endl << " bc_flag must be 0, ERROR in Update_Velo_By_Pressure() at Driver.cc() !!! " << endl ;
      exit(-1);
    } 
}

//************************************************************************************************************************
//function to update phi, with out flow BC, phi at t_n+1, i = Nx + 1 is given

void Update_phi_OutFlow(Node_Data& phi, const valarray<double>& x, const valarray<double>& Rph_np1){

  int i, j, Nx, Ny;

  Nx = phi.get_Nx();
  Ny = phi.get_Ny();

  if((Nx+1)*(Ny+1) != x.size() || Rph_np1.size() != Ny + 1)
    {
      cout << endl << " Dimension Incompatible in  Update_phi_OutFlow(), ERROR at line " << __LINE__ << " of file "
           << __FILE__ << endl;
      exit(-1);
    }

  for(j = 0; j <= Ny; j++)
    {
      for(i = 0; i <= Nx; i++)
	phi[i][j] = x[j*(Nx+1) + i];

      phi[Nx+1][j] = Rph_np1[j];
    }
}


//************************************************************************************************************************

void Update_c(Node_Data& c, const valarray<double>& x, int bc_flag){
  
  if(bc_flag == 0)    
    Update_Node_Data_P_x_D_N_y(c, x);
  else if(bc_flag == 1) 
    Update_Node_Data_N_x_D_N_y(c, x);
}

//************************************************************************************************************************

void  Update_c_OutFlow(Node_Data& c, const valarray<double>& x, const valarray<double>& Lc, const valarray<double>& Rc){

  int i, j, Nx, Ny;

  Nx = c.get_Nx();
  Ny = c.get_Ny();

  if(Nx*(Ny+1) != x.size() || Lc.size() != Ny + 1 || Rc.size() != Ny + 1)
    {
      cout << endl << " Dimension Incompatible in  Update_c_OutFlow(), ERROR at line " << __LINE__ << " of file "
           << __FILE__ << endl;
      exit(-1);
    }

  for(j = 0; j <= Ny; j++)
    {
      for(i = 1; i <= Nx; i++)
	c[i][j] = x[j*Nx + (i-1)];

      c[0][j] = Lc[j];
      c[Nx+1][j] = Rc[j];
    }
}

//************************************************************************************************************************

//function to update the variable diffusion coefficient of various species depending on the calcite volume fraction phi_c
void Update_Diffusion_Node(Node_Data& D_Ur_n, Node_Data& D_NH3_n, Node_Data& D_NH4_n, Node_Data& D_H2CO3_n,
                           Node_Data& D_HCO3_n, Node_Data& D_CO3_n, Node_Data& D_Ca_n, Node_Data& D_Cl_n, 
                           Node_Data& D_H_n, Node_Data& D_OH_n,  
                           Node_Data& D_c_n, const Node_Data& phic, const Control_Parameter& CP){

  double delta; //   phic_tol, 

//   phic_tol = CP.D_phic_tol;
  
  delta = CP.D_delta;

  Eff_Diff_Node_sis(D_Ur_n, phic, D_ur, delta); //, phic_tol);
  Eff_Diff_Node_sis(D_NH3_n, phic, D_nh3, delta); //, phic_tol, delta);
  Eff_Diff_Node_sis(D_NH4_n, phic, D_nh4, delta); //, phic_tol, delta);
  Eff_Diff_Node_sis(D_H2CO3_n, phic, D_h2co3, delta); //, phic_tol, delta);
  Eff_Diff_Node_sis(D_HCO3_n, phic, D_hco3, delta); //, phic_tol, delta);
  Eff_Diff_Node_sis(D_CO3_n, phic, D_co3, delta); //, phic_tol, delta);
  Eff_Diff_Node_sis(D_Ca_n, phic, D_ca, delta); //, phic_tol, delta);
  Eff_Diff_Node_sis(D_Cl_n, phic, D_cl, delta); //, phic_tol, delta);
  Eff_Diff_Node_sis(D_H_n, phic, D_h, delta); //, phic_tol, delta);
  Eff_Diff_Node_sis(D_OH_n, phic, D_oh, delta); //, phic_tol, delta);
  Eff_Diff_Node_sis(D_c_n, phic, Ds, delta); //, phic_tol, delta);
}
//************************************************************************************************************************
//function to average two Node_Data, the result is put into the second argument, and the first argument
//is constant
void Average_Node_Data(const Node_Data& v, Node_Data& u){

  int i, Nx;

  if(v.get_Nx() != u.get_Nx() || v.get_Ny() != u.get_Ny())
    {
      cout << endl << "Dimension Incompatible in Average_Node_Data() in Driver.cc !!! " << endl;
      exit(-1);
    }

  Nx = v.get_Nx();

  for(i = 0; i <= Nx + 1; i++)
    u[i] = .5*(u[i] + v[i]);
}

//************************************************************************************************************************

//function to calculate the time step dt at a given node point
//given velocity u, v, sound of speed a, CFL number, dx, dy 
double Explicit_dt(double u, double v, double a, double dx, double dy, double cfl){

  double dt = fabs(u)/dx + fabs(v)/dy + a*sqrt(1.0/dx/dx + 1.0/dy/dy);

  dt = cfl/dt;

  return dt;
}

//************************************************************************************************************************

//function to calculate the time step for the time integration, 
//which is the smallest dt from all the node points, note here we treat speed of sound as a constant
double Calculate_dt(const Node_Data& u, const Node_Data& v, double a, double dx, double dy, double cfl){

  double dt, dt_min;

  int i, j, Nx, Ny;

  dt_min = 100;

  Nx = u.get_Nx();
  Ny = u.get_Ny();

  for(j = 0; j <= Ny; j++)
    for(i = 0; i <= Nx + 1; i++)
      {
        dt = Explicit_dt(u[i][j], v[i][j], a, dx, dy, cfl);
        if(dt_min > dt)
	  dt_min = dt;
      }

  return dt_min;
}

//************************************************************************************************************************

double Cal_Ave_Viscosity(const Node_Data& phi_c, const Node_Data& phi_b, Node_Data& mu, const Control_Parameter CP){

  int i, j, Nx, Ny;

  double result, tmp;

  Nx = phi_c.get_Nx();
  Ny = phi_b.get_Ny();

  result = 0.0;

  //result is the biggest entry of phi
 
  if(CP.Ave_Vis_flag == 0) //simple average of the viscosity
    for(i = 0; i <= Nx + 1; i++)
      for(j = 0; j <= Ny; j++)
	{
	  tmp = phi_c[i][j]*eta_cac + phi_b[i][j]*eta_bio + (1.0 - phi_c[i][j] - phi_b[i][j])*eta_sol;

	  mu[i][j] = tmp;

	  if(result < tmp)
            result = tmp;
	}
  else   if(CP.Ave_Vis_flag == 1) //harmonic average of the viscosity
    for(i = 0; i <= Nx + 1; i++)
      for(j = 0; j <= Ny; j++)
	{
	  tmp = 1.0/(phi_c[i][j]/eta_cac + phi_b[i][j]/eta_bio + (1.0 - phi_c[i][j] - phi_b[i][j])/eta_sol);

	  mu[i][j] = tmp;

	  if(result < tmp)
	    result = tmp;
	}
  else
    {
      cout << endl << " CP. Ave_Vis_flag must be 0 or 1, error in line " << __LINE__ << " of file " << __FILE__ << endl;
      exit(-1);
    }

  return result;
}


//************************************************************************************************************************
void Balance_phi(Node_Data& phi_c, Node_Data& phi_b, Node_Data& phi_s){

    int i, j, Nx, Ny;

    Nx = phi_c.get_Nx();
    Ny = phi_c.get_Ny();

    double tmp;

    for(i = 0; i <= Nx+1; i++)
        for(j = 0; j <= Ny; j++)
            {
                if(phi_c[i][j] < 0.0)
                    phi_c[i][j] = 0.0;
                if(phi_b[i][j] < 0.0)
                    phi_b[i][j] = 0.0;

                tmp = phi_c[i][j] + phi_b[i][j];
                if(tmp > 1.0)
                    {
                        phi_c[i][j] /= tmp;
                        phi_b[i][j] /= tmp;
                    }

                phi_s[i][j] = 1.0 - phi_c[i][j] - phi_b[i][j];
            }
}


//************************************************************************************************************************
//function to calculate the value of a quantity q at the right boundary by the drift BC, i.e., convected by a knonwn drift velocity U
//Q_np1_Nxp1: value at t_n+1, at i = Nx + 1, 
//Q_n_Nxp1, Q_n_Nx are values at t_n, at i = Nx+1, Nx respectively
//D is the diffusion constant used at two end points of the boundary to remove the artificial singularity 
void Cal_Drift_R_BC(valarray<double>& Q_np1_Nxp1, const valarray<double>& U, const valarray<double>& Q_n_Nxp1, const valarray<double>& Q_n_Nx,
                    double dx, double dy, double dt, int Ny, double D){

  if(Q_np1_Nxp1.size() != Ny + 1 || U.size() != Ny + 1 || Q_n_Nxp1.size() != Ny + 1 || Q_n_Nx.size() != Ny + 1)
    {
      cout << endl << "Dimension incompatible in function Cal_Drift_R_BC(), at line " << __LINE__ << " of file "
           <<  __FILE__ << ", ERROR !!" << endl;
      exit(-1);
    }

  int j;

  double beta;

  for(j = 1; j <= Ny-1; j++)
    {
      beta = U[j]*dt/dx;

      Q_np1_Nxp1[j] = (1.0 - beta)*Q_n_Nxp1[j] + beta*Q_n_Nx[j];
    }  

  //treat the ends points at j = 0, Ny separately, with implicit diffusion in y-direction added
  beta =  2.0*D*dt/dy/dy;

  Q_np1_Nxp1[0] = (1.0/(1.0 + beta))*( (1.0 - U[0]*dt/dx)*Q_n_Nxp1[0] +  (U[0]*dt/dx)*Q_n_Nx[0] + beta*Q_np1_Nxp1[1]);
  
  Q_np1_Nxp1[Ny]= (1.0/(1.0 + beta))*( (1.0 - U[Ny]*dt/dx)*Q_n_Nxp1[Ny] +  (U[Ny]*dt/dx)*Q_n_Nx[Ny] + beta*Q_np1_Nxp1[Ny-1]);
}

//************************************************************************************************************************
//function to calculate Rph_cor, where phi_np1_Nxp2 = - phi_np1_nx + Rph_cor
//Rph_np1: phi value at t_n+1, i = Nx + 1, Lap_ph_np1_Nxp1: Laplacian of phi at t_n+1, i = Nx + 1
void Cal_Rph_cor(valarray<double>& Rph_cor, const valarray<double>& Rph_np1, const valarray<double>& Lap_ph_np1_Nxp1, double dx, double dy, int Ny){

  if(Rph_cor.size() != Ny + 1 || Rph_np1.size() != Ny + 1 || Lap_ph_np1_Nxp1.size() != Ny + 1)
    {
      cout << endl << "Dimension incompatible in function Cal_Rph_cor(), at line " << __LINE__ << " of file " 
           << __FILE__ << ", ERROR !!" << endl;
      exit(-1);
    }

  int j;

  double beta;

  beta = (dx*dx)/(dy*dy);

  for(j = 1; j <= Ny - 1; j++)
    {
      Rph_cor[j] = (2 + 2*beta)*Rph_np1[j] - beta*(Rph_np1[j+1] + Rph_np1[j-1]) + (dx*dx)*Lap_ph_np1_Nxp1[j];
    }

  j = 0;

  Rph_cor[j] = (2 + 2*beta)*Rph_np1[j] - beta*(Rph_np1[j+1] + Rph_np1[j+1]) + (dx*dx)*Lap_ph_np1_Nxp1[j];
   
  j = Ny;
    
  Rph_cor[j] = (2 + 2*beta)*Rph_np1[j] - beta*(Rph_np1[j-1] + Rph_np1[j-1]) + (dx*dx)*Lap_ph_np1_Nxp1[j];
}

//************************************************************************************************************************
//function to calculate Laplacian of phi at i = Nx, no flux BC in y-direction
void Cal_Lap_ph_Nx(valarray<double>& Lap_ph_Nx, const Node_Data& phi, double dx, double dy){

  int i, j, Nx, Ny;

  Nx = phi.get_Nx();
  Ny = phi.get_Ny();

  if(Lap_ph_Nx.size() != Ny + 1)
    {
      cout << endl << "Dimension incompatible in function Cal_Lap_ph_Nx(), at line " << __LINE__ << " of file " 
           << __FILE__ << ", ERROR !!" << endl;
      exit(-1);
    }     

  double Idx2, Idy2;

  i = Nx;

  Idx2 = 1.0/(dx*dx);
  Idy2 = 1.0/(dy*dy);

  for(j = 1; j <= Ny - 1; j++)
    {
      Lap_ph_Nx[j] = Idx2*(phi[i+1][j] - 2*phi[i][j] + phi[i-1][j]) + Idy2*(phi[i][j+1] - 2*phi[i][j] + phi[i][j-1]);
    }

  j = 0;

  Lap_ph_Nx[j] = Idx2*(phi[i+1][j] - 2*phi[i][j] + phi[i-1][j]) + Idy2*(phi[i][j+1] - 2*phi[i][j] + phi[i][j+1]);

  j = Ny;

  Lap_ph_Nx[j] = Idx2*(phi[i+1][j] - 2*phi[i][j] + phi[i-1][j]) + Idy2*(phi[i][j-1] - 2*phi[i][j] + phi[i][j-1]);
}


//************************************************************************************************************************
//function to calculate the drift velocity (x-component), may depend on time and space
void Cal_Drift_U(valarray<double>& U, double t_curr, double dy, int Ny, const Control_Parameter& CP){

  if(U.size() != Ny + 1)
    {
      cout << endl << "Dimension incompatible in function Cal_Drift_U(), at line " << __LINE__ << " of file " 
           << __FILE__ << ", ERROR !!" << endl;
      exit(-1);
    } 

  int j;

  double y, U_max, Ly;

  Ly = dy*Ny;

//   if(t_curr <= 1.0)
//     U_max = CP.Shear_v*t_curr;
//   else
    U_max = CP.Shear_v;

  for(j = 0; j <= Ny; j++)
    {
      y = j*dy; 

      U[j] = U_max*(4.0/Ly/Ly)*y*(Ly - y); 
    }
}

//************************************************************************************************************************
//function to calculate the value of nutrient c at the left boundary, may depend on time and space
void Cal_Lc(valarray<double>& Lc, double t_curr, double dy, int Ny, double In_value){

  if(Lc.size() != Ny + 1)
    {
      cout << endl << "Dimension incompatible in function Cal_Lc(), at line " << __LINE__ << " of file " 
           << __FILE__ << ", ERROR !!" << endl;
      exit(-1);
    } 

  int j;

  double y;

  for(j = 0; j <= Ny; j++)
    {
      y = j*dy;
      Lc[j] = In_value; // *4*y*(1.0 - y);
    }

}
//************************************************************************************************************************

//function to calculate the effluent of substance c (integrate it at x = L from y = 0 to y = 1) by trapezoidal rule
double Cal_Effluent(const Node_Data& c, double dy){

  int j, Nx, Ny;

  double result = 0.0;

  Nx = c.get_Nx();
  Ny = c.get_Ny();

  for(j = 1; j <= Ny - 1; j++)
    result += c[Nx+1][j];

  result += .5*(c[Nx+1][0] + c[Nx+1][Ny]);

  return result*dy;
}

//************************************************************************************************************************

//function to calculate the influent of substance c (integrate it at x = 0 from y = 0 to y = 1)  by trapezoidal rule
double Cal_Influent(const Node_Data& c, double dy){

  int j, Nx, Ny;

  double result = 0.0;

  Nx = c.get_Nx();
  Ny = c.get_Ny();

  for(j = 1; j <= Ny - 1; j++)
    result += c[0][j];

  result += .5*(c[0][0] + c[0][Ny]);

  return result*dy;
}

//************************************************************************************************************************
//the main driver function, integrate over time
void Time_Evolve(char* filename){

  int  Init_flag,  BC_flag, Grid_flag, GMRES_flag, GMRES_max_k, GMRES_num_restart, GMRES_num_k, F_Stride;

  size_t Nx, Ny, k;

  double dx, dy, T_final, CFL, GMRES_tol, GMRES_res, GMRES_Scale, Tfac, alpha, Lam_Helm_Velo, EP_source;

  double Efflu_curr[15];

  time_t start, end, lin_s, lin_e;
  double tdiff;

  char control_file[80];
  char global_file[80];
  char Inter_file[80];   //Saved Inter file contains all variables
  char Ini_file[13][80]; //Thirteen initial files, each for one variable
  char datacomment[250];
  char debfile[80];

  ifstream ifile(filename,ios::in);

  if(!ifile)
    {
        cout << endl << "Can not open input file " << filename << ", ERROR in Time_Eolve() at line " << __LINE__ << " of file " << __FILE__  << endl;
      exit(-1);
    } 

  //Read in the name of files containing the data

  ifile >> datacomment >> control_file
	>> datacomment >> global_file
	>> datacomment >> Inter_file
  	>> datacomment >> Ini_file[0]
  	>> datacomment >> Ini_file[1]
  	>> datacomment >> Ini_file[2]
  	>> datacomment >> Ini_file[3]
  	>> datacomment >> Ini_file[4]
  	>> datacomment >> Ini_file[5]
  	>> datacomment >> Ini_file[6]
 	>> datacomment >> Ini_file[7]
  	>> datacomment >> Ini_file[8]
  	>> datacomment >> Ini_file[9]
  	>> datacomment >> Ini_file[10]
  	>> datacomment >> Ini_file[11]
  	>> datacomment >> Ini_file[12];

         

  ifile.close();

  //read in the control parameters
  Control_Parameter CP;
  Global_Parameter  GP;

  Read_Control_Parameters(control_file, CP);

  dx = CP.dx;
  dy = CP.dy;
  T_final = CP.T_final;
  CFL = CP.CFL;
  GMRES_tol = CP.GMRES_tol;
  GMRES_max_k = CP.GMRES_max_k;  
  Nx = CP.Nx;
  k = CP.k;
  Ny = (1 << (k+1));  //Ny = m = 2^(k+1)
  Init_flag = CP.Init_flag;
  BC_flag = CP.BC_flag;
  Grid_flag = CP.Grid_flag;
  F_Stride = CP.F_Stride; //number of files skipped between the saved ones
  Tfac = CP.tfac;
  alpha = CP.alph;
 
  //read in the global parameters
  Read_Global_Parameters(global_file, GP, CP);

  //*************************DEBUG*****************************

  sprintf(debfile,"%sCon_para",SaveDir);
  CP.display(debfile);
  sprintf(debfile,"%sGlo_para",SaveDir);
  GP.display(debfile);
  Show_NonDimension_Para();
  sprintf(debfile,"%sNondim_Para",SaveDir);
  Save_NonDimension_Para(debfile);

  //*************************DEBUG*****************************

  //define the local  variables

  int i, j, n, Inter_file_NO;
  
  double dt_n, T_current;
  n = 0;
  Inter_file_NO = 0;
  T_current = 0.0;

  Node_Data tmp_Node(Nx, Ny);
  
  Node_Data un(Nx, Ny);
  Node_Data vn(Nx, Ny);

  valarray<double> Drift_U(0.0, Ny + 1);
  valarray<double> Rphb_np1_U(0.0, Ny + 1);
  valarray<double> Rphb_n_U(0.0, Ny + 1);  
  valarray<double> Rphc_np1_U(0.0, Ny + 1);
  valarray<double> Rphc_n_U(0.0, Ny + 1);  
  valarray<double> Rmu_np1(0.0, Ny + 1);
  valarray<double> Rmu_n(0.0, Ny + 1);  

  Node_Data phic_n(Nx, Ny);
  Node_Data phib_n(Nx, Ny);
  Node_Data new_phic(Nx, Ny);
  Node_Data new_phib(Nx, Ny);  

  valarray<double> Rph_np1(0.0, Ny + 1);
  valarray<double> Rph_cor(0.0, Ny + 1);
  valarray<double> Lap_phc_n_Nx(0.0, Ny + 1);
  valarray<double> Lap_phc_n_Nxp1(0.0, Ny + 1);
  valarray<double> Lap_phc_np1_Nxp1(0.0, Ny + 1);
  valarray<double> Lap_phb_n_Nx(0.0, Ny + 1);
  valarray<double> Lap_phb_n_Nxp1(0.0, Ny + 1);
  valarray<double> Lap_phb_np1_Nxp1(0.0, Ny + 1);

  Node_Data phis_np1(Nx, Ny);  
  Node_Data phis_n(Nx, Ny);

  Node_Data Ur_n(Nx, Ny);
  Node_Data NH3_n(Nx, Ny);
  Node_Data NH4_n(Nx, Ny);
  Node_Data H2CO3_n(Nx, Ny);
  Node_Data HCO3_n(Nx, Ny);
  Node_Data CO3_n(Nx, Ny);
  Node_Data Ca_n(Nx, Ny);
  Node_Data Cl_n(Nx, Ny);
  Node_Data H_n(Nx, Ny);
  Node_Data OH_n(Nx, Ny);
  Node_Data c_n(Nx, Ny);

  Node_Data new_Ur(Nx, Ny);  

  Node_Data D_Ur_n(Nx, Ny);
  Node_Data D_NH3_n(Nx, Ny);
  Node_Data D_NH4_n(Nx, Ny);
  Node_Data D_H2CO3_n(Nx, Ny);
  Node_Data D_HCO3_n(Nx, Ny);
  Node_Data D_CO3_n(Nx, Ny);
  Node_Data D_Ca_n(Nx, Ny);
  Node_Data D_Cl_n(Nx, Ny);
  Node_Data D_H_n(Nx, Ny);
  Node_Data D_OH_n(Nx, Ny);
  Node_Data D_c_n(Nx, Ny);

  Node_Data D_Ur_np1(Nx, Ny);
  Node_Data D_NH3_np1(Nx, Ny);
  Node_Data D_NH4_np1(Nx, Ny);
  Node_Data D_H2CO3_np1(Nx, Ny);
  Node_Data D_HCO3_np1(Nx, Ny);
  Node_Data D_CO3_np1(Nx, Ny);
  Node_Data D_Ca_np1(Nx, Ny);
  Node_Data D_Cl_np1(Nx, Ny);
  Node_Data D_H_np1(Nx, Ny);
  Node_Data D_OH_np1(Nx, Ny);
  Node_Data D_c_np1(Nx, Ny);

  //assign constant value to the diffusion coefficients, change it only if CP.D_flag = 1

  D_Ur_n = D_ur;
  D_NH3_n = D_nh3;
  D_NH4_n = D_nh4;
  D_H2CO3_n = D_h2co3;
  D_HCO3_n = D_hco3;
  D_CO3_n = D_co3;
  D_Ca_n = D_ca;
  D_Cl_n = D_cl;
  D_H_n = D_h;
  D_OH_n = D_oh;
  D_c_n = Ds;

  D_Ur_np1 = D_ur;
  D_NH3_np1 = D_nh3;
  D_NH4_np1 = D_nh4;
  D_H2CO3_np1 = D_h2co3;
  D_HCO3_np1 = D_hco3;
  D_CO3_np1 = D_co3;
  D_Ca_np1 = D_ca;
  D_Cl_np1 = D_cl;
  D_H_np1 = D_h;
  D_OH_np1 = D_oh;
  D_c_np1 = Ds;


  Node_Data EP(Nx, Ny);
  Node_Data NT(Nx, Ny);
  Node_Data CT(Nx, Ny);
  Node_Data CHA_SUM(Nx, Ny);

  //value at i = Nx + 2, out-flow BC
  valarray<double> Lc(0.0, Ny+1); 
//   valarray<double> Rc(0.0, Ny+1);
 
  valarray< valarray<double> > Rc_old;
  valarray< valarray<double> > Rc_new;

  Rc_old.resize(12);
  Rc_new.resize(12);

  for(i=0; i < Rc_old.size(); i++)
    {
      Rc_old[i].resize(Ny+1);
      Rc_new[i].resize(Ny+1);
    }

  
  Node_Data SS_n(Nx, Ny); //Saturation-State of CaCO3
  Node_Data Ion_Stren(Nx, Ny); //Ionic strength strength

  Node_Data s(Nx, Ny);
  Node_Data p(Nx, Ny);
                                                        
  double s_stg[(Nx+1)*Ny];   //pressure for half-staggered grid, must initialize it as 0

  for(i = 0; i < (Nx+1)*Ny; i++)
    s_stg[i] = 0.0;

  Node_Data mu(Nx, Ny);  //the average viscosity
 
  Node_Data rho_n(Nx, Ny);

  //variables for velocity update
  double Fx[(Nx+2)*(Ny+1)];
  double Fy[(Nx+2)*(Ny+1)];

  //the following variables are same for all BCs
  MultipleRHS d_press(Ny + 1, Nx + 2);
  double d_press_stg[(Nx+1)*Ny];
  valarray<double> p_x(0.0, (Nx+2)*(Ny-1));
  valarray<double> p_y(0.0, (Nx+2)*(Ny-1));
  Node_Data aux_s(Nx, Ny); 

  //variables for phi_c and phi_b update, they are same data type and can share one set of variables 
  PhiMat A_phi(Nx, Ny, BC_flag); //for No-Flux BC in x and y-direction
  valarray<double> F_phi(0.0, A_phi.get_dim());
  valarray<double> d_phi(0.0, A_phi.get_dim());  //vector contain the new solution of phi

  //variables for Ur, Ca, OH, CT, c update, they are same data type and can share one set of variables 
  c_Mat A_tran(Nx, Ny, 1);  //flag == 1, N_x_N_y
  valarray<double> F_tran(0.0, A_tran.get_dim());
  valarray<double> d_tran(0.0, A_tran.get_dim()); //vector contains new solution of Ur, Ca, OH, CT

  //auxiliary variables for solving the Helmholtz and Poisson equation
  double bda[Ny+1], bdb[Ny+1], bdc[Nx+2], bdd[Nx+2], d_press_reg[(Nx+2)*(Ny+1)];
  double v_bda, v_bdb, v_bdc, v_bdd; 

  c_Mat A_EP(Nx, Ny, 0); //matrix for electric potential, no-flux BC in x and y, flag == 1, N_x_N_y
  valarray<double> F_EP(0.0, A_EP.get_dim());
  valarray<double> d_EP(0.0, A_EP.get_dim());

  int ierror;

  double pertrb;

  //**********************************DEBUG*********************************************************

  ofstream debstat;
  sprintf(debfile,"%sdebstat.txt",SaveDir);
  debstat.open(debfile, ios::out);

  ofstream out_EPsource;
  sprintf(debfile,"%sEP_Source",SaveDir);
  out_EPsource.open(debfile, ios::out);

  //output effluent of various substance with time
  ofstream of_Effluent;
  sprintf(debfile,"%sEffluent",SaveDir);
  of_Effluent.open(debfile, ios::out);

   //************************************************************************************************


  //set the initial condition, set it by hand if Init_flag >= 0
  if(Init_flag >= 0)
    {
      Set_Initial_Condition(Ini_file, un, vn, phic_n, phib_n, Ur_n, NH3_n, NH4_n, H2CO3_n, HCO3_n, CO3_n, Ca_n, Cl_n, H_n, OH_n, c_n, dx, dy, CP);

      //set the value at i = Nx+2 the same as i = Nx+1, used by out-flow BC
      Rc_old[1] = Ur_n[Nx+1];
      Rc_old[2] = NH3_n[Nx+1]; 
      Rc_old[3] = H2CO3_n[Nx+1]; 
      Rc_old[4] = NH4_n[Nx+1]; 
      Rc_old[5] = OH_n[Nx+1]; 
      Rc_old[6] = HCO3_n[Nx+1]; 
      Rc_old[7] = CO3_n[Nx+1]; 
      Rc_old[8] = Ca_n[Nx+1]; 
      Rc_old[9] = Cl_n[Nx+1]; 
      Rc_old[10] = H_n[Nx+1]; 
      Rc_old[11] = c_n[Nx+1]; 



      if(CP.Difu_Init_flag == 1)
	{
	  Diffuse_Initial_Condition(phib_n, CP.N_dif_steps, dx, dy, .1*dx, 1.0, CP.eps_dif); //here dt = dx
	  Diffuse_Initial_Condition(phic_n, CP.N_dif_steps, dx, dy, .1*dx, 1.0, CP.eps_dif); //here dt = dx
	}

      Confine_Node_Data(phib_n, 0, 1);
      Confine_Node_Data(phic_n, 0, 1);

      //also initialize phis, rho_n
      Balance_phi(phic_n, phib_n, phis_n);      
      rho_n = rho_cac*phic_n + rho_bio*phib_n + rho_sol*phis_n;

      //Initialize Diffusion coefficients, do it only if CP.D_flag == 1
      if(CP.D_flag == 1)
	{
	  Update_Diffusion_Node(D_Ur_n, D_NH3_n, D_NH4_n, D_H2CO3_n, D_HCO3_n, D_CO3_n, D_Ca_n, D_Cl_n, D_H_n, D_OH_n, D_c_n, phic_n, CP);
	}

      //calculate Laplacian of phic and phib at i = Nx and Nx + 1
      Cal_Lap_ph_Nx(Lap_phc_n_Nx, phic_n, dx, dy);
      Lap_phc_n_Nxp1 = Lap_phc_n_Nx;
      Cal_Lap_ph_Nx(Lap_phb_n_Nx, phib_n, dx, dy);
      Lap_phb_n_Nxp1 = Lap_phb_n_Nx;

      //get total Nitrogen and Carbon
      NT = NH3_n + NH4_n; 
      CT = H2CO3_n + HCO3_n + CO3_n;
  
      //calculate the total charge at each position at current time, used to conserve charge while solving for [H+] 
      Cal_Charge_Sum(CHA_SUM, H_n, OH_n, Ca_n, Cl_n, NH4_n, HCO3_n, CO3_n);

      Write_All_Node_Data(SaveDir, Inter_file_NO, un, vn, phic_n, phib_n, Ur_n, NH3_n, NH4_n, H2CO3_n, HCO3_n, CO3_n, Ca_n, Cl_n, H_n, OH_n, NT, CT, EP, SS_n,  Ion_Stren, c_n, 
                          T_current, CP);

      eta_ave = Cal_Ave_Viscosity(phic_n, phib_n, mu, CP);

      Rphb_n_U = phib_n[Nx+1];
      Rphc_n_U = phic_n[Nx+1];
      Rmu_n = mu[Nx+1];
 
      Inter_file_NO++;
    }
  else if(Init_flag == -1)
    {
      Read_Inter_Data(Inter_file, Nx, Ny, un, vn, p, s, s_stg, phic_n, phib_n, Ur_n, NH3_n, NH4_n, H2CO3_n, HCO3_n, CO3_n, Ca_n, Cl_n, H_n, OH_n, c_n, 
                      T_current, Inter_file_NO, n);

      //set the value at i = Nx+2 the same as i = Nx+1, used by out-flow BC
      Rc_old[1] = Ur_n[Nx+1];
      Rc_old[2] = NH3_n[Nx+1]; 
      Rc_old[3] = H2CO3_n[Nx+1]; 
      Rc_old[4] = NH4_n[Nx+1]; 
      Rc_old[5] = OH_n[Nx+1]; 
      Rc_old[6] = HCO3_n[Nx+1]; 
      Rc_old[7] = CO3_n[Nx+1]; 
      Rc_old[8] = Ca_n[Nx+1]; 
      Rc_old[9] = Cl_n[Nx+1]; 
      Rc_old[10] = H_n[Nx+1]; 
      Rc_old[11] = c_n[Nx+1]; 

      //initialize other variables
      Balance_phi(phic_n, phib_n, phis_n);      
      rho_n = rho_cac*phic_n + rho_bio*phib_n + rho_sol*phis_n;

      //Initialize Diffusion coefficients, do it only if CP.D_flag == 1
      if(CP.D_flag == 1)
	{
	  Update_Diffusion_Node(D_Ur_n, D_NH3_n, D_NH4_n, D_H2CO3_n, D_HCO3_n, D_CO3_n, D_Ca_n, D_Cl_n, D_H_n, D_OH_n, D_c_n, phic_n, CP);
	}

      //calculate Laplacian of phic and phib at i = Nx and Nx + 1
      Cal_Lap_ph_Nx(Lap_phc_n_Nx, phic_n, dx, dy);
      Lap_phc_n_Nxp1 = Lap_phc_n_Nx;

      Cal_Lap_ph_Nx(Lap_phb_n_Nx, phib_n, dx, dy);
      Lap_phb_n_Nxp1 = Lap_phb_n_Nx;

      //get total Nitrogen and Carbon
      NT = NH3_n + NH4_n; 
      CT = H2CO3_n + HCO3_n + CO3_n;
 

      //calculate the total charge at each position at current time, used to conserve charge while solving for [H+] 
      Cal_Charge_Sum(CHA_SUM, H_n, OH_n, Ca_n, Cl_n, NH4_n, HCO3_n, CO3_n);

      cout << endl << endl << " Resume computation from previously saved data, Inter_file_NO = " << Inter_file_NO 
           << ", current time = " << T_current <<  endl << endl;

      eta_ave = Cal_Ave_Viscosity(phic_n, phib_n, mu, CP);

      Rphb_n_U = phib_n[Nx+1];
      Rphc_n_U = phic_n[Nx+1];
      Rmu_n = mu[Nx+1];

      Inter_file_NO++;
    }

  else
      {
        cout << endl << " Invalid value for Init_flag, it must be bigger than or equal to -1. ERROR in Time_Evolve(), at line " << __LINE__ 
             << " of file " << __FILE__ << endl;
        exit(-1);
      }
   
  
  //now do the time evolution
  while(T_current < T_final) 
    {
      //currently at time step n
      n++;

      dt_n = Tfac*dy;   //set dt = dx

      //Calculate the concentration of species after the fast reactions
      Cal_Fast_Species_Node(H_n, NT, K_nh, K_w, CT, K_acid_1, K_acid_2, Ca_n, Cl_n, CHA_SUM, NH3_n, NH4_n, H2CO3_n, HCO3_n, CO3_n, OH_n, CP.Ion_Act_flag);

//       //Calculate [OH-]
//       OH_n = K_w/H_n;

      //******************************************** DEBUG ********************************************
// #ifdef  DEBUG
//       sprintf(debfile,"%sHF_%d",SaveDir,Inter_file_NO);
//       Write_Node_Data(debfile, H_n);

//       sprintf(debfile,"%sNH3F_%d",SaveDir,Inter_file_NO);
//       Write_Node_Data(debfile, NH3_n);

//       sprintf(debfile,"%sNH4F_%d",SaveDir,Inter_file_NO);
//       Write_Node_Data(debfile, NH4_n);

//       sprintf(debfile,"%sH2CO3F_%d",SaveDir,Inter_file_NO);
//       Write_Node_Data(debfile, H2CO3_n);

//       sprintf(debfile,"%sHCO3F_%d",SaveDir,Inter_file_NO);
//       Write_Node_Data(debfile, HCO3_n);

//       sprintf(debfile,"%sCO3F_%d",SaveDir,Inter_file_NO);
//       Write_Node_Data(debfile, CO3_n);
// #endif

      //******************************************** DEBUG ********************************************

      //now calculate the electric potential, always compute it even if doesn't include its effect in the diffusion equation
      if(CP.EP_flag == 1) //include the effect of electric potential
	{
	  Build_EP_Interior(A_EP, F_EP, NH4_n, HCO3_n, CO3_n, Ca_n, Cl_n, H_n, OH_n, D_NH4_n, D_HCO3_n, D_CO3_n, D_Ca_n, 
			    D_Cl_n, D_H_n, D_OH_n, phis_n, dx, dy, dt_n, a_EPU, EP_source, Inter_file_NO);

	  if((n-1)%F_Stride == 0)
	    {
	      out_EPsource << Inter_file_NO << "  " <<  EP_source << endl;
	    }

	  GMRES_Scale = A_EP.get_max_element();
	  if((n-1)%F_Stride == 0) //out result every F_Stride steps
	    { 
	      debstat<< endl << "GMRES solver for Electric potential at n =  " << n << ", t = " << T_current 
		     <<", max(A) = " << GMRES_Scale << ", max(F) = " << F_EP.max() << ", min(F) = " << F_EP.min() << endl;

	      //********************************** DEBUG *****************************************************************
	      // #ifdef DEBUG
	      // 	  sprintf(debfile,"%sF_EP_%d",SaveDir,Inter_file_NO);
	      // 	  Write_Vector(debfile, F_EP, Nx);
	      // #endif
	      //********************************** DEBUG *****************************************************************
	    }

	  A_EP.Scale();
	  F_EP /= GMRES_Scale;

	  time(&lin_s);
	  GMRES(A_EP, F_EP, d_EP, GMRES_tol, GMRES_max_k, GMRES_flag, GMRES_res, GMRES_num_restart, GMRES_num_k);
	  time(&lin_e);

	  if((n-1)%F_Stride == 0) //out result every F_Stride steps
	    {
	      tdiff = difftime(lin_e, lin_s);
	      debstat<< endl << " Got out of  GMRES solver for Electric potential,"
		     << " GMRES_res = " << GMRES_res << ", GMRES_num_restart = " << GMRES_num_restart  
		     << " GMRES_num_k = " << GMRES_num_k  << ", time used = " << tdiff << " seconds " << endl; 
	      debstat<< endl << "max(d_EP) =  " << d_EP.max() 
		     << ", min(d_EP) = " << d_EP.min() << endl;
	      debstat<< endl << endl << "_________________________________________________________________________________"
		     << endl << endl;
	    }
  
	  Update_EP(EP, d_EP);
	} 

      //calculate the saturation state of CaCO3, this is done after the fast reactions, with or without the ionic strength taken into account

      if(CP.Ion_Act_flag == 1)  //only do this if ionic-strength flag is on
	{
	  Cal_Ionic_Strength_Node(Ion_Stren, H_n, OH_n, Ca_n, Cl_n, NH4_n, HCO3_n, CO3_n);
	}
      SS_n = Saturation_State(Ca_n, CO3_n, Ion_Stren, CP.Ion_Act_flag);


      time(&start);

      //calculate the average viscosity
      eta_ave = Cal_Ave_Viscosity(phic_n, phib_n, mu, CP);

      //u has homogeneous Neumann BC at x = 0 and x = L, un[Nx+2] = un[Nx]
      Cal_Drift_R_BC(Rphc_np1_U, un[Nx], Rphc_n_U, phic_n[Nx+1], dx, dy, dt_n, Ny, 0.0);
      Cal_Drift_R_BC(Rphb_np1_U, un[Nx], Rphb_n_U, phib_n[Nx+1], dx, dy, dt_n, Ny, 0.0);
      Cal_Drift_R_BC(Rmu_np1, un[Nx], Rmu_n, mu[Nx+1], dx, dy, dt_n, Ny, 0.0);

      //Calculate the flux due to constant pressure drop
      Q_flux = Cal_Influent(un, dy);

      //Calculate the pressure gradient
      G_pres = 12.0*(1.0 - L_theta)*Q_flux/L_theta/Re_s - 8.0/L_theta/Re_s;


//       //********************DEBUG PURPOSE*****************************************************
// // #ifdef DEBUG
// //       sprintf(debfile,"%sRU_%d",SaveDir,Inter_file_NO);
// //       Write_Vector(debfile, Drift_U, Ny+1);
// // #endif
//       //********************DEBUG PURPOSE*****************************************************


      //update velocity
      //buid the momentum equation according to BC_flag
     if(Grid_flag == 0) //regular grid
       Build_NSE_RHS_REG_GRID(Fx, Fy, rho_n, mu, Rmu_np1, Rphb_np1_U, Rphc_np1_U, un, vn, s, phib_n, phic_n, dx, dy, dt_n, CP);  
     else if(Grid_flag == 1)  //half-staggered grid
       Build_NSE_RHS_STG_GRID(Fx, Fy, rho_n, mu, Rmu_np1, Rphb_np1_U, Rphc_np1_U, un, vn, s_stg, phib_n, phic_n, dx, dy, dt_n, CP);

//       //******************************************** DEBUG ********************************************
// #ifdef  DEBUG
//       sprintf(debfile,"%sFx_%d",SaveDir,Inter_file_NO);
//       Write_Array(debfile, Fx, (Nx+2)*(Ny+1), Nx+2);

//       sprintf(debfile,"%sFy_%d",SaveDir,Inter_file_NO);
//       Write_Array(debfile, Fy, (Nx+2)*(Ny+1), Nx+2);
// #endif
//       //******************************************** DEBUG ********************************************

     //solve the momentum equation by the FISHPACK Helmholtz solver
     Lam_Helm_Velo = -rho_n[0][0]/eta_ave/dt_n;
      
     time(&lin_s);

     for(int j = 0; j <= Ny; j++)
       {
	 bda[j] = 0.0;
	 bdb[j] = 0.0;
       } 

     //No-flux in both x and y direction
     my_hwscrt(0, dx*(Nx+1), Nx+1, 3, bda, bdb, 0, dy*Ny, Ny, 1, &v_bdc, &v_bdd, Lam_Helm_Velo, Fx, Nx+2, pertrb, ierror);
     if(ierror != 0)
         {
             cout << endl << " Error in Helmholtz solver for x-velocity u at iteration n = " << n << ", ierror = " << ierror << endl;
             exit(-1);
         }

     my_hwscrt(0, dx*(Nx+1), Nx+1, 3, bda, bdb, 0, dy*Ny, Ny, 1, &v_bdc, &v_bdd, Lam_Helm_Velo, Fy, Nx+2, pertrb, ierror);
     if(ierror != 0)
         {
             cout << endl << " Error in Helmholtz solver for y-velocity v at iteration n = " << n << ", ierror = " << ierror << ", pertrb = " << pertrb << endl;
             exit(-1);
         }


      time(&lin_e);

      //update u and v
      Node_Data_Array_Conversion(un, Fx, 2);     
      Node_Data_Array_Conversion(vn, Fy, 2);     



      //******************************************** DEBUG ********************************************
// #ifdef  DEBUG
//       sprintf(debfile,"%sui_%d",SaveDir,Inter_file_NO);
//       Write_Node_Data(debfile, un);

//       sprintf(debfile,"%svi_%d",SaveDir,Inter_file_NO);
//       Write_Node_Data(debfile, vn);
// #endif
      //******************************************** DEBUG ********************************************

      //********************DEBUG PURPOSE*****************************************************
      if((n-1)%F_Stride == 0) //out result every F_Stride steps
	{
	  tdiff = difftime(lin_e, lin_s);

          debstat << endl << "**********************************************************************************" 
                  << endl << " Get out of NSE solver of velocity at n =  " << n << ", t = " << T_current 
                  << endl << "time used = " << tdiff << " seconds, "
	          << " max(du) =  " << un.max() << ", min(du) = " << un.min() 
                  << " max(dv) =  " << vn.max() << ", min(dv) = " << vn.min()
                  << endl;
	}
      //****************************DEBUG********************************************************

       //pressure update, Note here we assume density rho is a constant and use the fast Poisson solver
      //it may need modification for variable density!!!!!!!!!!  here d_press = pressure_phi/rho
      if(Grid_flag == 0)  //regular grid
	Build_Poisson_RHS_REG_GRID((-1.0)*un, (-1.0)*vn, dx, dy, d_press, BC_flag);  //d_press = -div(U)
      else if(Grid_flag == 1) //half-staggered grid
	Build_Poisson_RHS_STG_GRID((-1.0)*un, (-1.0)*vn, dx, dy, d_press_stg, BC_flag); //d_press_stg = -div(U)

      //********************DEBUG PURPOSE*****************************************************
// #ifdef DEBUG
//       sprintf(debfile,"%sp_stg_RHS_%d",SaveDir,Inter_file_NO);
//       Write_Array(debfile, d_press_stg, (Nx+1)*Ny, Nx+1);
//  #endif
      //********************DEBUG PURPOSE*****************************************************


      //update s 
      if(Grid_flag == 0)
	s = s + d_press;
      else if(Grid_flag == 1)
	Add_Array(s_stg, d_press_stg, (Nx+1)*Ny);
      
      //call the fishpack subroutine to solve the poisson equation, put it in a block to avoid name confliction
      {
          //set the homogeneous Neumann BC
          for(int i = 0; i <= Nx+1; i++)
              {
                  bdc[i] = 0.0;
                  bdd[i] = 0.0;
              }
          for(int j = 0; j <= Ny; j++)
              {
                  bda[j] = 0.0;
                  bdb[j] = 0.0;
              }

          //convert the RHS to an array
          if(Grid_flag == 0)   //regular grid
              {
                  MultipleRHS_Array_Conversion(d_press_reg, d_press, 1);

                  //no-flux BC in x-direction
                  my_hwscrt(0.0, dx*(Nx+1), Nx+1, 1, &v_bda, &v_bdb, 0.0, dy*Ny, Ny, 3, bdc, bdd, 0.0, d_press_reg, Nx+2, pertrb, ierror);

                  // convert the RHS to MultipleRHS
                  MultipleRHS_Array_Conversion(d_press_reg, d_press, 2);
              }
          else if(Grid_flag == 1) //half-staggered grid
              {
                  //no-flux BC in x-direction
                  my_hstcrt(0.0, dx*(Nx+1), Nx+1, 1, bda, bdb, 0.0, dy*Ny, Ny, 3, bdc, bdd, 0.0, d_press_stg, Nx+1, pertrb, ierror);
              }

          if(ierror != 0) 
              {
                  cout << endl << "Error in the poisson solver from fishpack for pressure at iteration n = " << n <<", ierror = " << ierror << ", pertrb = " << pertrb << endl;
              }
      }

      if(Grid_flag == 0)
	d_press.To_Gradient(dx, dy, p_x, p_y, BC_flag);
      else if(Grid_flag == 1)
	STG_pressure_gradient(Nx, Ny, dx, dy, d_press_stg, p_x, p_y, BC_flag);

      //********************DEBUG PURPOSE*****************************************************
// #ifdef DEBUG
//       sprintf(debfile,"%sp_stg_%d",SaveDir,Inter_file_NO);
//       Write_Array(debfile, d_press_stg, (Nx+1)*Ny, Nx+1);
//       sprintf(debfile,"%sp_x_%d",SaveDir,Inter_file_NO);
//       Write_Vector(debfile, p_x, Nx+2);
//       sprintf(debfile,"%sp_y_%d",SaveDir,Inter_file_NO);
//       Write_Vector(debfile, p_y, Nx+2);
// #endif
      //********************DEBUG PURPOSE*****************************************************

      Update_Velo_By_Pressure(un, p_x, BC_flag);
      Update_Velo_By_Pressure(vn, p_y, BC_flag);

      Rphb_np1_U = Rphb_n_U;
      Rphc_np1_U = Rphc_n_U;
      Rmu_np1 = Rmu_n;

      //*************************************DEBUG*********************************************
      if((n-1)%F_Stride == 0) //out result every F_Stride steps
	{
	  debstat<< endl << "Pressure update: max(p_x) = " << p_x.max() << ", min(p_x) = " << p_x.min() 
		 << ", max(p_y) = " << p_y.max() << ", min(p_y) = " << p_y.min() 
		 << ", time used = " << tdiff << " seconds " << endl;
    
	}
      //*****************************************************************************************

      //update phi_c, the volume fraction of CaCO3

      //calculate phi at t_n+1, i = Nx + 1
      Cal_Drift_R_BC(Rph_np1, un[Nx+1], phic_n[Nx+1], phic_n[Nx], dx, dy, dt_n, Ny, 0.0);

      //calculate Laplacian phi at t_n+1 at i = Nx + 1
      Cal_Drift_R_BC(Lap_phc_np1_Nxp1, un[Nx+1], Lap_phc_n_Nxp1, Lap_phc_n_Nx, dx, dy, dt_n, Ny, 0.0);

      //calculate Rph_cor, here phi_np1_Nx + 2 = - phi_np1_Nx + Rph_cor
      Cal_Rph_cor(Rph_cor, Rph_np1, Lap_phc_np1_Nxp1, dx, dy, Ny);

      //no-flux BC in x and y-direction
      Build_Phi_c_F_x_N_y_MCH(A_phi, F_phi, Rph_np1, Rph_cor, un, vn, phic_n, phib_n, phis_n, SS_n, dx, dy, dt_n, alpha);

      //get the max-element of A_phi to scale A and F_phi
      GMRES_Scale = A_phi.get_max_element();

      //********************DEBUG PURPOSE*****************************************************
      if((n-1)%F_Stride == 0) //out result every F_Stride steps
          {
	      debstat << endl << "GMRES solver for phi_c at n =  " << n << ",  t = " << T_current 
		      <<", max(A_phi) = " << GMRES_Scale  << ", max(F_phi) = " << F_phi.max() << ", min(F_phi) = " << F_phi.min()  << endl;
          }
      //**************************************************************************************

      //scale A and Fx, Fy
      A_phi.Scale();
      F_phi /= GMRES_Scale;

//       //*************************************DEBUG*********************************************
//       if((n-1)%F_Stride == 0) //out result every F_Stride steps
// 	{
//           for(i = 0; i <= Nx+1; i++)
// 	    for(j = 0; j <= Ny; j++)
// 	      tmp_Node[i][j] = chem_s(phic_n[i][j], phib_n[i][j], phis_n[i][j]);
    
//  	  sprintf(debfile,"%schems_%d",SaveDir,Inter_file_NO);
// 	  Write_Node_Data(debfile, tmp_Node);

//           for(i = 0; i <= Nx+1; i++)
// 	    for(j = 0; j <= Ny; j++)
// 	      tmp_Node[i][j] = chem_b(phic_n[i][j], phib_n[i][j], phis_n[i][j]);
    
//  	  sprintf(debfile,"%schemb_%d",SaveDir,Inter_file_NO);
// 	  Write_Node_Data(debfile, tmp_Node);

//           for(i = 0; i <= Nx+1; i++)
// 	    for(j = 0; j <= Ny; j++)
// 	      tmp_Node[i][j] = chem_c(phic_n[i][j], phib_n[i][j], phis_n[i][j]);
    
//  	  sprintf(debfile,"%schemc_%d",SaveDir,Inter_file_NO);
// 	  Write_Node_Data(debfile, tmp_Node);


// 	  sprintf(debfile,"%sF_phc_%d",SaveDir,Inter_file_NO);
// 	  Write_Vector(debfile, F_phi, Nx+2);
// 	}
//       //*************************************DEBUG*********************************************


      time(&lin_s);

      GMRES(A_phi, F_phi, d_phi, GMRES_tol, GMRES_max_k,  GMRES_flag, GMRES_res, GMRES_num_restart, GMRES_num_k);

      time(&lin_e);

      //********************DEBUG PURPOSE*****************************************************
      if((n-1)%F_Stride == 0) //out result every F_Stride steps
	{
	  tdiff = difftime(lin_e, lin_s);
	  debstat<< endl << " Got out of  GMRES solver for phi_c, GMRES_res = " << GMRES_res << ", GMRES_num_restart = " << GMRES_num_restart  
		 << " GMRES_num_k = " << GMRES_num_k << ", time used = " << tdiff << " seconds " << endl; 
	  debstat<< endl << " max(d_phi_c) =  " << d_phi.max() 
		 << ", min(d_phi_c) = " << d_phi.min() <<  endl;
	  debstat<< endl << endl << "_________________________________________________________________________________"
		 << endl << endl;
	}
      //************************************************************************************** 

      //put new result in a temporary place, update with phi_b together
      Update_phi_OutFlow(new_phic, d_phi, Rph_np1);

      //update related quantities
      Lap_phc_n_Nxp1 = Lap_phc_np1_Nxp1;
    
      Cal_Lap_ph_Nx(Lap_phc_n_Nx, new_phic, dx, dy);

      //update phi_b, the volume fraction of biofilm

      //calculate phi at t_n+1, i = Nx + 1
      Cal_Drift_R_BC(Rph_np1, un[Nx+1], phib_n[Nx+1], phib_n[Nx], dx, dy, dt_n, Ny, 0.0);

      //calculate Laplacian phi at t_n+1 at i = Nx + 1
      Cal_Drift_R_BC(Lap_phb_np1_Nxp1, un[Nx+1], Lap_phb_n_Nxp1, Lap_phb_n_Nx, dx, dy, dt_n, Ny, 0.0);

      //calculate Rph_cor, here phi_np1_Nx + 2 = - phi_np1_Nx + Rph_cor
      Cal_Rph_cor(Rph_cor, Rph_np1, Lap_phb_np1_Nxp1, dx, dy, Ny);


      //no-flux BC in x and y-direction
      Build_Phi_b_F_x_N_y_MCH(A_phi, F_phi, Rph_np1, Rph_cor, un, vn, phib_n, phic_n, phis_n, c_n, dx, dy, dt_n, alpha);

      //get the max-element of A_phi to scale A and F_phi
      GMRES_Scale = A_phi.get_max_element();

      //********************DEBUG PURPOSE*****************************************************
      if((n-1)%F_Stride == 0) //out result every F_Stride steps
          {
	      debstat << endl << "GMRES solver for phi_b at n =  " << n << ",  t = " << T_current 
		      <<", max(A_phi) = " << GMRES_Scale  << ", max(F_phi) = " << F_phi.max() << ", min(F_phi) = " << F_phi.min()  << endl;
          }
      //**************************************************************************************

      //scale A and Fx, Fy
      A_phi.Scale();
      F_phi /= GMRES_Scale;


      time(&lin_s);  

      GMRES(A_phi, F_phi, d_phi, GMRES_tol, GMRES_max_k,  GMRES_flag, GMRES_res, GMRES_num_restart, GMRES_num_k);

      time(&lin_e);

      //********************DEBUG PURPOSE*****************************************************
      if((n-1)%F_Stride == 0) //out result every F_Stride steps
          {
	    tdiff = difftime(lin_e, lin_s);
	    debstat<< endl << " Got out of  GMRES solver for phi_b, GMRES_res = " << GMRES_res << ", GMRES_num_restart = " << GMRES_num_restart  
		   << " GMRES_num_k = " << GMRES_num_k << ", time used = " << tdiff << " seconds " << endl; 
	    debstat<< endl << " max(d_phi_b) =  " << d_phi.max() 
		   << ", min(d_phi_b) = " << d_phi.min() <<  endl;
	    debstat<< endl << endl << "_________________________________________________________________________________"
		   << endl << endl;
          }
      //************************************************************************************** 

      Update_phi_OutFlow(phib_n, d_phi, Rph_np1);

      //update related quantities
      Lap_phb_n_Nxp1 = Lap_phb_np1_Nxp1;
    
      Cal_Lap_ph_Nx(Lap_phb_n_Nx, phib_n, dx, dy);

      //update phi_c now
      phic_n = new_phic;

      //get phis at t_n+1
      Balance_phi(phic_n, phib_n, phis_np1);   

      //update diffusion coefficient, do it only if CP.D_flag == 1
      if(CP.D_flag == 1)
	{
	  Update_Diffusion_Node(D_Ur_np1, D_NH3_np1, D_NH4_np1, D_H2CO3_np1, D_HCO3_np1, D_CO3_np1, D_Ca_np1, D_Cl_np1, 
				D_H_np1, D_OH_np1, D_c_np1, phic_n, CP);
	}


      //******************************************** DEBUG ********************************************
// #ifdef DEBUG
//       sprintf(debfile,"%sDur_%d",SaveDir,Inter_file_NO);
//       Write_Node_Data(debfile, D_Ur_np1);

//       sprintf(debfile,"%sDnh4_%d",SaveDir,Inter_file_NO);
//       Write_Node_Data(debfile, D_NH4_np1);

//       sprintf(debfile,"%sDco3_%d",SaveDir,Inter_file_NO);
//       Write_Node_Data(debfile, D_CO3_np1);

//       sprintf(debfile,"%sDca_%d",SaveDir,Inter_file_NO);
//       Write_Node_Data(debfile, D_Ca_np1);

//       sprintf(debfile,"%sDh_%d",SaveDir,Inter_file_NO);
//       Write_Node_Data(debfile, D_H_np1);
// #endif
      //******************************************** DEBUG ********************************************


      //now update the variables governed by advection-diffusion-reaction equations, 
      //Ur, NH3, NH4, H2CO3, HCO3, CO3, Ca, Cl, c
     //---------------------------------------------------------------------------------------------------------------- 

//       //calculte the value of Ur at the left boundary at t_n+1
//       Cal_Lc(Lc, T_current, dy, Ny, CP.Ini_Ur);

      //calculate Ur at t_n+1, i = Nx + 1
      Cal_Drift_R_BC(Rc_new[1], un[Nx+1], Rc_old[1], Ur_n[Nx+1], dx, dy, dt_n, Ny, D_ur);

      //update [Ur]
      Build_Slow_Process_FluxBC(A_tran, F_tran, 1, D_Ur_n, D_Ur_np1, Z_ur, Ur_n, Rc_new[1], phib_n, Ur_n, SS_n, c_n, 
				phis_np1, phis_n, un, vn, EP, a_EPU, dx, dy, dt_n, CP.Ini_Ur, CP);


      GMRES_Scale = A_tran.get_max_element();
      if((n-1)%F_Stride == 0) //out result every F_Stride steps
	{
	  debstat<< endl << "GMRES solver for Ur at n =  " << n << ", t = " << T_current 
		 <<", max(A) = " << GMRES_Scale << ", max(F) = " << F_tran.max() << ", min(F) = " << F_tran.min() << endl;
	}

      A_tran.Scale();
      F_tran /= GMRES_Scale;

      time(&lin_s);
      GMRES(A_tran, F_tran, d_tran, GMRES_tol, GMRES_max_k, GMRES_flag, GMRES_res, GMRES_num_restart, GMRES_num_k);
      time(&lin_e);

      if((n-1)%F_Stride == 0) //out result every F_Stride steps
          {
	      tdiff = difftime(lin_e, lin_s);
	      debstat<< endl << " Got out of  GMRES solver for Ur,"
		     << " GMRES_res = " << GMRES_res << ", GMRES_num_restart = " << GMRES_num_restart  
		     << " GMRES_num_k = " << GMRES_num_k  << ", time used = " << tdiff << " seconds " << endl; 
              debstat<< endl << "max(d_Ur) =  " << d_tran.max() 
                     << ", min(d_Ur) = " << d_tran.min() << endl;
	      debstat<< endl << endl << "_________________________________________________________________________________"
		     << endl << endl;
          }
	  
      Update_Node_Data_N_x_N_y(new_Ur, d_tran);



     //---------------------------------------------------------------------------------------------------------------- 

//       //calculte the value of NH3 at the left boundary at t_n+1
//       Cal_Lc(Lc, T_current, dy, Ny, CP.Ini_NH3);

      //calculate NH3 at t_n+1, i = Nx + 1
      Cal_Drift_R_BC(Rc_new[2], un[Nx+1], Rc_old[2], NH3_n[Nx+1], dx, dy, dt_n, Ny, D_nh3);

      //update [NH3]
      Build_Slow_Process_FluxBC(A_tran, F_tran, 2, D_NH3_n, D_NH3_np1, Z_nh3, NH3_n, Rc_new[2], phib_n, Ur_n, SS_n, c_n, 
				phis_np1, phis_n, un, vn, EP, a_EPU, dx, dy, dt_n, CP.Ini_NH3, CP);


      GMRES_Scale = A_tran.get_max_element();
      if((n-1)%F_Stride == 0) //out result every F_Stride steps
	{
	  debstat<< endl << "GMRES solver for NH3 at n =  " << n << ", t = " << T_current 
		 <<", max(A) = " << GMRES_Scale << ", max(F) = " << F_tran.max() << ", min(F) = " << F_tran.min() << endl;
	}

      A_tran.Scale();
      F_tran /= GMRES_Scale;

      time(&lin_s);
      GMRES(A_tran, F_tran, d_tran, GMRES_tol, GMRES_max_k, GMRES_flag, GMRES_res, GMRES_num_restart, GMRES_num_k);
      time(&lin_e);

      if((n-1)%F_Stride == 0) //out result every F_Stride steps
          {
	      tdiff = difftime(lin_e, lin_s);
	      debstat<< endl << " Got out of  GMRES solver for NH3,"
		     << " GMRES_res = " << GMRES_res << ", GMRES_num_restart = " << GMRES_num_restart  
		     << " GMRES_num_k = " << GMRES_num_k  << ", time used = " << tdiff << " seconds " << endl; 
              debstat<< endl << "max(d_NH3) =  " << d_tran.max() 
                     << ", min(d_NH3) = " << d_tran.min() << endl;
	      debstat<< endl << endl << "_________________________________________________________________________________"
		     << endl << endl;
          }
	  
      Update_Node_Data_N_x_N_y(NH3_n, d_tran);

      //make sure NH3 is in [0, Ini_NH3]
      Confine_Node_Data(NH3_n, 0, 1.0);

     //---------------------------------------------------------------------------------------------------------------- 

//      //calculte the value of NH4 at the left boundary at t_n+1
//       Cal_Lc(Lc, T_current, dy, Ny, CP.Ini_NH4);

      //calculate NH4 at t_n+1, i = Nx + 1
      Cal_Drift_R_BC(Rc_new[4], un[Nx+1], Rc_old[4], NH4_n[Nx+1], dx, dy, dt_n, Ny, D_nh4);

      //update [NH4]
      Build_Slow_Process_FluxBC(A_tran, F_tran, 4, D_NH4_n, D_NH4_np1, Z_nh4, NH4_n, Rc_new[4], phib_n, Ur_n, SS_n, c_n, 
				phis_np1, phis_n, un, vn, EP, a_EPU, dx, dy, dt_n, CP.Ini_NH4, CP);


      GMRES_Scale = A_tran.get_max_element();
      if((n-1)%F_Stride == 0) //out result every F_Stride steps
	{
	  debstat<< endl << "GMRES solver for NH4 at n =  " << n << ", t = " << T_current 
		 <<", max(A) = " << GMRES_Scale << ", max(F) = " << F_tran.max() << ", min(F) = " << F_tran.min() << endl;
	}

      A_tran.Scale();
      F_tran /= GMRES_Scale;

      time(&lin_s);
      GMRES(A_tran, F_tran, d_tran, GMRES_tol, GMRES_max_k, GMRES_flag, GMRES_res, GMRES_num_restart, GMRES_num_k);
      time(&lin_e);

      if((n-1)%F_Stride == 0) //out result every F_Stride steps
          {
	      tdiff = difftime(lin_e, lin_s);
	      debstat<< endl << " Got out of  GMRES solver for NH4,"
		     << " GMRES_res = " << GMRES_res << ", GMRES_num_restart = " << GMRES_num_restart  
		     << " GMRES_num_k = " << GMRES_num_k  << ", time used = " << tdiff << " seconds " << endl; 
              debstat<< endl << "max(d_NH4) =  " << d_tran.max() 
                     << ", min(d_NH4) = " << d_tran.min() << endl;
	      debstat<< endl << endl << "_________________________________________________________________________________"
		     << endl << endl;
          }
	  
      Update_Node_Data_N_x_N_y(NH4_n, d_tran);

      //make sure NH4 is in [0, Ini_NH4]
      Confine_Node_Data(NH4_n, 0, 1.0);

     //---------------------------------------------------------------------------------------------------------------- 

//      //calculte the value of H2CO3 at the left boundary at t_n+1
//       Cal_Lc(Lc, T_current, dy, Ny, CP.Ini_H2CO3);

      //calculate H2CO3 at t_n+1, i = Nx + 1
      Cal_Drift_R_BC(Rc_new[3], un[Nx+1], Rc_old[3], H2CO3_n[Nx+1], dx, dy, dt_n, Ny, D_h2co3);

      //update [H2CO3]
      Build_Slow_Process_FluxBC(A_tran, F_tran, 3, D_H2CO3_n, D_H2CO3_np1, Z_h2co3, H2CO3_n, Rc_new[3], phib_n, Ur_n, SS_n, c_n, 
				phis_np1, phis_n, un, vn, EP, a_EPU, dx, dy, dt_n, CP.Ini_H2CO3, CP);


      GMRES_Scale = A_tran.get_max_element();
      if((n-1)%F_Stride == 0) //out result every F_Stride steps
	{
	  debstat<< endl << "GMRES solver for H2CO3 at n =  " << n << ", t = " << T_current 
		 <<", max(A) = " << GMRES_Scale << ", max(F) = " << F_tran.max() << ", min(F) = " << F_tran.min() << endl;
	}

      A_tran.Scale();
      F_tran /= GMRES_Scale;

      time(&lin_s);
      GMRES(A_tran, F_tran, d_tran, GMRES_tol, GMRES_max_k, GMRES_flag, GMRES_res, GMRES_num_restart, GMRES_num_k);
      time(&lin_e);

      if((n-1)%F_Stride == 0) //out result every F_Stride steps
          {
	      tdiff = difftime(lin_e, lin_s);
	      debstat<< endl << " Got out of  GMRES solver for H2CO3,"
		     << " GMRES_res = " << GMRES_res << ", GMRES_num_restart = " << GMRES_num_restart  
		     << " GMRES_num_k = " << GMRES_num_k  << ", time used = " << tdiff << " seconds " << endl; 
              debstat<< endl << "max(d_H2CO3) =  " << d_tran.max() 
                     << ", min(d_H2CO3) = " << d_tran.min() << endl;
	      debstat<< endl << endl << "_________________________________________________________________________________"
		     << endl << endl;
          }
	  
      Update_Node_Data_N_x_N_y(H2CO3_n, d_tran);

      //make sure H2CO3 is in [0, Ini_H2CO3]
      Confine_Node_Data(H2CO3_n, 0, 1.0);

     //---------------------------------------------------------------------------------------------------------------- 

//      //calculte the value of HCO3 at the left boundary at t_n+1
//       Cal_Lc(Lc, T_current, dy, Ny, CP.Ini_HCO3);

      //calculate HCO3 at t_n+1, i = Nx + 1
      Cal_Drift_R_BC(Rc_new[6], un[Nx+1], Rc_old[6], HCO3_n[Nx+1], dx, dy, dt_n, Ny, D_hco3);

      //update [HCO3]
      Build_Slow_Process_FluxBC(A_tran, F_tran, 6, D_HCO3_n, D_HCO3_np1, Z_hco3, HCO3_n, Rc_new[6], phib_n, Ur_n, SS_n, c_n, 
				phis_np1, phis_n, un, vn, EP, a_EPU, dx, dy, dt_n, CP.Ini_HCO3, CP);


      GMRES_Scale = A_tran.get_max_element();
      if((n-1)%F_Stride == 0) //out result every F_Stride steps
	{
	  debstat<< endl << "GMRES solver for HCO3 at n =  " << n << ", t = " << T_current 
		 <<", max(A) = " << GMRES_Scale << ", max(F) = " << F_tran.max() << ", min(F) = " << F_tran.min() << endl;
	}

      A_tran.Scale();
      F_tran /= GMRES_Scale;

      time(&lin_s);
      GMRES(A_tran, F_tran, d_tran, GMRES_tol, GMRES_max_k, GMRES_flag, GMRES_res, GMRES_num_restart, GMRES_num_k);
      time(&lin_e);

      if((n-1)%F_Stride == 0) //out result every F_Stride steps
          {
	      tdiff = difftime(lin_e, lin_s);
	      debstat<< endl << " Got out of  GMRES solver for HCO3,"
		     << " GMRES_res = " << GMRES_res << ", GMRES_num_restart = " << GMRES_num_restart  
		     << " GMRES_num_k = " << GMRES_num_k  << ", time used = " << tdiff << " seconds " << endl; 
              debstat<< endl << "max(d_HCO3) =  " << d_tran.max() 
                     << ", min(d_HCO3) = " << d_tran.min() << endl;
	      debstat<< endl << endl << "_________________________________________________________________________________"
		     << endl << endl;
          }
	  
      Update_Node_Data_N_x_N_y(HCO3_n, d_tran);

      //make sure HCO3 is in [0, Ini_HCO3]
      Confine_Node_Data(HCO3_n, 0, 1.0);

     //---------------------------------------------------------------------------------------------------------------- 

//      //calculte the value of CO3 at the left boundary at t_n+1
//       Cal_Lc(Lc, T_current, dy, Ny, CP.Ini_CO3);

      //calculate CO3 at t_n+1, i = Nx + 1
      Cal_Drift_R_BC(Rc_new[7], un[Nx+1], Rc_old[7], CO3_n[Nx+1], dx, dy, dt_n, Ny, D_co3);

      //update [CO3]
      Build_Slow_Process_FluxBC(A_tran, F_tran, 7, D_CO3_n, D_CO3_np1, Z_co3, CO3_n, Rc_new[7], phib_n, Ur_n, SS_n, c_n, 
				phis_np1, phis_n, un, vn, EP, a_EPU, dx, dy, dt_n, CP.Ini_CO3, CP);


      GMRES_Scale = A_tran.get_max_element();
      if((n-1)%F_Stride == 0) //out result every F_Stride steps
	{
	  debstat<< endl << "GMRES solver for CO3 at n =  " << n << ", t = " << T_current 
		 <<", max(A) = " << GMRES_Scale << ", max(F) = " << F_tran.max() << ", min(F) = " << F_tran.min() << endl;
	}

	  //********************************** DEBUG *****************************************************************
#ifdef DEBUG
	  sprintf(debfile,"%sF_CO3_%d",SaveDir,Inter_file_NO);
	  Write_Vector(debfile, F_tran, Nx+2);
#endif
	  //********************************** DEBUG *****************************************************************

      A_tran.Scale();
      F_tran /= GMRES_Scale;

      time(&lin_s);
      GMRES(A_tran, F_tran, d_tran, GMRES_tol, GMRES_max_k, GMRES_flag, GMRES_res, GMRES_num_restart, GMRES_num_k);
      time(&lin_e);

      if((n-1)%F_Stride == 0) //out result every F_Stride steps
          {
	      tdiff = difftime(lin_e, lin_s);
	      debstat<< endl << " Got out of  GMRES solver for CO3,"
		     << " GMRES_res = " << GMRES_res << ", GMRES_num_restart = " << GMRES_num_restart  
		     << " GMRES_num_k = " << GMRES_num_k  << ", time used = " << tdiff << " seconds " << endl; 
              debstat<< endl << "max(d_CO3) =  " << d_tran.max() 
                     << ", min(d_CO3) = " << d_tran.min() << endl;
	      debstat<< endl << endl << "_________________________________________________________________________________"
		     << endl << endl;
          }
	  
      Update_Node_Data_N_x_N_y(CO3_n, d_tran);

      //make sure CO3 is in [0, Ini_CO3]
      Confine_Node_Data(CO3_n, 0, 1.0);

     //---------------------------------------------------------------------------------------------------------------- 

//      //calculte the value of Ca at the left boundary at t_n+1
//       Cal_Lc(Lc, T_current, dy, Ny, CP.Ini_Ca);

      //calculate Ca at t_n+1, i = Nx + 1
      Cal_Drift_R_BC(Rc_new[8], un[Nx+1], Rc_old[8], Ca_n[Nx+1], dx, dy, dt_n, Ny, D_ca);

      //update [Ca]
      Build_Slow_Process_FluxBC(A_tran, F_tran, 8, D_Ca_n, D_Ca_np1, Z_ca, Ca_n, Rc_new[8], phib_n, Ur_n, SS_n, c_n, 
				phis_np1, phis_n, un, vn, EP, a_EPU, dx, dy, dt_n, CP.Ini_Ca, CP);


      GMRES_Scale = A_tran.get_max_element();
      if((n-1)%F_Stride == 0) //out result every F_Stride steps
	{
	  debstat<< endl << "GMRES solver for Ca at n =  " << n << ", t = " << T_current 
		 <<", max(A) = " << GMRES_Scale << ", max(F) = " << F_tran.max() << ", min(F) = " << F_tran.min() << endl;

	  //********************************** DEBUG *****************************************************************
#ifdef DEBUG
	  sprintf(debfile,"%sF_Ca_%d",SaveDir,Inter_file_NO);
	  Write_Vector(debfile, F_tran, Nx+2);
#endif
	  //********************************** DEBUG *****************************************************************
	}

      A_tran.Scale();
      F_tran /= GMRES_Scale;

      time(&lin_s);
      GMRES(A_tran, F_tran, d_tran, GMRES_tol*1e-6, GMRES_max_k, GMRES_flag, GMRES_res, GMRES_num_restart, GMRES_num_k);
      time(&lin_e);

      if((n-1)%F_Stride == 0) //out result every F_Stride steps
          {
	      tdiff = difftime(lin_e, lin_s);
	      debstat<< endl << " Got out of  GMRES solver for Ca,"
		     << " GMRES_res = " << GMRES_res << ", GMRES_num_restart = " << GMRES_num_restart  
		     << " GMRES_num_k = " << GMRES_num_k  << ", time used = " << tdiff << " seconds " << endl; 
              debstat<< endl << "max(d_Ca) =  " << d_tran.max() 
                     << ", min(d_Ca) = " << d_tran.min() << endl;
	      debstat<< endl << endl << "_________________________________________________________________________________"
		     << endl << endl;
          }
	  
      Update_Node_Data_N_x_N_y(Ca_n, d_tran);

      //make sure Ca is in [0, Ini_Ca]
      Confine_Node_Data(Ca_n, 0, 1.0); //CP.Ini_Ca);

     //---------------------------------------------------------------------------------------------------------------- 

//      //calculte the value of Cl at the left boundary at t_n+1
//       Cal_Lc(Lc, T_current, dy, Ny, CP.Ini_Cl);

      //calculate Cl at t_n+1, i = Nx + 1
      Cal_Drift_R_BC(Rc_new[9], un[Nx+1], Rc_old[9], Cl_n[Nx+1], dx, dy, dt_n, Ny, D_cl);

      //update [Cl]
      Build_Slow_Process_FluxBC(A_tran, F_tran, 9, D_Cl_n, D_Cl_np1, Z_cl, Cl_n, Rc_new[9], phib_n, Ur_n, SS_n, c_n, 
				phis_np1, phis_n, un, vn, EP, a_EPU, dx, dy, dt_n, CP.Ini_Cl, CP);


      GMRES_Scale = A_tran.get_max_element();
      if((n-1)%F_Stride == 0) //out result every F_Stride steps
	{
	  debstat<< endl << "GMRES solver for Cl at n =  " << n << ", t = " << T_current 
		 <<", max(A) = " << GMRES_Scale << ", max(F) = " << F_tran.max() << ", min(F) = " << F_tran.min() << endl;

	  //********************************** DEBUG *****************************************************************
// #ifdef DEBUG
// 	  sprintf(debfile,"%sF_Cl_%d",SaveDir,Inter_file_NO);
// 	  Write_Vector(debfile, F_tran, Nx);
// 	  sprintf(debfile,"%sRc_Cl_%d",SaveDir,Inter_file_NO);
// 	  Write_Vector(debfile, Rc, Ny+1);
// 	  sprintf(debfile,"%sLc_Cl_%d",SaveDir,Inter_file_NO);
// 	  Write_Vector(debfile, Lc, Ny+1);
// #endif
	  //********************************** DEBUG *****************************************************************
	}

      A_tran.Scale();
      F_tran /= GMRES_Scale;

      time(&lin_s);
      GMRES(A_tran, F_tran, d_tran, GMRES_tol*1e-6, GMRES_max_k, GMRES_flag, GMRES_res, GMRES_num_restart, GMRES_num_k);
      time(&lin_e);

      if((n-1)%F_Stride == 0) //out result every F_Stride steps
          {
	      tdiff = difftime(lin_e, lin_s);
	      debstat<< endl << " Got out of  GMRES solver for Cl,"
		     << " GMRES_res = " << GMRES_res << ", GMRES_num_restart = " << GMRES_num_restart  
		     << " GMRES_num_k = " << GMRES_num_k  << ", time used = " << tdiff << " seconds " << endl; 
              debstat<< endl << "max(d_Cl) =  " << d_tran.max() 
                     << ", min(d_Cl) = " << d_tran.min() << endl;
	      debstat<< endl << endl << "_________________________________________________________________________________"
		     << endl << endl;
          }
	  
      Update_Node_Data_N_x_N_y(Cl_n, d_tran);

      //make sure Cl is in [0, Ini_Cl]
      Confine_Node_Data(Cl_n, 0, 1.0); //CP.Ini_Cl);

     //---------------------------------------------------------------------------------------------------------------- 

//      //calculte the value of [H+] at the left boundary at t_n+1
//       Cal_Lc(Lc, T_current, dy, Ny, CP.Ini_H);

      //calculate [H+] at t_n+1, i = Nx + 1
      Cal_Drift_R_BC(Rc_new[10], un[Nx+1], Rc_old[10], H_n[Nx+1], dx, dy, dt_n, Ny, D_h);


// 	  //********************************** DEBUG *****************************************************************
// #ifdef DEBUG
// 	  sprintf(debfile,"%sRc_H_%d",SaveDir,Inter_file_NO);
// 	  Write_Vector(debfile, Rc_new[10], Ny+1);
// #endif
// 	  //********************************** DEBUG *****************************************************************

      //update [H+]
      Build_Slow_Process_FluxBC(A_tran, F_tran, 10, D_H_n, D_H_np1, Z_h, H_n, Rc_new[10], phib_n, Ur_n, SS_n, c_n, 
				phis_np1, phis_n, un, vn, EP, a_EPU, dx, dy, dt_n, CP.Ini_H, CP);


      GMRES_Scale = A_tran.get_max_element();
      if((n-1)%F_Stride == 0) //out result every F_Stride steps
	{
	  debstat<< endl << "GMRES solver for H at n =  " << n << ", t = " << T_current 
		 <<", max(A) = " << GMRES_Scale << ", max(F) = " << F_tran.max() << ", min(F) = " << F_tran.min() << endl;

	  //********************************** DEBUG *****************************************************************
// #ifdef DEBUG
// 	  sprintf(debfile,"%sF_H_%d",SaveDir,Inter_file_NO);
// 	  Write_Vector(debfile, F_tran, Nx);
// 	  sprintf(debfile,"%sRc_H_%d",SaveDir,Inter_file_NO);
// 	  Write_Vector(debfile, Rc, Ny+1);
// 	  sprintf(debfile,"%sLc_H_%d",SaveDir,Inter_file_NO);
// 	  Write_Vector(debfile, Lc, Ny+1);
// #endif
	  //********************************** DEBUG *****************************************************************
	}

      A_tran.Scale();
      F_tran /= GMRES_Scale;

      time(&lin_s);
      GMRES(A_tran, F_tran, d_tran, GMRES_tol*1e-6, GMRES_max_k, GMRES_flag, GMRES_res, GMRES_num_restart, GMRES_num_k);
      time(&lin_e);

      if((n-1)%F_Stride == 0) //out result every F_Stride steps
          {
	      tdiff = difftime(lin_e, lin_s);
	      debstat<< endl << " Got out of  GMRES solver for H,"
		     << " GMRES_res = " << GMRES_res << ", GMRES_num_restart = " << GMRES_num_restart  
		     << " GMRES_num_k = " << GMRES_num_k  << ", time used = " << tdiff << " seconds " << endl; 
              debstat<< endl << "max(d_H) =  " << d_tran.max() 
                     << ", min(d_H) = " << d_tran.min() << endl;
	      debstat<< endl << endl << "_________________________________________________________________________________"
		     << endl << endl;
          }
	  
      Update_Node_Data_N_x_N_y(H_n, d_tran);

      //make sure Cl is in [0, Ini_Cl]
      Confine_Node_Data(H_n, 0, 1.0); 

     //---------------------------------------------------------------------------------------------------------------- 

//      //calculte the value of [OH-] at the left boundary at t_n+1
//       Cal_Lc(Lc, T_current, dy, Ny, CP.Ini_OH);

      //calculate [OH-] at t_n+1, i = Nx + 1
      Cal_Drift_R_BC(Rc_new[5], un[Nx+1], Rc_old[5], OH_n[Nx+1], dx, dy, dt_n, Ny, D_oh);

// 	  //********************************** DEBUG *****************************************************************
// #ifdef DEBUG
// 	  sprintf(debfile,"%sRc_OH_%d",SaveDir,Inter_file_NO);
// 	  Write_Vector(debfile, Rc_new[5], Ny+1);
// #endif
// 	  //********************************** DEBUG *****************************************************************

      //update [OH-]
      Build_Slow_Process_FluxBC(A_tran, F_tran, 5, D_OH_n, D_OH_np1, Z_oh, OH_n, Rc_new[5], phib_n, Ur_n, SS_n, c_n, 
				phis_np1, phis_n, un, vn, EP, a_EPU, dx, dy, dt_n, CP.Ini_OH, CP);


      GMRES_Scale = A_tran.get_max_element();
      if((n-1)%F_Stride == 0) //out result every F_Stride steps
	{
	  debstat<< endl << "GMRES solver for OH at n =  " << n << ", t = " << T_current 
		 <<", max(A) = " << GMRES_Scale << ", max(F) = " << F_tran.max() << ", min(F) = " << F_tran.min() << endl;

	  //********************************** DEBUG *****************************************************************
// #ifdef DEBUG
// 	  sprintf(debfile,"%sF_OH_%d",SaveDir,Inter_file_NO);
// 	  Write_Vector(debfile, F_tran, Nx);
// 	  sprintf(debfile,"%sRc_OH_%d",SaveDir,Inter_file_NO);
// 	  Write_Vector(debfile, Rc, Ny+1);
// 	  sprintf(debfile,"%sLc_OH_%d",SaveDir,Inter_file_NO);
// 	  Write_Vector(debfile, Lc, Ny+1);
// #endif
	  //********************************** DEBUG *****************************************************************
	}

      A_tran.Scale();
      F_tran /= GMRES_Scale;

      time(&lin_s);
      GMRES(A_tran, F_tran, d_tran, GMRES_tol*1e-6, GMRES_max_k, GMRES_flag, GMRES_res, GMRES_num_restart, GMRES_num_k);
      time(&lin_e);

      if((n-1)%F_Stride == 0) //out result every F_Stride steps
          {
	      tdiff = difftime(lin_e, lin_s);
	      debstat<< endl << " Got out of  GMRES solver for OH,"
		     << " GMRES_res = " << GMRES_res << ", GMRES_num_restart = " << GMRES_num_restart  
		     << " GMRES_num_k = " << GMRES_num_k  << ", time used = " << tdiff << " seconds " << endl; 
              debstat<< endl << "max(d_OH) =  " << d_tran.max() 
                     << ", min(d_OH) = " << d_tran.min() << endl;
	      debstat<< endl << endl << "_________________________________________________________________________________"
		     << endl << endl;
          }
	  
      Update_Node_Data_N_x_N_y(OH_n, d_tran);

      //make sure Cl is in [0, Ini_Cl]
      Confine_Node_Data(OH_n, 0, 1.0); 

     //---------------------------------------------------------------------------------------------------------------- 

      //Update [Ur], make sure same Ur are used in updating all components above
      Ur_n = new_Ur;
      //make sure Ur is in [0, Ini_Ur]
      Confine_Node_Data(Ur_n, 0, CP.Ini_Ur);

     //---------------------------------------------------------------------------------------------------------------- 

      if(T_current < CP.cutoff_c_time)
	{      
// 	  //calculte the value of c at the left boundary at t_n+1
// 	  Cal_Lc(Lc, T_current, dy, Ny, CP.Ini_c);

	  //calculate c at t_n+1, i = Nx + 1
	  Cal_Drift_R_BC(Rc_new[11], un[Nx+1], Rc_old[11], c_n[Nx+1], dx, dy, dt_n, Ny, Ds);

	  //update [c]
	  Build_Slow_Process_FluxBC(A_tran, F_tran, 11, D_c_n, D_c_np1, Z_c, c_n, Rc_new[11], phib_n, Ur_n, SS_n, c_n, 
				    phis_np1, phis_n, un, vn, EP, a_EPU, dx, dy, dt_n, CP.Ini_c, CP);


	  GMRES_Scale = A_tran.get_max_element();
	  if((n-1)%F_Stride == 0) //out result every F_Stride steps
	    {
	      debstat<< endl << "GMRES solver for c at n =  " << n << ", t = " << T_current 
		     <<", max(A) = " << GMRES_Scale << ", max(F) = " << F_tran.max() << ", min(F) = " << F_tran.min() << endl;
	    }

	  A_tran.Scale();
	  F_tran /= GMRES_Scale;

	  time(&lin_s);
	  GMRES(A_tran, F_tran, d_tran, GMRES_tol, GMRES_max_k, GMRES_flag, GMRES_res, GMRES_num_restart, GMRES_num_k);
	  time(&lin_e);

	  if((n-1)%F_Stride == 0) //out result every F_Stride steps
	    {
	      tdiff = difftime(lin_e, lin_s);
	      debstat<< endl << " Got out of  GMRES solver for c,"
		     << " GMRES_res = " << GMRES_res << ", GMRES_num_restart = " << GMRES_num_restart  
		     << " GMRES_num_k = " << GMRES_num_k  << ", time used = " << tdiff << " seconds " << endl; 
              debstat<< endl << "max(d_c) =  " << d_tran.max() 
                     << ", min(d_c) = " << d_tran.min() << endl;
	      debstat<< endl << endl << "_________________________________________________________________________________"
		     << endl << endl;
	    }
	  
	  Update_Node_Data_N_x_N_y(c_n, d_tran);

	  //make sure c is in [0, Ini_c]
	  Confine_Node_Data(c_n, 0, CP.Ini_c);
	}
      else  //turn off the growth, don't calculate c
	{
	  ;
	}

      //update phis for next timestep
      phis_n = phis_np1;

      //update the out-flow BC
      Rc_old = Rc_new;

      if(CP.D_flag == 1)
	{
	  D_Ur_n = D_Ur_np1;
	  D_NH3_n = D_NH3_np1;
	  D_NH4_n = D_NH4_np1;
	  D_H2CO3_n = D_H2CO3_np1;
	  D_HCO3_n = D_HCO3_np1;
	  D_CO3_n = D_CO3_np1;
	  D_Ca_n = D_Ca_np1;
	  D_Cl_n = D_Cl_np1;
	  D_H_n = D_H_np1;
	  D_OH_n = D_OH_np1;
	  D_c_n = D_c_np1;
	}

      //get total Nitrogen and Carbon and [H+]
      NT = NH3_n + NH4_n; 
      CT = H2CO3_n + HCO3_n + CO3_n;
   
      //calculate the total charge at each position at current time, used to conserve charge while solving for [H+] 
      Cal_Charge_Sum(CHA_SUM, H_n, OH_n, Ca_n, Cl_n, NH4_n, HCO3_n, CO3_n);    
     
      //update the time information
 
      T_current += dt_n;

      time(&end);
      tdiff = difftime(end, start);

      if((n-1)%F_Stride == 0) //out result every F_Stride steps
	{
	  Write_All_Node_Data(SaveDir, Inter_file_NO, un, vn, phic_n, phib_n, Ur_n, NH3_n, NH4_n, H2CO3_n, HCO3_n, CO3_n, Ca_n, Cl_n, H_n, OH_n,
                              NT, CT, EP, SS_n, Ion_Stren, c_n, T_current, CP);


	  sprintf(debfile,"%sP_%d", SaveDir, Inter_file_NO);          
          if(Grid_flag == 0)
            Write_Node_Data(debfile,s);
          else if(Grid_flag == 1)
            Write_Array(debfile, s_stg, (Nx+1)*Ny, Nx+1);

          //calculate the effluent of different substances
          Efflu_curr[0] = Cal_Effluent(phib_n,dy);
          Efflu_curr[1] = Cal_Effluent(Ur_n,dy);
          Efflu_curr[2] = Cal_Effluent(NH3_n,dy);
          Efflu_curr[3] = Cal_Effluent(H2CO3_n,dy);
          Efflu_curr[4] = Cal_Effluent(NH4_n,dy);
          Efflu_curr[5] = Cal_Effluent(OH_n,dy);
          Efflu_curr[6] = Cal_Effluent(HCO3_n,dy);
          Efflu_curr[7] = Cal_Effluent(CO3_n,dy);
          Efflu_curr[8] = Cal_Effluent(Ca_n,dy);
          Efflu_curr[9] = Cal_Effluent(Cl_n,dy);
          Efflu_curr[10] = Cal_Effluent(H_n,dy);
          Efflu_curr[11] = Cal_Effluent(c_n,dy);
          Efflu_curr[12] = Cal_Effluent(phic_n,dy);
          Efflu_curr[13] = Q_flux; // flux due to pressure drop 
          Efflu_curr[14] = G_pres; // pressure gradient in the computation domain

	  of_Effluent << Inter_file_NO << " ";
 
          for(int Ej = 0; Ej <= 13; Ej++)
	    of_Effluent << Efflu_curr[Ej] << " ";
  
          of_Effluent << Efflu_curr[14] << endl;

	  //write the save intermediate data code here
	  if((n-1)%(20*F_Stride) == 0) //out result every F_Stride steps
	    {
	      sprintf(Inter_file, "%sInter_%d.data",SaveDir,Inter_file_NO);
	      Save_Inter_Data(Inter_file, un, vn, p, s, s_stg, phic_n, phib_n, Ur_n, NH3_n, NH4_n, H2CO3_n, HCO3_n, CO3_n,
                              Ca_n, Cl_n, H_n, OH_n, c_n, T_current, Inter_file_NO, n);
	    }

	  //Display current time and number of time steps as well as solution information on screen
          Screen_Display(n, tdiff, dt_n, eta_ave,  T_current, dx, dy, un, vn, phic_n, phib_n, Ur_n, NH3_n, NH4_n, H2CO3_n,
                         HCO3_n, CO3_n, Ca_n, Cl_n, H_n, OH_n, SS_n, c_n);

	  Inter_file_NO++;

 	}

      rho_n = rho_cac*phic_n + rho_bio*phib_n + rho_sol*phis_n;
 
    }
     

  //**********************************DEBUG*********************************************************
  debstat.close();
  out_EPsource.close();
  of_Effluent.close();
  Show_NonDimension_Para();
  //************************************************************************************************
 
     
}
  
