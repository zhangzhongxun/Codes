#include "MatVec.h"
#include "IO.h"

//define a class contains all the parameters used in the 1D diffusion simulation
class DiffPara_1D {

public:

  DiffPara_1D() {} //default constructor does nothing
  ~DiffPara_1D() {} //destructor does nothing

  double dx, D1, D2, Z1, Z2, Tfac, alpha, T_final, Ini_S1, Ini_reg;

  int N, EP_flag, Init_flag, F_Stride;

  char savedir[100];

  void display(char*) const;
};


//************************************************************************************************************************

//function to read in control parameters
void Read_1D_Para(char* filename, DiffPara_1D& CP){

  ifstream ifile(filename, ios::in);

  if( !ifile)
    {
      cout << " \n\t Can not open the control parameters file " << filename << " !! " << endl;
      exit(-1);
    }

  char datacomment[200];

  double dx, dt, D1, D2, Z1, Z2, alpha, T_final, Ini_S1;

  int EP_flag, Init_flag;

  char savedir[100];


  ifile >> datacomment >> CP.dx 
        >> datacomment >> CP.N
 //        >> datacomment >> CP.dt
        >> datacomment >> CP.D1
        >> datacomment >> CP.D2
        >> datacomment >> CP.Z1
        >> datacomment >> CP.Z2
        >> datacomment >> CP.Tfac
        >> datacomment >> CP.alpha
        >> datacomment >> CP.T_final
        >> datacomment >> CP.F_Stride
        >> datacomment >> CP.Ini_S1
        >> datacomment >> CP.Ini_reg
        >> datacomment >> CP.EP_flag
        >> datacomment >> CP.Init_flag
        >> datacomment >> CP.savedir;

  ifile.close();
}

//************************************************************************************************************************

//function to write out control parameters
void DiffPara_1D::display(char* filename) const {

  ofstream ofile(filename,ios::out);
  if(!ofile)
    {
      cout << endl << " Can not open output file " << filename <<", ERROR in  DiffPara_1D::display(), line " << __LINE__ << " of file " << __FILE__ << endl;
      exit(-1);
    }

 ofile << endl << endl << "********************* The Raw Control Parameters of 1D Diffusion ***************************" << endl << endl;

 ofile  << endl << " dx = " <<  dx 
        << endl << " N = " <<  N
 //        << endl << " dt = " <<  dt
        << endl << " D1 = " <<  D1
        << endl << " D2 = " <<  D2
        << endl << " Z1 = " <<  Z1
        << endl << " Z2 = " <<  Z2
        << endl << " Tfac = " << Tfac 
        << endl << " alpha = " <<  alpha
        << endl << " T_final =  " <<  T_final
        << endl << " F_Stride =  " <<  F_Stride
        << endl << " Ini_S1 = " <<  Ini_S1
        << endl << " Ini_reg = " <<  Ini_reg
        << endl << " EP_flag = " <<  EP_flag
        << endl << " Init_flag = " <<  Init_flag
        << endl << " savedir = " <<  savedir;

  ofile.close();
}

//************************************************************************************************************************

//function to set the initial profile of the chemical concentration
  void Set_Ini_1D(valarray<double>& S1, valarray<double>& S2, const DiffPara_1D& CP1D){

    double dx, x, Ini_S1, Ini_reg, mu, sig, Z1, Z2;
    int Init_flag, n, j;

    dx = CP1D.dx;
    Ini_S1 = CP1D.Ini_S1;
    Ini_reg = CP1D.Ini_reg;
    Z1 = CP1D.Z1;
    Z2 = CP1D.Z2;

    mu = 0.5; //mean of the Gaussian distribution
    sig = 0.05; //std of the Gaussian distribution

    Init_flag = CP1D.Init_flag;

    if(S1.size() != S2.size())
      {
	cout << endl << "Size incompatible in Set_Ini_1D, line " << __LINE__ << " of file " << __FILE__ << endl;
	exit(-1);
      }

    n = S1.size();

    for(j = 0; j < n; j++)
      {
	x = dx*j;
        
        if(Init_flag == 0) //step initial value
	  {
	    S1[j] = (x < 0.5) ? Ini_S1+Ini_reg : Ini_reg;
	    S2[j] = -Z1*S1[j]/Z2;
	  }
	else if(Init_flag == 1) // Gaussian initial value
	  {
	    S1[j] = Ini_S1*exp(-(x - mu)*(x - mu)/(2*sig*sig)) + Ini_reg;
	    S2[j] = -Z1*S1[j]/Z2;
	  }
	else
	  {
	    cout << endl << "Initial flag must be either 0 or 1, here it is " << Init_flag << ", error in line " << __LINE__ << " of file " << __FILE__ << endl;
            exit(-1);
	  }
      }
  }


//************************************************************************************************************************



void Build_S(TriDiag& A, valarray<double>& rhs, const valarray<double>& S, const valarray<double>& EP, int EP_flag, double D, double Z,
             double alpha, double dx, double dt){

  if(A.get_n() != rhs.size())
    {
      cout << endl << "Dimension incompatible in Build_S(), line " << __LINE__ << " of file " << __FILE__ << endl;
      exit(-1);
    }  

  int n, i, j;

  double EP_ind;

  EP_ind = (EP_flag == 0) ? 0.0 : 1.0;

  n = A.get_n();

  valarray<double> d(0.0, n);
  valarray<double> e(0.0, n);
  valarray<double> f(0.0, n);

  //Interior nodes first
  for(j = 1; j < n-1; j++) 
    {
      d[j] = 1.0/dt + 2.0*alpha*D/dx/dx - EP_ind*0.5*D*Z*(EP[j+1] - 2.0*EP[j] + EP[j-1])/dx/dx;
      e[j] = -alpha*D/dx/dx + EP_ind*.5*D*Z*(EP[j] - EP[j-1])/dx/dx;
      f[j] = -alpha*D/dx/dx - EP_ind*.5*D*Z*(EP[j+1] - EP[j])/dx/dx;

      rhs[j] = S[j]/dt + (1.0 - alpha)*D*(S[j+1] - 2.0*S[j] + S[j-1])/dx/dx;
    }
  //Left and Right boundary node, no-flux BC
  j = 0;

  d[j] = 1.0/dt + 2.0*alpha*D/dx/dx - EP_ind*0.5*D*Z*(EP[j+1] - 2.0*EP[j] + EP[j+1])/dx/dx;
  e[j] = -alpha*D/dx/dx + EP_ind*.5*D*Z*(EP[j] - EP[j+1])/dx/dx;
  f[j] = -alpha*D/dx/dx - EP_ind*.5*D*Z*(EP[j+1] - EP[j])/dx/dx;
  rhs[j] = S[j]/dt + (1.0 - alpha)*D*(S[j+1] - 2.0*S[j] + S[j+1])/dx/dx;

  j = n-1;

  d[j] = 1.0/dt + 2.0*alpha*D/dx/dx - EP_ind*0.5*D*Z*(EP[j-1] - 2.0*EP[j] + EP[j-1])/dx/dx;
  e[j] = -alpha*D/dx/dx + EP_ind*.5*D*Z*(EP[j] - EP[j-1])/dx/dx;
  f[j] = -alpha*D/dx/dx - EP_ind*.5*D*Z*(EP[j-1] - EP[j])/dx/dx;
  rhs[j] = S[j]/dt + (1.0 - alpha)*D*(S[j-1] - 2.0*S[j] + S[j-1])/dx/dx;

  //Assign values of A
  A.get_d() = d; //diagonal 
  
  for(j = 0; j < n-1; j++)
    {
      A.get_e()[j] = e[j+1];
      A.get_f()[j] = f[j];
    }
  
  A.get_f()[0] += e[0];
  A.get_e()[n-2] += f[n-1]; 
}


//************************************************************************************************************************

void Build_EP(TriDiag& A, valarray<double>& rhs, const valarray<double>& S1, const valarray<double>& S2, double D1, double Z1, double D2, double Z2,
              double dx, double dt){

 if(A.get_n() != rhs.size())
    {
      cout << endl << "Dimension incompatible in Build_S(), line " << __LINE__ << " of file " << __FILE__ << endl;
      exit(-1);
    }  

  int n, i, j;

  n = A.get_n();

  valarray<double> d(0.0, n);
  valarray<double> e(0.0, n);
  valarray<double> f(0.0, n);

  //extend S_i by no-flux BC
  valarray<double> S1_E(0.0, n+2);
  valarray<double> S2_E(0.0, n+2);
 
  for(j = 1; j < n+1; j++)
    {
      S1_E[j] = S1[j-1];
      S2_E[j] = S2[j-1];
    }

  S1_E[0] = S1[1];
  S1_E[n+1] = S1[n-2];

  S2_E[0] = S2[1];
  S2_E[n+1] = S2[n-2];

  //cell center value of S_i by average
  valarray<double> S1_a(0.0, n+1);
  valarray<double> S2_a(0.0, n+1);

  for(j = 0; j < n+1; j++)
    {
      S1_a[j] = .5*(S1_E[j] + S1_E[j+1]);
      S2_a[j] = .5*(S2_E[j] + S2_E[j+1]);
    }

  //Interior nodes first
  for(j = 0; j < n; j++) 
    {
      d[j] = - ( Z1*Z1*D1*(S1_a[j] + S1_a[j+1]) +  Z2*Z2*D2*(S2_a[j] + S2_a[j+1]) );
      e[j] =  Z1*Z1*D1*S1_a[j] + Z2*Z2*D2*S2_a[j];
      f[j] =  Z1*Z1*D1*S1_a[j+1] + Z2*Z2*D2*S2_a[j+1];

      rhs[j] = - ( D1*Z1*(S1_E[j+2] - 2.0*S1_E[j+1] + S1_E[j]) + D2*Z2*(S2_E[j+2] - 2.0*S2_E[j+1] + S2_E[j]) )
	       - ( Z1*S1[j] + Z2*S2[j] )*dx*dx/dt;
    }
 
  //Assign values of A
  A.get_d() = d; //diagonal 
  
  for(j = 0; j < n-1; j++)
    {
      A.get_e()[j] = e[j+1];
      A.get_f()[j] = f[j];
    }
  
  A.get_f()[0] += e[0];
  A.get_e()[n-2] += f[n-1]; 
}


//************************************************************************************************************************

int main(int argc, char* argv[]){

  if(argc < 2)
    {
      cout << endl << "Add the data file after command " << argv[0] << ", Data Missed !!!" << endl;
      exit(-1);
    }

  cout << endl << endl << "**************************** Running    \" " << argv[0] << "  " << argv[1] 
       << " \" ********************************* " << endl << endl;

  time_t start, end;

  double tdiff;

  int hours, minutes, seconds;

  //start clock
  time(&start);

  int N, i, j, n, Inter_file_NO, EP_flag, F_Stride;

  double dx, dt, T_final, T_current, alpha, D1, D2, Z1, Z2, total_S1, total_S2, shift;

  DiffPara_1D CP1D;

  char SaveDir[100];
  char debfile[100];

  //Read parameters
  Read_1D_Para(argv[1], CP1D);

  N = CP1D.N;
  dx = CP1D.dx;
  dt = dx*CP1D.Tfac;
  D1 = CP1D.D1;
  D2 = CP1D.D2;
  Z1 = CP1D.Z1;
  Z2 = CP1D.Z2;
  alpha = CP1D.alpha;
  T_final = CP1D.T_final;
  F_Stride = CP1D.F_Stride;
  strcpy(SaveDir, CP1D.savedir);
  EP_flag = CP1D.EP_flag;

  valarray<double> S1(0.0, N);
  valarray<double> S2(0.0, N);
  valarray<double> EP(0.0, N);

  TriDiag A(N, 0); //Tridiagonal matrix, Dirichelet BC for S (due to the time dependence)
  TriDiag AE(N, 1); //Tridiagonal matrix, no-flux BC for the Electrical potential
  valarray<double> rhs(0.0, N);
  MultipleRHS b(1, N);
  
  n = 0;
  Inter_file_NO = 0;
  T_current = 0.0;

  //set up initial data
  Set_Ini_1D(S1, S2, CP1D);

  //save the total value of S1 and S2 (conserved due to no-flux BC)
  total_S1 = S1.sum();
  total_S2 = S2.sum();

  //save initial data
  sprintf(debfile,"%sS1_%d", SaveDir, Inter_file_NO);
  Write_Vector(debfile, S1, N);
  sprintf(debfile,"%sS2_%d", SaveDir, Inter_file_NO);
  Write_Vector(debfile, S2, N);
  Inter_file_NO++;

  //now do the time evolution
  while(T_current < T_final) 
    {
      //currently at time step n
      n++;

      //calculate the electrical potential if EP_flag == 1
      if(EP_flag == 1)
	{
	  Build_EP(AE, rhs, S1, S2, D1, Z1, D2, Z2, dx, dt);

          //assign valarray type variable rhs to MultipleRHS type variable b
          b[0] = rhs;
          //solve the tridiagonal system
          AE.ChaseSolve(b);
          //save the solution 
          EP = b[0];
	}

      //update S1 and S2
      Build_S(A, rhs, S1, EP, EP_flag, D1, Z1, alpha, dx, dt);

      //assign valarray type variable rhs to MultipleRHS type variable b
      b[0] = rhs;
      //solve the tridiagonal system
      A.ChaseSolve(b);
          //save the solution, shift the value so total amount is conserved
      S1 = b[0];
//       shift = (total_S1 - S1.sum())/N;
//       S1 += shift;

      Build_S(A, rhs, S2, EP, EP_flag, D2, Z2, alpha, dx, dt);
      //assign valarray type variable rhs to MultipleRHS type variable b
      b[0] = rhs;
      //solve the tridiagonal system
      A.ChaseSolve(b);
      //save the solution, shift the value so total amount is conserved
      S2 = b[0];
//       shift = (total_S2 - S2.sum())/N;
//       S2 += shift;

      T_current += dt;

      if((n-1)%F_Stride == 0) //out result every F_Stride steps
	{
	  if(EP_flag == 1)
	    {
	      sprintf(debfile,"%sEP_%d", SaveDir, Inter_file_NO);
	      Write_Vector(debfile, EP, N);
	    }	 
	  sprintf(debfile,"%sS1_%d", SaveDir, Inter_file_NO);
	  Write_Vector(debfile, S1, N);
	  sprintf(debfile,"%sS2_%d", SaveDir, Inter_file_NO);
	  Write_Vector(debfile, S2, N);
	  Inter_file_NO++;
	}
    }

  time(&end);

  tdiff = difftime(end, start);
  seconds = (int)(tdiff);

  hours = seconds/3600;
  seconds %= 3600;
  minutes = seconds/60;
  seconds %= 60;
  
  cout << "\n\t\t\t--------------------------------\n";
  cout << "\t\t\t TOTAL TIME ELAPSED: ";
  cout.width(2); cout.fill('0'); cout << hours   << ':';
  cout.width(2); cout.fill('0'); cout << minutes << ':';
  cout.width(2); cout.fill('0'); cout << seconds;
  cout << "\n\t\t\t--------------------------------\n\n\n";

//   int m, i, j;

//   double c;

//   char debfile[80];

//   m = 10;

//   valarray<double> vd(-2,m);
//   valarray<double> ve(1,m-1);

//   TriDiag A(m,1);

//   A.get_d() = vd;
//   A.get_e() = ve;
//   A.get_f() = ve;

//   A.get_e()[m-2] = 2.0;
//   A.get_f()[0] = 2.0;

//   MultipleRHS b(1,m);

//   for(i = 0; i < m; i++)
//     b[0][i] = 4*pow(-1.0,i);

//   sprintf(debfile,"/var/tmp/1D_Diff/rhs");
//   Write_MultipleRHS(debfile,b);
  
//   A.ChaseSolve(b);

//   sprintf(debfile,"/var/tmp/1D_Diff/solu_0");
//   Write_MultipleRHS(debfile,b);

//   c = -b[0].sum()/m;

//   b[0] += c;

//   sprintf(debfile,"/var/tmp/1D_Diff/solu_c");
//   Write_MultipleRHS(debfile,b);
}  



