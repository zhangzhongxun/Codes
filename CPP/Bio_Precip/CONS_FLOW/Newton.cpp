#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <valarray>
#include <limits>

using namespace std;

//Nonlinear function for [H+] (F([H+] = 0, derived from charge neutral and chemical reaction equilibrium, 
//need to be solved by Newton's method)
double NLF_H(double H, double NT, double K3, double Kw, double CT, double K1, double K2, double Ca, double Cl){

  double result, tmp0, tmp1, tmp2;

  tmp0 = 1.0 + K3*H; 
  tmp1 = H/K1 + 1.0 + K2/H;
  tmp2 = H*H/K1/K2 + H/K2 + 1.0;

  result = NT*K3*H/tmp0 + H - Kw/H - CT/tmp1 - 2.0*CT/tmp2 + 2.0*Ca - Cl;

  return result;
}

//Derivative of the previous nonlinear function (NLF_H), used in Newton's method
double Deri_NLF_H(double H, double NT, double K3, double Kw, double CT, double K1, double K2, double Ca, double Cl){

  double result, tmp0, tmp1, tmp2;
  
  tmp0 = 1.0 + K3*H; 
  tmp1 = H/K1 + 1.0 + K2/H;
  tmp2 = H*H/K1/K2 + H/K2 + 1.0;

  result = NT*K3/tmp0/tmp0 + 1.0 + Kw/H/H + CT*(1.0/K1 - K2/H/H)/tmp1/tmp1 + 2.0*CT*(2.0*H/K1/K2 + 1.0/K2)/tmp2/tmp2;

  return result;
}

//Function to find the root of NLF_H, parameter H holds initial guess at input, and the root as output
//Return 0 if successful, 1 not successful

int Newton_Root(double& H, double NT, double K3, double Kw, double CT, double K1, double K2, double Ca, double Cl, 
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
      tmp1 = NLF_H(x_n, NT, K3, Kw, CT, K1, K2, Ca, Cl);
      tmp2 = Deri_NLF_H(x_n, NT, K3, Kw, CT, K1, K2, Ca, Cl); 
 
      if(fabs(tmp1) < deps) //already find a solution
	{
	  H = x_n;
          N_iter = 0;
          cout << endl << "Exit before iteration! " << endl;
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
 
      tmp1 = NLF_H(x_n, NT, K3, Kw, CT, K1, K2, Ca, Cl);
      tmp2 = Deri_NLF_H(x_n, NT, K3, Kw, CT, K1, K2, Ca, Cl); 

      for(n = 0; n < 20; n++)
	{
	  x_np1 = x_n - tmp1/tmp2;

	  if(fabs((x_np1 - x_n)/x_n) < 1e-8) //relative error small, find a solution!
	    {
	      break;
	    }

	  x_n = x_np1;

	  tmp1 = NLF_H(x_n, NT, K3, Kw, CT, K1, K2, Ca, Cl);
	  tmp2 = Deri_NLF_H(x_n, NT, K3, Kw, CT, K1, K2, Ca, Cl); 
	}

      ini_N -= 1.0;
    }while(x_np1 <= 0);

  N_iter = n;
  H = x_np1;
  tmp1 = NLF_H(x_np1, NT, K3, Kw, CT, K1, K2, Ca, Cl);
  cout << endl << "the residual is " << tmp1 << ", x_np1 = " << x_np1 << endl ;
 
  if(n < 20)
    return 0;
  else 
    return 1;
}

//function to calculate the concentration of species participating fast reactions
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

int main(){

  double H, NT, K3, Kw, CT, K1, K2, Ca, Cl, FH, DFH, ADFH, dH, NH3, NH4, H2CO3, HCO3, CO3;

  int flag, N;

  //here use the unit milli-mole per liter, then PH = -log([H+]*1e-3)

  H = 1e-6;
  NT = 1.0e-5; //4e-4;
  K3 = 1e9; //1e-9;
  Kw = 1e-14; //1e-14;
  CT = 1e-5; //0e-7;
  K1 = 1e-6; //1e-6;
  K2 = 1e-10; //1e-10;
  Ca = 1e-6; //1e-5;
  Cl = 2e-6; //2e-5;
  dH = 1e-10;

//   FH = NLF_H(H, NT, K3, Kw, CT, K1, K2, Ca, Cl);
//   DFH = Deri_NLF_H(H, NT, K3, Kw, CT, K1, K2, Ca, Cl);
//   ADFH = (NLF_H(H + .5*dH, NT, K3, Kw, CT, K1, K2, Ca, Cl) - NLF_H(H-.5*dH, NT, K3, Kw, CT, K1, K2, Ca, Cl))/dH;

//   cout << endl << "At H = " << H << ", F(H) = " << FH << ", and the derivative is F'(H) = " << DFH
//        << endl << "the approximation of the derivative with dH = " << dH << " is " << ADFH << endl;

//   cout << endl << "The machine epsilon for double is " << numeric_limits<double>::epsilon() << endl;

//   cout << endl << "(1.2e-100 - 1e-100)/1e-100 = " << (1.2e-100 - 1e-100)/1e-100 << endl;

  flag = Newton_Root(H, NT, K3, Kw, CT, K1, K2, Ca, Cl, N);

  cout << endl << "The Newton's method returned after " << N << " steps, flag = " << flag 
       << ", H = " << setprecision(14) << H << endl;

  cout << endl << "F(H) = " << NLF_H(H, NT, K3, Kw, CT, K1, K2, Ca, Cl) 
       << ", the PH is " << -log10(H) << endl;
    
  Cal_Fast_Species(H, NT, K3, Kw, CT, K1, K2, NH3, NH4, H2CO3, HCO3, CO3);

  cout << endl << " NH3 = " << NH3 << ", NH4 = " << NH4 << " NH3 + NH4 = " << NH3 + NH4 << ", NT = " << NT
       << endl << " H2CO3 = " << H2CO3 << ", HCO3 = " << HCO3 << ", CO3 = " << CO3 << ", H2CO3 + HCO3 + CO3 = "
       << H2CO3 + HCO3 + CO3 << ", CT = " << CT << endl;
       
  
}
