//function to test the calculation of species concentration after fast reaction

#include "Build.h"
#include "DataStruct.h"
#include "IO.h"

int main(int argc, char* argv[]){

  if(argc < 2)
    {
      cout << endl << "Add the data file after command " << argv[0] << ", Data Missed !!!" << endl;
      exit(-1);
    }

  char Ini_file[5][80];
  char datacomment[100];
  char debfile[100];
  char SaveD[100] = "/var/tmp/TFast/";

  ifstream ifile(argv[1],ios::in);

  if(!ifile)
    {
      cout << endl << "Can not open input file " << argv[1] << ", ERROR at line " << __LINE__ << " of file " << __FILE__  << endl;
      exit(-1);
    } 

  ifile >> datacomment >> Ini_file[0]
 	>> datacomment >> Ini_file[1]
  	>> datacomment >> Ini_file[2]
  	>> datacomment >> Ini_file[3]
	>> datacomment >> Ini_file[4];
		 
  ifile.close();

  int Nx, Ny;

  double K1, K2, K3, Kw;

  K1 = 1e-6;
  K2 = 1e-10;
  K3 = 1e9;
  Kw = 1e-14;

  Nx = 255;
  Ny = 64;

  Node_Data H_n(Nx, Ny);
  Node_Data NT(Nx, Ny);
  Node_Data CT(Nx, Ny);
  Node_Data Ca(Nx, Ny);
  Node_Data Cl(Nx, Ny);
  Node_Data NH3(Nx, Ny);
  Node_Data NH4(Nx, Ny);
  Node_Data H2CO3(Nx, Ny);
  Node_Data HCO3(Nx, Ny);
  Node_Data CO3(Nx, Ny);

  ifile.open(Ini_file[0],ios::in);
  if(!ifile)
    {
      cout << endl << " Can not open input file " << Ini_file[0] << ", ERROR at line " << __LINE__ << " of file " << __FILE__ << endl;
      exit(-1);
    }

  ifile >> H_n;
  ifile.close();

  ifile.open(Ini_file[1],ios::in);
  if(!ifile)
    {
      cout << endl << " Can not open input file " << Ini_file[1] << ", ERROR at line " << __LINE__ << " of file " << __FILE__ << endl;
      exit(-1);
    }

  ifile >> NT;
  ifile.close();


  ifile.open(Ini_file[2],ios::in);
  if(!ifile)
    {
      cout << endl << " Can not open input file " << Ini_file[2] << ", ERROR at line " << __LINE__ << " of file " << __FILE__ << endl;
      exit(-1);
    }

  ifile >> CT;
  ifile.close();



  ifile.open(Ini_file[3],ios::in);
  if(!ifile)
    {
      cout << endl << " Can not open input file " << Ini_file[3] << ", ERROR at line " << __LINE__ << " of file " << __FILE__ << endl;
      exit(-1);
    }

  ifile >> Ca;
  ifile.close();



  ifile.open(Ini_file[4],ios::in);
  if(!ifile)
    {
      cout << endl << " Can not open input file " << Ini_file[4] << ", ERROR at line " << __LINE__ << " of file " << __FILE__ << endl;
      exit(-1);
    }

  ifile >> Cl;
  ifile.close();

  //Calculate the concentration of species after the fast reactions
  Cal_Fast_Species_Node(H_n, NT, K3, Kw, CT, K1, K2, Ca, Cl, NH3, NH4, H2CO3, HCO3, CO3);


  sprintf(debfile,"%sH",SaveD);
  Write_Node_Data(debfile, H_n);

  sprintf(debfile,"%sNH3",SaveD);
  Write_Node_Data(debfile, NH3);


  sprintf(debfile,"%sNH4",SaveD);
  Write_Node_Data(debfile, NH4);

  sprintf(debfile,"%sH2CO3",SaveD);
  Write_Node_Data(debfile, H2CO3);


  sprintf(debfile,"%sHCO3",SaveD);
  Write_Node_Data(debfile, HCO3);

  sprintf(debfile,"%sCO3",SaveD);
  Write_Node_Data(debfile, CO3);

  sprintf(debfile,"%sNT",SaveD);
  Write_Node_Data(debfile, NT);

  sprintf(debfile,"%sCT",SaveD);
  Write_Node_Data(debfile, CT);


  sprintf(debfile,"%sCa",SaveD);
  Write_Node_Data(debfile, Ca);

  sprintf(debfile,"%sCl",SaveD);
  Write_Node_Data(debfile, Cl);

}
