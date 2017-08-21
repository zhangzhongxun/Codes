//Fuctions to set the initial conditions of phi_b, phi_c, Ur, Ca, OH, CT, c to be 
//a given constant, used for supplement when some of the above variables should be read in 
//from previous running

#include <iostream>
#include <fstream>
#include <cmath>
#include <valarray>
#include <iomanip>

using namespace std;

int main(int argc, char* argv[]){

     if(argc < 2)
    {
      cout << endl << "Add the data file after command " << argv[0] << ", Data Missed !!!" << endl;
      exit(-1);
    }

     ifstream ifile(argv[1],ios::in);

     if(!ifile)
       {
	 cout << endl << "Can not open input file" << argv[1] << ", Error in at line " << __LINE__ << " of file " << __FILE__ << endl;
	 exit(-1);
       }

     int Nx, Ny, i, j;

     double phb, phc, Ur, Ca, OH, CT, c;

     char datacomment[100];

     ifile >> datacomment >> Nx
           >> datacomment >> Ny
           >> datacomment >> phc
           >> datacomment >> phb
           >> datacomment >> Ur
           >> datacomment >> Ca
           >> datacomment >> OH
           >> datacomment >> CT
           >> datacomment >> c;

     ifile.close();

     ofstream ofile;
 
     ofile.open("ini__phc", ios::out);
     if(!ofile)
       {
	 cout << endl << "Can not open file ini__phc, ERROR at line " << __LINE__ << " of file " << __FILE__ << endl;
	 exit(-1);
       }

     for(j = 0; j <= Ny; j++)
       {
	 for(i = 0; i <= Nx + 1; i++)
           ofile << phc   << "  ";
	 ofile << endl;
       }
 
     ofile << endl;

     ofile.close();

     ofile.open("ini__phb", ios::out);
     if(!ofile)
       {
	 cout << endl << "Can not open file ini__phb, ERROR at line " << __LINE__ << " of file " << __FILE__ << endl;
	 exit(-1);
       }

     for(j = 0; j <= Ny; j++)
       {
	 for(i = 0; i <= Nx + 1; i++)
           ofile << phb   << "  ";
	 ofile << endl;
       }
 
     ofile << endl;

     ofile.close();


     ofile.open("ini__Ur", ios::out);
     if(!ofile)
       {
	 cout << endl << "Can not open file ini__Ur, ERROR at line " << __LINE__ << " of file " << __FILE__ << endl;
	 exit(-1);
       }

     for(j = 0; j <= Ny; j++)
       {
	 for(i = 0; i <= Nx + 1; i++)
           ofile << Ur   << "  ";
	 ofile << endl;
       }
 
     ofile << endl;

     ofile.close();


     ofile.open("ini__OH", ios::out);
     if(!ofile)
       {
	 cout << endl << "Can not open file ini__OH, ERROR at line " << __LINE__ << " of file " << __FILE__ << endl;
	 exit(-1);
       }

     for(j = 0; j <= Ny; j++)
       {
	 for(i = 0; i <= Nx + 1; i++)
           ofile << OH   << "  ";
	 ofile << endl;
       }
 
     ofile << endl;

     ofile.close();


     ofile.open("ini__Ca", ios::out);
     if(!ofile)
       {
	 cout << endl << "Can not open file ini__Ca, ERROR at line " << __LINE__ << " of file " << __FILE__ << endl;
	 exit(-1);
       }

     for(j = 0; j <= Ny; j++)
       {
	 for(i = 0; i <= Nx + 1; i++)
           ofile << Ca   << "  ";
	 ofile << endl;
       }
 
     ofile << endl;

     ofile.close();


     ofile.open("ini__CT", ios::out);
     if(!ofile)
       {
	 cout << endl << "Can not open file ini__CT, ERROR at line " << __LINE__ << " of file " << __FILE__ << endl;
	 exit(-1);
       }

     for(j = 0; j <= Ny; j++)
       {
	 for(i = 0; i <= Nx + 1; i++)
           ofile << CT   << "  ";
	 ofile << endl;
       }
 
     ofile << endl;

     ofile.close();

    

     ofile.open("ini__c", ios::out);
     if(!ofile)
       {
	 cout << endl << "Can not open file ini__c, ERROR at line " << __LINE__ << " of file " << __FILE__ << endl;
	 exit(-1);
       }

     for(j = 0; j <= Ny; j++)
       {
	 for(i = 0; i <= Nx + 1; i++)
           ofile << c   << "  ";
	 ofile << endl;
       }
 
     ofile << endl;

     ofile.close();


}











