#include "Driver.h"

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

  time(&start);


  Time_Evolve(argv[1]);

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
}
