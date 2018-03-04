// Trial runs for the Fitter
#include "Fitter.h"
#include <iostream>

using namespace std;
using namespace orbits;

int main(int argc,
	 char *argv[])
{
  try {
    Ephemeris eph;

    Fitter fit(eph);

    ifstream ifs(argv[1]);
    fit.readObservations(ifs);
    ifs.close();
    
    fit.chooseFrame();
  } catch (std::runtime_error& e) {
    quit(e);
  }
  exit(0);
}
