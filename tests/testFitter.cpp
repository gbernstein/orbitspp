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
    fit.readMPCObservations(ifs);
    ifs.close();

    /**/cerr << "Ready to choose" << endl;
    
    fit.chooseFrame(-1);

    /**
    fit.setLinearOrbit();
    cerr << "ABG:\n" << fit.abg << endl;
    cerr << "distance: " << 1./fit.abg[ABG::G] << endl;
    **/
  } catch (std::runtime_error& e) {
    quit(e);
  }
  exit(0);
}
