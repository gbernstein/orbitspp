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

    //Fitter fit(eph, Gravity::BARY);
    Fitter fit(eph, Gravity::GIANTS);
    fit.setBindingConstraint(1.);
    
    ifstream ifs(argv[1]);
    fit.readMPCObservations(ifs);
    ifs.close();

    fit.chooseFrame(-1);

    fit.setLinearOrbit();
    cerr << "First ABG:";
    fit.abg.writeTo(cerr);
    cerr << "distance: " << 1./fit.abg[ABG::G] << endl;

    fit.newtonFit();
    fit.printResiduals(cerr);
    cerr << "Final ABG:\n";
    fit.abg.writeTo(cerr);
    cerr << endl << "Chisq: " << fit.getChisq() << endl;
    cerr << "distance: " << 1./fit.abg[ABG::G] << endl;
    fit.printCovariance(cerr);
    cerr << endl;
  } catch (std::runtime_error& e) {
    quit(e);
  }
  exit(0);
}
