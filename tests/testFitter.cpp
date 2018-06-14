// Trial runs for the Fitter
#include "Fitter.h"
#include "Elements.h"

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

    //**fit.chooseFrame(-1);
    //**fit.chooseFrame(0);
    /**{
      string line = "56877.1624158 -3:03:03.9501 -41:23:03.0142 0.1 807";
      auto obs = mpc2Observation(MPCObservation(line), eph);
      astrometry::Orientation orient(obs.radec);
      orient.alignToEcliptic();
      double tdb0 = 14.;
      astrometry::CartesianICRS origin = obs.observer;
      Frame f(origin,orient,tdb0);
      fit.setFrame(f);
      } /**/

    fit.setLinearOrbit();
    cerr << "First ABG:";
    fit.abg.writeTo(cerr);
    cerr << "distance: " << 1./fit.abg[ABG::G] << endl;

    /** {
      fit.abg[0] = -0.0125341;
      fit.abg[1] = 0.0131069;
      fit.abg[2] = 0.0276466;
      fit.abg[3] = 0.0209088;
      fit.abg[4] = -0.0218552;
      fit.abg[5] = 0.00432035;
      cerr << "Fake ABG:";
      fit.abg.writeTo(cerr);
  } /**/
    
    fit.newtonFit();
    fit.printResiduals(cerr);
    cerr << "Final ABG:\n";
    fit.abg.writeTo(cerr);
    cerr << endl << "Chisq: " << fit.getChisq() << endl;
    cerr << "distance: " << 1./fit.abg[ABG::G] << endl;
    fit.printCovariance(cerr);
    cerr << endl;

    // Calculate orbital elements:
    cerr << "Elements: " << endl;
    cerr << fit.getElements();
    cerr << endl;
    cerr << "Covariance: " << endl;
    cerr << fit.getElementCovariance();
    cerr << endl;
    
    
  } catch (std::runtime_error& e) {
    quit(e);
  }
  exit(0);
}
