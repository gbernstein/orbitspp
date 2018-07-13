// Run a simple fit on data given in degrees, from CTIO
#include "Fitter.h"
#include "Elements.h"
#include "StringStuff.h"

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

    string buffer;
    while (stringstuff::getlineNoComment(cin, buffer)) {
      std::istringstream iss(buffer);
      double mjd, ra, dec;
      iss >> mjd >> ra >> dec;
      Observation obs;
      obs.radec = astrometry::SphericalICRS(ra*DEGREE, dec*DEGREE);
      obs.tdb = eph.mjd2tdb(mjd);
      obs.observer = eph.observatory(807, obs.tdb);
      obs.cov(0,0) = obs.cov(1,1) = pow(0.015*ARCSEC,2.);
      obs.cov(1,0) = obs.cov(0,1) = 0.;
      fit.addObservation(obs);
    }

    fit.chooseFrame(-1);
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

    cerr << "Reference frame: " << endl;
    Frame::writeHeader(cerr);
    fit.getFrame().write(cerr);
    
    fit.setLinearOrbit();
    cerr << "First ABG:" << endl;
    ABG::writeHeader(cerr);
    cerr << fit.getABG(true) << endl;

    cerr << "distance: " << 1./fit.getABG(true)[ABG::G] << endl;

    /** {
      fit.getABG(true)[0] = -0.0125341;
      fit.getABG(true)[1] = 0.0131069;
      fit.getABG(true)[2] = 0.0276466;
      fit.getABG(true)[3] = 0.0209088;
      fit.getABG(true)[4] = -0.0218552;
      fit.getABG(true)[5] = 0.00432035;
      cerr << "Fake ABG:";
      fit.getABG().writeTo(cerr);
  } /**/
    
    fit.newtonFit();
    fit.printResiduals(cerr);
    cerr << "Final ABG:\n";
    ABG::writeHeader(cerr);
    cerr << fit.getABG(true) << endl;
    {
      ABGCovariance cov = fit.getInvCovarABG().inverse();
      cov.write(cerr);
    }
    cerr << "Chisq: " << fit.getChisq() << endl;
    cerr << "distance: " << 1./fit.getABG()[ABG::G] << endl;

    // Calculate orbital elements:
    cerr << "Elements: " << endl;
    Elements::writeHeader(cerr);
    cerr << fit.getElements() << endl;
    {
      auto cov =  fit.getElementCovariance();
      cov.write(cerr);
    }
    
    
  } catch (std::runtime_error& e) {
    quit(e);
  }
  exit(0);
}
