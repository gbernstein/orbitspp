// Compare predictions of orbits by me to those from Gerdes/pyephem

#include <map>
#include <cctype>
#include "AstronomicalConstants.h"
#include "Astrometry.h"
#include "OrbitTypes.h"
#include "Elements.h"
#include "Ephemeris.h"
#include "Trajectory.h"
#include "Fitter.h"

#include "StringStuff.h"

using namespace std;
using namespace orbits;
using namespace astrometry;

string usage =
  "FitBulk: fit orbits to a collection of observations of multiple\n"
  "objects from a fixed site.\n"
  "<id>  <mjd>  <ra>  <dec>\n"
  "on each line (blank lines and # lines are ignored).\n"
  "Usage:\n"
  "FitBulk [obscode] [distance] [distance_range]\n"
  " obscode: MPC observatory code for observations (default: 807=CTIO)\n"
  " distance: optional central value for distance prior (default: none)\n"
  " distance_range: fractional 1-sigma width of distance prior, if used (default: 0.1)\n"
  " Input:   lines containing\n"
  "          <id>  <mjd>  <ra>  <dec> <sigma_in_arcsec>\n"
  "          ra and dec are in decimal degrees. Blank lines and # lines are ignored.\n"
  " Output: Orbit fit info per unique ID";
  

int
main(int argc,
     char *argv[])
{
  if (argc>5 || (argc>1 && !std::isdigit(argv[1][0]))) {
    cerr << usage << endl;
    exit(1);
  }
  try {
    // Create ephemeris
    Ephemeris eph;

    // Parse command line args
    int obscode = 807;
    if (argc>1)
      obscode = atoi(argv[1]);
    double gammaPrior = 0.;
    if (argc>2)
      gammaPrior = 1. / atof(argv[2]);
    double gammaPriorSigma = 0.1 * gammaPrior;
    if (argc>3)
      gammaPrior = atof(argv[3]) * gammaPrior;
    
    // Read all observations
    string line;
    typedef multimap<int, Observation*> ObsMap;
    ObsMap obsmap;
    
    while (stringstuff::getlineNoComment(cin,line)) {
      int id;
      float mjd, ra, dec, sigma;
      istringstream iss(line);
      if (!(iss >> id >> mjd >> ra >> dec >> sigma)) {
	cerr << "Bad input line: <" << line << ">" << endl;
	exit(1);
      }

      auto obs = new Observation;
      obs->radec = SphericalICRS(ra*DEGREE, dec*DEGREE);
      obs->tdb = eph.mjd2tdb(mjd);
      obs->observer = eph.observatory(obscode,obs->tdb);
      obs->cov(0,0) = obs->cov(1,1)= pow(sigma*ARCSEC,2);
      obs->cov(0,1) = obs->cov(1,0) = 0.;

      obsmap.insert(ObsMap::value_type(id,obs));
    }

    // Now we're going to enter a loop where we fit each object.

    for (auto obsptr = obsmap.begin(); obsptr != obsmap.end(); ) {
      // Make a new fitter and give it the observations for this ID
      Fitter fit(eph, Gravity::GIANTS);
      if (gammaPrior > 0.)
	fit.setGammaConstraint(gammaPrior,gammaPriorSigma);
      
      int id = obsptr->first;
      for (auto ptr = obsptr; ptr!=obsmap.upper_bound(id); ++ptr) {
	fit.addObservation(*(ptr->second));
      }
      obsptr = obsmap.upper_bound(id);

      // Execute a fit
      fit.chooseFrame(1);
      fit.setLinearOrbit();
      fit.setBindingConstraint(2.); //****
      try {
	fit.newtonFit();
      } catch (std::runtime_error& m) {
	cout << "ID " << id << " FAILURE***" << endl;
	try { fit.newtonFit(0.01, true);} catch (std::runtime_error& m) {}
      }
      cout << "ID " << id
	   << " chisq " << fit.getChisq()
	   << " distance " << 1./fit.getABG()[ABG::G]
	   << " energy " <<  (pow(fit.getABG()[ABG::ADOT],2.)+
			     pow(fit.getABG()[ABG::BDOT],2.)+
			      pow(fit.getABG()[ABG::GDOT],2.)) / (2*GM * pow(fit.getABG()[ABG::G],3.)) - 1
			     << endl;
      fit.printResiduals(cout);
      cout << endl;
    }

  
  } catch (std::runtime_error& m) {
    quit(m);
  }
  exit(0);
}
