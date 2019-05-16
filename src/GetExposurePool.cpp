#include <iostream>
#include "StringStuff.h"
#include "Astrometry.h"
#include "Exposures.h"
#include "Pset.h"

const string usage=
  "GetExposurePool: Return list of all DES exposures that could\n"
  "have observed a TNO in a given distance (gamma) range which\n"
  "was within given radius of a reference position at a given\n"
  "time.\n"
  "\n"
  "Usage: GetExposurePool <ra0> <dec0> <search radius> <tdb0> <gamma> <dgamma>\n"
  "                       [-parameter value] ...\n"
  "Inputs:\n"
  "   <ra0>,<dec0> = RA,Dec of reference point (degrees)\n"
  "   <search radius> = max distance from reference point at\n"
  "            reference time (degrees)\n"
  "   <tdb0> = reference time, (years past J2000.0)\n"
  "   <gamma> is central value of gamma=1/distance\n"
  "   <dgamma> is half-width of allowed range of gamma\n"
  "   [-parameter value] are any number of parameter/value pairs.\n"
  "            Possible parameters are listed below.\n"
  "            Note that exposure and ephemeris file locations\n"
  "            are read from environment variables if not specified.\n"
  "Output:\n"
  "   First line recaps inputs\n"
  "   Succeeding lines are <expnum> <dt> <astrometric?> for each possible\n"
  "   exposure, where\n"
  "   <expnum> = DECam exposure number\n"
  "   <dt>     = time from tdb0 (in years)\n"
  "   <astrometric> = 1 if exposure is part of astrometric solution, 0 otherwise";

using namespace std;
using namespace orbits;


int main(int argc,
	 char *argv[])
{
  string ephemerisPath;
  string exposurePath;
  int obsCode;
  Pset parameters;

  double ra0, dec0, tdb0;
  double searchRadius=0.;
  double gamma0=0., dGamma=0.;
  try {
    {
      const int def=PsetMember::hasDefault;
      const int low=PsetMember::hasLowerBound;
      const int up=PsetMember::hasUpperBound;
      const int lowopen = low | PsetMember::openLowerBound;
      const int upopen = up | PsetMember::openUpperBound;

      parameters.addMember("ephemerisFile",&ephemerisPath, def,
			   "SPICE file (null=>environment)", "");
      parameters.addMember("exposureFile",&exposurePath, def,
			   "DECam exposure table (null=>environment)", "");
      parameters.addMember("obscode",&obsCode, def | low,
			   "Observatory code (default: CTIO)", 807, 0);
    }
    parameters.setDefault();
    if (argc<7 || string(argv[1])=="-h" || string(argv[1])=="--help") {
      cout << usage << endl;
      cout << "\nParameters:" << endl;
      parameters.dump(cerr);
      exit(1);
    }

    ra0 = atof(argv[1]);
    dec0 = atof(argv[2]);
    searchRadius = atof(argv[3]);
    tdb0 = atof(argv[4]);
    gamma0 = atof(argv[5]);
    dGamma = atof(argv[6]);

    // Read any arguments from command line
    parameters.setFromArguments(argc-6, argv+6);

    Ephemeris ephem(ephemerisPath);

    // Build the reference frame
    Frame frame;
    {
      astrometry::SphericalICRS pole(ra0*DEGREE, dec0*DEGREE);
      frame.orient = astrometry::Orientation(pole);
      frame.orient.alignToEcliptic();
      frame.origin = ephem.observatory(obsCode, tdb0);
      frame.tdb0 = tdb0;
    }

    // Output configuration info:
    cout << "# RA0  Dec0  searchRadius TDB0 gamma0  dGamma" << endl;
    cout << fixed << setprecision(7) << ra0
	 << " " << showpos << dec0
	 << " " << noshowpos << setprecision(3) << searchRadius
	 << " " << setprecision(7) << tdb0
	 << " " << setprecision(6) << gamma0
	 << " " << dGamma
	 << endl;
    cout << "# expnum  TDB-TDB0" << endl;

    ExposureTable et(exposurePath);

    // Find exposures that potentially contain the object
    auto possibleExposures = et.getPool(frame, ephem,
					gamma0, dGamma, searchRadius*DEGREE,
					false); // Get non-astrometric exposures too

    // Output results
    for (auto eptr : possibleExposures) {
      cout << noshowpos << setw(6) << eptr->expnum
	   << " " << setprecision(6) << showpos << setw(9) << eptr->tobs
	   << " " << noshowpos << eptr->astrometric << endl;
      delete eptr;
    }
    
  } catch (std::runtime_error& e) {
    quit(e);
  }
  exit(0);
}

