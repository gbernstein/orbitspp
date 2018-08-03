// Testing selection of exposures from list that might contain orbits
// Right now just prints out a list of the accepted exposures' times and positions

#include "Ephemeris.h"
#include "AstronomicalConstants.h"
#include "Exposures.h"
#include <iostream>
#include <iomanip>
using namespace std;

using namespace orbits;
using namespace astrometry;

int
main(int argc,
     char *argv[])
{
  bool fail = false;  // Flag any failures

  if (argc!=3) {
    cerr << "Usage: testExposureSelection <gamma0> <dGamma>\n"
      " <gamma0> is middle of range of gamma for possible objects\n"
      " <dGamma> is full width of gamma range\n"
      " Returns: Table of <expnum> <dt> <x> <y> relative to exposure 387104"
	 << endl;
    exit(1);
  }
  double gamma0 = atof(argv[1]);
  double dGamma = atof(argv[2]);
    
  try {
    orbits::Ephemeris ephem;
  
    // Get exposures connected to some particular
    // place/time: - use center of exposure 387104
    double tdb0 = ephem.mjd2tdb(57004.18198925);
    SphericalICRS pole(31.38881405*DEGREE, -15.42717461*DEGREE);
    Orientation orient(pole);
    orient.alignToEcliptic();
    CartesianICRS origin(0.1591462 ,  0.8910222 ,  0.38608784);
    
    Frame frame(origin, orient, tdb0);
    
    ExposureTable et;
    auto answer = et.getPool(frame, ephem,
			     gamma0, dGamma, 0.5,
			     false);   // Accept non-astrometric exposures

    for (auto& expo : answer) {
      cout << expo->expnum << " " << expo->tdb-tdb0
	   << " " << expo->axis[0] << " " << expo->axis[1]
	   << " " << expo->astrometric
	   << endl;
    }
  } catch (std::runtime_error& e) {
    quit(e);
  }
  
  exit(fail ? 1 : 0);
}

