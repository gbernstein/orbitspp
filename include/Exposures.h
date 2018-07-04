// Manipulate information on DECam exposures
#ifndef EXPOSURES_H
#define EXPOSURES_H

#include "Astrometry.h"
#include "OrbitTypes.h"
#include "PlaneGeometry.h"
#include "Ephemeris.h"

namespace orbits {

  class Exposure {
  public:
    int expnum;
    double mjd;
    double tdb;
    astrometry::Vector3 earth;  // Observatory position in Frame
    Point axis;   // Optical axis coords, in Frame gnomonic, degrees

    double detectionDensity; // Density of transients, per sq degree

    // ??? Add per-CCD info on boundaries, detection lists.
  };

  extern
  std::vector<Exposure>
  selectExposures(const Frame& frame,   // Starting coordinates, time
		  const Ephemeris& ephem,  // Ephemeris to use
		  double gamma0,        // Center and width of range 
		  double dGamma,        // of gamma to cover
		  double searchRadius,  // Range of starting coords to cover
		  string exposureFile="data/y4a1.exposure.positions.fits",  // File of exposure data
		  double fieldRadius = 1.1); // Circumscribed field radius (degrees)
  // ?? Allow exclusion of filters??
  // ?? Add CCD corners, detections ??
  

} // end namespace orbits
#endif 
