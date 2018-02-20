// Interface to JPL Horizons ephemeris via the CSPICE code
#ifndef EPHEMERIS_H
#define EPHEMERIS_H

#include "Astrometry.h"
#include "OrbitTypes.h"

namespace orbits {
  class SpiceError: public std::runtime_error {
  public:
    SpiceError(const string &m=""): 
      std::runtime_error("SPICE Error: " +m) {}
  };
    
  class Ephemeris {
  public:
    Ephemeris(const string& kernelFile="");   // Will read specified kernel, or a default
    // No copying since cspice has its own internal state
    Ephemeris(const Ephemeris& rhs) =delete;
    void operator=(const Ephemeris& rhs) =delete;

    // Get barycentric body geometric position
    astrometry::CartesianICRS position(int body,
				       double tdb) const;
    // Get barycentric body position and velocity
    State state(int body,
		double tdb) const;

    // Get barycentric coordinate of Earth-fixed observatory
    // from geodetic lat, lon, elevation
    astrometry::CartesianICRS observatory(double lon,
					  double lat,
					  double elev,
					  double tdb) const;
    // Or from an MPC observatory code
    astrometry::CartesianICRS observatory(int obscode,
					  double tdb) const;

    // Return barycentric dynamical time (TDB) from UTC.
    // The former is the time argument of planetary equations of motion and
    // of the DE430 ephemeris.  The latter roughly tracks Earth rotation, so
    // differs by leap seconds, and by some gravitational redshift terms.
    // Can take as input either a string in valid CSPICE format, or my UT object.
    // TDB is defined here as Julian years past J2000 TDB.
    double utc2tdb(const string utc) const;
    double utc2tdb(const astrometry::UT& utc) const;
    // Convert UTC-based JD or MJD to TDB.
    double jd2tdb(double jd) const;
    double mjd2tdb(double mjd) const;
    
    // Here are standard object ID's we will use:
    const int SSB = 0; //Solar system barycenter
    const int SUN = 10; // Heliocenter
    const int EARTH = 399; // Geocenter
    const int JUPITER = 5; // Jupiter system barycenter
    const int SATURN = 6; // Saturn barycenter
    const int URANUS = 7; // Uranus barycenter
    const int NEPTUNE = 8; // Neptune barycenter
    
  private:
    static bool init;  // We will only allow one Ephemeris to exist
    // Some static constants we can get from spice
    static double earthRadius;
    static double earthFlatten;
  };

} // end namespace orbits


#endif
