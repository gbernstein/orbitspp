// Interface to JPL Horizons ephemeris via the CSPICE code
#ifndef EPHEMERIS_H
#define EPHEMERIS_H

#include "Astrometry.h"
#include "OrbitTypes.h"
#include <map>
#include <vector>
#include "SharedLUT.h"

namespace orbits {
  // Here are standard spice object ID's we will use:
  const int SSB = 0; //Solar system barycenter
  const int SUN = 10; // Heliocenter
  const int EARTH = 399; // Geocenter
  const int JUPITER = 5; // Jupiter system barycenter
  const int SATURN = 6; // Saturn barycenter
  const int URANUS = 7; // Uranus barycenter
  const int NEPTUNE = 8; // Neptune barycenter

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

    ~Ephemeris() { 
      // Delete the caches 
      for (auto& pr : caches)
    	delete pr.second;
    }
    
    // Get barycentric body geometric position
    // Last option false to force call to SPICE.
    astrometry::CartesianICRS position(int body,
				       double tdb,
				       bool useCache=true) const;
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
    double tdb2jd(double tdb) const;
    double tdb2mjd(double tdb) const;
    
    // Cache the positions of chosen body at chosen interval/resolution,
    // so it will be obtained without having to single-thread through SPICE.
    void cachePositions(int body, double dt=10.*DAY) {
      caches[body] = new Cache();
      caches[body]->dt = dt;
    }
    
  private:
    static bool init;  // We will only allow one Ephemeris to exist
    // Some static constants we can get from spice
    static double earthRadius;
    static double earthFlatten;

    // Maintain an internal table of observatories
    struct ObsInfo {
      double lon;  // geodetic Lon, lat in radians
      double lat;
      double elev; // geodetic elevation, meters
    };
    mutable std::map<int, ObsInfo> obsTable;

    struct Cache {
      double dt;
      mutable SharedLUT<Vector3, Eigen::aligned_allocator<Vector3>> lut;
      Vector3 operator[](int i) const {return lut[i];}
    };
      
    mutable std::map<int, Cache*> caches;
  };

  

  // Function to return a standard observation given MPC info
  extern
  Observation mpc2Observation(const MPCObservation& mpc, const Ephemeris& ephem);
    
} // end namespace orbits


#endif
