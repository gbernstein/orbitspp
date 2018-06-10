/** Classes/structures for orbit fitting **/
#ifndef ORBITTYPES_H
#define ORBITTYPES_H

#include "AstronomicalConstants.h"

#include "LinearAlgebra.h"
#include "Astrometry.h"

namespace orbits {

  class State {
    /* State vector for solar system body */
  public:
    astrometry::CartesianICRS x;   // Position, AU, barycentric, nominally ICRS
    astrometry::CartesianICRS v;   // Velocity, AU per Julian year
    double tdb;                    // TDB for which it applies, Julian yrs since J2000
  };

  class ABG: public linalg::SVector<double,6> {
    /* Alpha-beta-gamma basis used for Bernstein-Khushalani */
  public:
    // Indices for components:
    static const int A=0;
    static const int B=1;
    static const int G=2;
    static const int ADOT=3;
    static const int BDOT=4;
    static const int GDOT=5;
    typedef linalg::SVector<double,6> Base;
    ABG(const State& s) {
      (*this)[A] = s.x[0]/s.x[2];
      (*this)[B] = s.x[1]/s.x[2];
      (*this)[G] = 1./s.x[2];
      (*this)[ADOT] = s.v[0]/s.x[2];
      (*this)[BDOT] = s.v[1]/s.x[2];
      (*this)[GDOT] = s.v[2]/s.x[2];
    }
    ABG(): linalg::SVector<double,6>(0.) { (*this)[G]=0.03;}
    
  };

  class Elements: public linalg::SVector<double,6> {
    // Keplerian elements for elliptical orbit
  public:
    static const int A=0;
    static const int E=1;
    static const int I=2;
    static const int LAN=3; // Longitude of ascending node
    static const int AOP=4; // Argument of perihelion
    static const int TOP=5; // Time of perihelion passage
  };

  class Observation {
    // General observer-frame information about a detection
  public:
    EIGEN_NEW
    Observation() {}
    astrometry::SphericalICRS radec;
    double tdb;  //  TDB of observation, referenced to J2000
    astrometry::Matrix22 cov; // RA/Dec error covariance matrix, radians
    astrometry::CartesianICRS observer;  // ICRS barycentric position of observer
  };
    
  class MPCObservation {
    // Information that comes from a standard MPC line
  public:
    MPCObservation() {}
    MPCObservation(const string& line); // Parse MPC style input string
    astrometry::SphericalICRS radec;
    double mjd;  // UTC MJD of observation
    double sigma; // Position uncertainty per component, arcsec
    int obscode; // MPC observatory code
  };

  typedef linalg::SMatrix<double,6,6> ABGCovar;
  
  class Frame: public astrometry::ReferenceFrame {
    // Coordinate reference frame, including
    // 3d coords of origin,
    // orientation of axes,
    // and TDB of temporal origin.
  public:
    Frame(const astrometry::CartesianCoords& origin,
	  const astrometry::Orientation& orient,
	  const double tdb): ReferenceFrame(origin,orient),
			     tdb0(tdb) {}
    Frame(): ReferenceFrame(), tdb0(0.) {}
    double tdb0;  // Time coordinate origin
  };

} // end namespace orbits

#endif
