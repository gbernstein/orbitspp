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
    ABG(): Base(0.) { (*this)[G]=0.03;}

    void getState(double t,
		  astrometry::Vector3& x,
		  astrometry::Vector3& v) const {
      // Calculate state vector in reference system at given time for inertial motion
      x[0] = (*this)[ADOT]*t + (*this)[A];
      x[1] = (*this)[BDOT]*t + (*this)[B];
      x[2] = (*this)[GDOT]*t + 1.;
      x /= (*this)[G];
      v[0] = (*this)[ADOT] / (*this)[G];
      v[1] = (*this)[BDOT] / (*this)[G];
      v[2] = (*this)[GDOT] / (*this)[G];
      return;
    }

    astrometry::DMatrix getXYZ(const astrometry::DVector& t) const;
      // Calculate inertial XYZ positions at an array of times, return n x 3 matrix

    void writeTo(std::ostream& os, int precision=6) const;
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
    // This inherited class adds bulk coordinate transforms via Eigen.
  public:
    Frame(const astrometry::CartesianCoords& origin,
	  const astrometry::Orientation& orient,
	  const double tdb): ReferenceFrame(origin,orient),
			     tdb0(tdb) {}
    Frame(): ReferenceFrame(), tdb0(0.) {}
    double tdb0;  // Time coordinate origin

    // Transform positions, or a 3xN array of positions,
    // to or from this Frame from or to ICRS.
    // isVelocity=true will 
    astrometry::Vector3 toICRS(const astrometry::Vector3& x,
			       bool isVelocity=false) const;
    astrometry::Vector3 fromICRS(const astrometry::Vector3& x,
			       bool isVelocity=false) const;
    astrometry::DMatrix toICRS(const astrometry::DMatrix& x,
			       bool isVelocity=false) const;
    astrometry::DMatrix fromICRS(const astrometry::DMatrix& x,
			       bool isVelocity=false) const;
  };

} // end namespace orbits

#endif
