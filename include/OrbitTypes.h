/** Classes/structures for orbit fitting **/
#ifndef ORBITTYPES_H
#define ORBITTYPES_H

#include "AstronomicalConstants.h"

#include "LinearAlgebra.h"
#include "Astrometry.h"

namespace orbits {

  // Bring standard vectors, matrices into namespace
  typedef linalg::SVector<double,2> Vector2;
  typedef linalg::SVector<double,3> Vector3;
  typedef linalg::SVector<double,6> Vector6;
  typedef linalg::SMatrix<double,2,2> Matrix22;
  typedef linalg::SMatrix<double,3,3> Matrix33;
  typedef linalg::SMatrix<double,6,6> Matrix66;
  typedef linalg::Vector<double> DVector;
  typedef linalg::Vector<bool> BVector;
  typedef linalg::Matrix<double> DMatrix;

  // ??? Introduce half-dynamic MatrixN2, MatrixN3?
  
  class State {
    /* State vector for solar system body */
  public:
    astrometry::CartesianICRS x;   // Position, AU, barycentric, nominally ICRS
    astrometry::CartesianICRS v;   // Velocity, AU per Julian year
    double tdb;                    // TDB for which it applies, Julian yrs since J2000
  };

  class ABG: public Vector6 {
    /* Alpha-beta-gamma basis used for Bernstein-Khushalani */
  public:
    // Indices for components:
    static const int A=0;
    static const int B=1;
    static const int G=2;
    static const int ADOT=3;
    static const int BDOT=4;
    static const int GDOT=5;
    typedef Vector6 Base;
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
		  Vector3& x,
		  Vector3& v) const {
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

    Matrix66 getStateDerivatives() const;
    // Return d(state vector) / d (ABG)
    // State is in reference system, at reference epoch

    DMatrix getXYZ(const DVector& t) const;
    // Calculate inertial XYZ positions at an array of times, return n x 3 matrix

    // Read and write
    std::ostream& write(std::ostream& os, int precision=7) const;
    std::istream& read(std::istream& is);
    static std::ostream& writeHeader(std::ostream& os, int precision=7); // Write column header
  };

  extern std::ostream& operator<<(std::ostream& os, const ABG& abg);
  extern std::istream& operator>>(std::istream& is, ABG& abg);

  class ABGCovariance: public Matrix66 {
  public:
    ABGCovariance(const Matrix66& rhs): Matrix66(rhs) {}
    ABGCovariance() = default;
    std::ostream& write(std::ostream& os, int precision=6) const;
    std::istream& read(std::istream& is);
  };

  class Elements: public Vector6 {
    // Keplerian elements for elliptical orbit
  public:
    static const int A=0;
    static const int E=1;
    static const int I=2;
    static const int LAN=3; // Longitude of ascending node
    static const int AOP=4; // Argument of perihelion
    static const int TOP=5; // Time of perihelion passage
    // Read and write - angles in degrees for I/O; precision is basically degrees of posn.
    std::ostream& write(std::ostream& os, int precision=7) const;
    std::istream& read(std::istream& is);
    static std::ostream& writeHeader(std::ostream& os, int precision=7); // Write column header
  };

  extern std::ostream& operator<<(std::ostream& os, const Elements& el);
  extern std::istream& operator>>(std::istream& is, Elements& el);

  class ElementCovariance: public Matrix66 {
  public:
    ElementCovariance(const Matrix66& rhs): Matrix66(rhs) {}
    ElementCovariance() = default;
    std::ostream& write(std::ostream& os, int precision=6) const;
    std::istream& read(std::istream& is);
  };

  class Observation {
    // General observer-frame information about a detection
  public:
    EIGEN_NEW
    Observation() {}
    astrometry::SphericalICRS radec;
    double tdb;  //  TDB of observation, referenced to J2000
    Matrix22 cov; // RA/Dec error covariance matrix, radians
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
    // isVelocity=true will skip translation of origin
    Vector3 toICRS(const Vector3& x,
		   bool isVelocity=false) const;
    Vector3 fromICRS(const Vector3& x,
		     bool isVelocity=false) const;
    DMatrix toICRS(const DMatrix& x,
		  bool isVelocity=false) const;
    DMatrix fromICRS(const DMatrix& x,
		    bool isVelocity=false) const;

    // Write/read one-line specification.
    std::ostream& write(std::ostream& os) const;
    std::istream& read(std::istream& os);
    static std::ostream& writeHeader(std::ostream& os);
    
  };

} // end namespace orbits

#endif
