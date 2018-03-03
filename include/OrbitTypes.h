/** Classes for orbit fitting **/
#ifndef ORBITTYPES_H
#define ORBITTYPES_H

#include "AstronomicalConstants.h"

#include "LinearAlgebra.h"
#include "Astrometry.h"

namespace orbits {

  class State {
  public:
    astrometry::CartesianICRS x;
    astrometry::CartesianICRS v;
    double tdb;
  };

  class ABG: public linalg::SVector<double,6> {
  public:
    const int A=0;
    const int B=1;
    const int G=2;
    const int ADOT=3;
    const int BDOT=4;
    const int GDOT=5;
    ABG(const State& s) {} //???
  };

  class Elements: public linalg::SVector<double,6> {
  public:
    static const int A=0;
    static const int E=1;
    static const int I=2;
    static const int LAN=3; // Longitude of ascending node
    static const int AOP=4; // Argument of perihelion
    static const int TOP=5; // Time of perihelion passage
  };

  class MPCObservation {
  public:
    MPCObservation() {}
    MPCObservation(const string& line); // Parse MPC style input string
    astrometry::SphericalICRS radec;
    double mjd;  // UTC MJD of observation
    double sigma; // Position uncertainty per component
    int obscode; // MPC observatory code
  };

  typedef linalg::SMatrix<double,6,6> ABGCovar;
  
  class Frame: public astrometry::ReferenceFrame {
  public:
    double tdb0;  // Time coordinate origin
  };

} // end namespace orbits

#endif
