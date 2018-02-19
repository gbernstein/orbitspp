/** Classes for orbit fitting **/
#ifndef ORBITTYPES_H
#define ORBITTYPES_H

#include "AstronomicalConstants.h"

#include "LinearAlgebra.h"
#include "Astrometry.h"

namespace orbits {

  class State {
  public:
    astrometry::CartesianEcliptic x;
    astrometry::CartesianEcliptic v;
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
    const int A=0;
    const int E=1;
    const int I=2;
    const int LAN=3; // Longitude of ascending node
    const int AOP=4; // Argument of perihelion
    const int TOP=5; // Time of perihelion passage
  };

  class Observation {
  public:
    astrometry::Vector2 theta;    // Angular coordinates, in chosen projection
    astrometry::Matrix22 invcov;  // Inverse covariance of measurement error
    double tdb;  // TDB of observation
    astrometry::Vector3 xe;  // Position of observatory in chosen frame
    int obscode; // MPC observatory code
  };

  class ABGCovar: public linalg::SMatrix<double,6,6> {
  };

  class Frame: public astrometry::ReferenceFrame {
  public:
    double tdb0;  // Time coordinate origin
  };

} // end namespace orbits

#endif
