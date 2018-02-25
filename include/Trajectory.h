// Class representing target body trajectory, including
// integration of gravitational dynamics
#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include "Astrometry.h"
#include "OrbitTypes.h"
#include "Ephemeris.h"
#include "LinearAlgebra.h"

namespace orbits {

  enum Gravity {INERTIAL, BARY, GIANTS};
  // INERTIAL = no gravity (inertial)
  // BARY = gravity from point mass at solar system barycenter with mass
  //       of all 8 planets
  // GIANTS = gravity from giant planets considered distinctly, terrestrials
  //       are placed with the Sun.

  class Trajectory {
  public:
    Trajectory(const Ephemeris& ephem_,
	       const State& s0,
	       double dt_=20*DAY,
	       Gravity grav_ = GIANTS);
    // dt will be time step for integrators.  All positions are in
    // ICRS barycentric coordinates, all times are TDB, units are AU and
    // years.

    linalg::Matrix<double> position(const linalg::Vector<double>& tdb) const;
    // Return 3 x N matrix of object positions at TDB's given in input array.
    // Will assume that times are in ascending order in array.

    astrometry::CartesianICRS position(double tdb) const;
    // Return position at a single time.

    EIGEN_NEW
    
  private:
    const Ephemeris& ephem;
    astrometry::Vector3 x0;  // Initial position
    astrometry::Vector3 v0;  // Initial velocity
    double tdb0;
    double dt;
    Gravity grav;
    // Caches for positions and velocities at time steps
    // after and before initial state.
    // And accelerations, to make quadratic interp easy.
    mutable std::vector<astrometry::Vector3> xfwd;
    mutable std::vector<astrometry::Vector3> xbwd;
    mutable std::vector<astrometry::Vector3> vfwd;
    mutable std::vector<astrometry::Vector3> vbwd;
    mutable std::vector<astrometry::Vector3> afwd;
    mutable std::vector<astrometry::Vector3> abwd;
    mutable astrometry::Vector3 vnext_fwd;
    mutable astrometry::Vector3 vnext_bwd;
    
    // Acceleration calculator - return a*dt
    astrometry::Vector3 deltaV(const astrometry::Vector3& x, double tdb) const;
    // Extend integration to include this time range
    void span(double tdbMin, double tdbMax) const;
  };

} // end namespace orbits
#endif // TRAJECTORY_H
