// Class representing target body trajectory, including
// integration of gravitational dynamics
#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include "Astrometry.h"
#include "OrbitTypes.h"
#include "Ephemeris.h"
#include "LinearAlgebra.h"
#include "Elements.h"

// Need to following to make a vector containing structures with fixed-size
// Eigen arrays
#ifdef USE_EIGEN
#include <Eigen/StdVector>
#endif

namespace orbits {

  enum Gravity {INERTIAL, BARY, GIANTS, WRONG};
  // INERTIAL = no gravity (inertial)
  // BARY = gravity from point mass at solar system barycenter with mass
  //       of all 8 planets
  // GIANTS = gravity from giant planets considered distinctly, terrestrials
  //       are placed with the Sun.
  // WRONG = trajectory is an ellipse relative to the Sun (not barycenter).
  //       This is not a valid model but is what PyEphem does.  The
  //       code for this will be slow, redoing the ellipse calculations
  //       every time.

  class Trajectory {
  public:
    Trajectory(const Ephemeris& ephem_,
	       const State& s0,
	       Gravity grav_ = GIANTS,
	       double dt_=20*DAY);
    // dt will be time step for integrators.  All positions are in
    // ICRS barycentric coordinates, all times are TDB, units are AU and
    // years.

    DMatrix position(const DVector& tdb,
		     DMatrix* velocity=nullptr) const;
    // Return Nx3 matrix of object positions at TDB's given in input array.
    // If a pointer to velocity array is given, it is resized to Nx3 and filled with v.

    astrometry::CartesianICRS position(double tdb,
				       astrometry::CartesianICRS* velocity=nullptr) const;
    // Return position at a single time. Velocity is put into 2nd arg if given.

    // Return observed astrometric position from observer position/time.
    astrometry::SphericalICRS observe(double tdbObserve,
				      const astrometry::CartesianICRS& observer);
    // Return observed astrometric position in bulk.
    // Input and output matrices are Nx3.  Output is ICRS direction cosines.
    DMatrix observe(const DVector& tdbObserve,
		    const DMatrix& observer);
    EIGEN_NEW
    
  private:
    const Ephemeris& ephem;
    Vector3 x0;  // Initial position
    Vector3 v0;  // Initial velocity
    double tdb0;
    double dt;
    Gravity grav;
    Elements el;            //  For the WRONG model.
    // Need to use special form for vectors containing static-size
    // Eigen structures.
#ifdef USE_EIGEN
    typedef std::vector<Vector3, Eigen::aligned_allocator<Vector3> > v3vector;
#else
    typedef std::vector<Vector3> v3vector;
#endif
    // Caches for positions and velocities at time steps
    // after and before initial state.
    // And accelerations, to make quadratic interp easy.
    mutable v3vector xfwd;
    mutable v3vector xbwd;
    mutable v3vector vfwd;
    mutable v3vector vbwd;
    mutable v3vector afwd;
    mutable v3vector abwd;
    mutable Vector3 vnext_fwd;
    mutable Vector3 vnext_bwd;
    
    // Acceleration calculator - return a*dt
    Vector3 deltaV(const astrometry::Vector3& x, double tdb) const;
    // Extend integration to include this time
    void span(double tdb) const;
  };

} // end namespace orbits
#endif // TRAJECTORY_H
