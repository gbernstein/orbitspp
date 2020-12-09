// Class representing target body trajectory, including
// integration of gravitational dynamics
#ifndef TRAJECTORY_H
#define TRAJECTORY_H

#include "Astrometry.h"
#include "OrbitTypes.h"
#include "Ephemeris.h"
#include "LinearAlgebra.h"
#include "Elements.h"
#include "SharedLUT.h"

// Need to following to make a vector containing structures with fixed-size
// Eigen arrays
#ifdef USE_EIGEN
#include <Eigen/StdVector>
#endif

// If we are multithreaded, we will need to guard our LUTs.
#ifdef _OPENMP
#include <omp.h>
#endif

namespace orbits {

  enum Gravity {INERTIAL, BARY, GIANTS};
  // INERTIAL = no gravity (inertial)
  // BARY = gravity from point mass at solar system barycenter with mass
  //       of all 8 planets
  // GIANTS = gravity from giant planets considered distinctly, terrestrials
  //       are placed with the Sun.

  struct Circumstances {
    // Structure to hold object distance/phase information
  public:
    double obsDistance;     // Distance from Observer, in AU
    double sunDistance;     // Distance from Sun, in AU
    double phaseAngle;      // Angle between Observer-Target and Sun-Target vectors
    double phasePA;         // position angle (N thru E) of bright hemisphere
  };

  class Trajectory {
  public:
    Trajectory(const Ephemeris& ephem_,
	       const State& s0,
	       Gravity grav_ = GIANTS,
	       double dt_=20*DAY);
    // dt will be time step for integrators.  All positions are in
    // ICRS barycentric coordinates, all times are TDB, units are AU and
    // years.

    // Copy constructor won't copy the LUT
    Trajectory(const Trajectory& rhs); 

    DMatrix position(const DVector& tdb,
		     DMatrix* velocity=nullptr) const;
    // Return Nx3 matrix of object positions at TDB's given in input array.
    // If a pointer to velocity array is given, it is resized to Nx3 and filled with v.

    astrometry::CartesianICRS position(double tdb,
				       astrometry::CartesianICRS* velocity=nullptr) const;
    // Return position at a single time. Velocity is put into 2nd arg if given.

    // Return observed astrometric position from observer position/time.
    astrometry::SphericalICRS observe(double tdbObserve,
				      const astrometry::CartesianICRS& observer) const;
    // Return observed astrometric position in bulk.
    // Input and output matrices are Nx3.  Output is ICRS direction cosines.
    DMatrix observe(const DVector& tdbObserve,
		    const DMatrix& observer) const;

    // Return observational phase/circumstances for observation at given time/posn
    Circumstances getCircumstances(double tdb,
				   const Vector3& observer) const;

    EIGEN_NEW
    
  private:
    const Ephemeris& ephem;
    Vector3 x0;  // Initial position
    Vector3 v0;  // Initial velocity
    double tdb0;
    double dt;
    Gravity grav;

    // Caches for positions, velocities, and accelerations at time steps
    // after and before initial state.
    // Store each direction as a "SharedLUT" which will manage the object
    // to permit multiple-thread reading and limit appending to a single
    // thread at a time.
    mutable SharedLUT<Matrix33,Eigen::aligned_allocator<Matrix33>> xvaLUT;
    
    // Acceleration calculator - return a*dt
    Vector3 deltaV(const astrometry::Vector3& x, double tdb) const;
    // Extend integration to include this time
    void span(double tdb) const;

  };


} // end namespace orbits
#endif // TRAJECTORY_H
