// Declarations of functions for conversion between state vector and orbital elements
#ifndef ELEMENTS_H
#define ELEMENTS_H

#include "OrbitTypes.h"

namespace orbits {

  // Return orbital elements in J2000 ecliptic
  // given an ICRS state vector.
  // heliocentric flag tells whether state is WRT Sun (in which
  // case solar mass is used for elements) or WRT barycenter
  // (in which case full solar system mass is used)
  extern
  Elements getElements(const State& s, bool heliocentric=false);

  // Get state from orbital elements
  // tdb is dynamical barycentric time, in years past J2000 TDB.
  extern 
  State getState(const Elements& el, double tdb, bool heliocentric=false);

  // Get derivative matrix of State->Elements
  extern
  Matrix66
  getElementDerivatives( const orbits::State& xv, bool heliocentric=false);
  
} // end namespace orbits

#endif // ELEMENTS_H
