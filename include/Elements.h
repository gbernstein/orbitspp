// Declarations of functions for conversion between state vector and orbital elements
#ifndef ELEMENTS_H
#define ELEMENTS_H

#include "OrbitTypes.h"
#include "Ephemeris.h"

namespace orbits {

  /**  Return orbital elements in J2000 ecliptic
       given an ICRS state vector.  The state vector is **always**
       barycentric.  The heliocentric flag tells whether elements
       that are returned are osculating heliocentric (solar mass
       is used and solar motion is removed from state) or 
       barycentric (in which case full solar system mass is used).
       If heliocentric is chosen, we need an Ephemeris.
  **/
  extern
  Elements getElements(const State& s,
		       bool heliocentric=false, const Ephemeris* ephem=nullptr);

  /** Get state from orbital elements
      tdb is dynamical barycentric time, in years past J2000 TDB.
      Again elements are either barycentric, total SS mass, or
      they are osculating heliocentric with solar mass.  The
      State returned is **always** barycentric ICRS. 
      If heliocentric is chosen, we need an Ephemeris.
  **/
  extern 
  State getState(const Elements& el, double tdb,
		 bool heliocentric=false, const Ephemeris* ephem=nullptr);

  // Get derivative matrix of State->Elements - only available as barycentric
  extern
  Matrix66
  getElementDerivatives(const orbits::State& xv,
			bool heliocentric=false, const Ephemeris* ephem=nullptr);
  
  extern
  Elements
  heliocentric2barycentric(const Elements& el, double tdb, const Ephemeris& ephem);

  extern
  Elements
  barycentric2heliocentric(const Elements& el, double tdb, const Ephemeris& ephem);
  
} // end namespace orbits

#endif // ELEMENTS_H
