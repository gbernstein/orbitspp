// Declarations of functions to produce outputs in the old ASCII formats from C code
#ifndef LEGACY_H
#define LEGACY_H

#include "OrbitTypes.h"
#include "Elements.h"

namespace orbits {

  extern void
  writeOldAEI(string filename, const Elements& elements,
	      const ElementCovariance& cov,
	      double tdb0,
	      Ephemeris& eph);

  extern void
  writeOldABG(string filename, const ABG& abg,
	      const ABGCovariance& abgcov,
	      const Frame& frame,
	      Ephemeris& eph);
}

#endif  //LEGACY_H
