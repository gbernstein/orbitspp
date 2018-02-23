// Testing program for state<->element conversion.
// ephemeris via cspice.

#include "Elements.h"
#include "Ephemeris.h"
#include <iostream>
#include <iomanip>
using namespace std;

using orbits::Elements;

int
main(int argc,
     char *argv[])
{
  bool fail = false;  // Flag any failures

  try {
    orbits::Ephemeris eph;
    double tdb = eph.utc2tdb("2018 Jan 1 00:00:00");
  
    // Make an elliptical orbit in ICRS plane with peri at
    // x=0, so 90 deg from ascending node in ecliptic coords
    double peri = 1.5;
    double ee = 0.21; // ellipticity to produce
    double aa = peri / (1-ee); // semi-major we should get
    orbits::State s;
    s.x[0] = s.x[2] = 0.;
    s.x[1] = peri;
    s.v[0] = -sqrt(SolarSystemGM * (1+ee) / peri);
    s.v[1] = s.v[2] = 0.;
    s.tdb = tdb;

    orbits::Elements e = orbits::getElements(s);
    cout << "Reconstructed elements: " << endl;
    cout << "  A, E, I: " << e[Elements::A] << " " << e[Elements::E]
	 << " " << e[Elements::I]/DEGREE << endl;
    cout << "  LAN, AOP: " << e[Elements::LAN]/DEGREE << " " << e[Elements::AOP]/DEGREE
	 << " TOP " << e[Elements::TOP] << endl;
    
    // Check for failures
    const double tolerance = 1e-8;
    if ( abs(e[Elements::A]-aa) > tolerance*peri ||
	 abs(e[Elements::E]-ee) > tolerance ||
	 abs(e[Elements::I] -EclipticInclination ) > tolerance ||
	 abs(e[Elements::LAN]-PI) > tolerance ||
	 abs(e[Elements::AOP] - 1.5*PI) > tolerance ||
	 abs(e[Elements::TOP] - tdb) > 0.1*SECOND) {
      cout << "***FAILURE: orbit not as expected" << endl;
      fail = true;
    }

    // Convert back to state:
    orbits::State s2 = getState(e, tdb);
    s2.x -= s.x;
    s2.v -= s.v;
    cout << "x diff " << s2.x << " v " << s2.v << endl;
    {
      double tolerance = METER;
      auto vv = s2.x.getVector();
      if (vv.dot(vv) > tolerance*tolerance) {
	cout << "***FAILURE: position not recovered" << endl;
	fail = true;
      }
      tolerance = METER / YEAR;
      vv = s2.v.getVector();
      if (vv.dot(vv) > tolerance*tolerance) {
	cout << "***FAILURE: position not recovered" << endl;
	fail = true;
      }
    }
    
    cout << "--------- d(elements) / d(state) ---------" << endl;
    // Now get the state somewhere away from perihelion and check derivatives
    // against finite differences.
    tdb += 0.5*YEAR;
    s = getState(e, tdb);
    s.v[2] += 0.1;  // Need to move away from LAN=180 to avoid bad derivs
    auto derivs = aei_derivs(s);
    cout << derivs << endl;

    for (int i=0; i<6; i++) {
      const double dx=0.01;
      orbits::State splus(s);
      orbits::State sminus(s);
      if (i<3) {
	splus.x[i] += dx/2.;
	sminus.x[i] -= dx/2.;
      } else {
	splus.v[i-3] += dx/2.;
	sminus.v[i-3] -= dx/2.;
      }
	
      auto eplus = orbits::getElements(splus);
      auto eminus = orbits::getElements(sminus);
      for (int j=0; j<6; j++) {
	if ( abs( (eplus[j]-eminus[j])/(dx*derivs(j,i))-1.) > 0.001) {
	  cout << "****FAILURE: deriv("<<j<<","<<i<<") mismatch.  Numerical: " 
	       << (eplus[j]-eminus[j])/dx
	       << " Analytic: " << derivs(j,i)
	       << " diff " << (eplus[j]-eminus[j])/(dx*derivs(j,i))-1
	       << endl;
	  fail = true;
	}
      }
    }

  } catch (std::runtime_error& e) {
    cerr << e.what() << endl;
    exit(1);
  }
  
  exit(fail ? 1 : 0);
}

