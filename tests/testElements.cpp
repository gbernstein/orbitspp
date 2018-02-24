// Testing program for state<->element conversion.
// ephemeris via cspice.

#include "Elements.h"
#include "Ephemeris.h"
#include <iostream>
#include <iomanip>
using namespace std;

using orbits::Elements;

void spinIt(orbits::State& s, double theta) {
  // Rotate state vector about ecliptic z axis CCW by theta
  astrometry::Matrix33 spin(0.);
  spin(0,0) = cos(theta);
  spin(1,1) = cos(theta);
  spin(2,2) = 1.;
  spin(0,1) = -sin(theta);
  spin(1,0) = +sin(theta);

  astrometry::CartesianEcliptic tx(spin * astrometry::CartesianEcliptic(s.x).getVector());
  s.x.convertFrom(tx);
  astrometry::CartesianEcliptic tv(spin * astrometry::CartesianEcliptic(s.v).getVector());
  s.v.convertFrom(tv);
}

void printElements(const orbits::Elements& e) {
  cout << " A: " << e[Elements::A]
       << " E: " << e[Elements::E]
       << " I: " << e[Elements::I]/DEGREE
       << " LAN: " << e[Elements::LAN]/DEGREE
       << " AOP: " << e[Elements::AOP]/DEGREE
       << " TOP: " << e[Elements::TOP]
       << endl;
}

bool testDerivs(const orbits::State& s) {
  bool fail = false;
  auto derivs = orbits::aei_derivs(s);
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
  return fail;
}

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
    cout << "-------------- Analytic elements ---------- " << endl;
    printElements(e);
    
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

    cout << endl;
    cout << "------------ Convert elements back to state vector -----:" << endl;
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

    cout << endl;
    cout << "------------ QB1 heliocentric elements: " << endl;
    s.x[0] = 4.080592008663736E+01;
    s.x[1] = 2.262269771879630E+00;
    s.x[2] = 1.103045097983797E+00;
    s.v[0] = -1.607010770801691E-04/DAY;
    s.v[1] =  2.504643884337053E-03/DAY;
    s.v[2] =  1.201121986118657E-03/DAY;
    s.tdb = eph.utc2tdb("1995-Apr-18 00:00:00") - 61.185610*SECOND;
    auto qb = getElements(s, true);
    printElements(qb);
    {
      /** Compare to these elements from horizons:
	  EC= .06971140563997051  QR= 40.87890604759809   TP= 2448965.9524885667      
	  OM= 359.3985763074872   W=  .7782031838266231   IN= 2.182573356133807       
      **/
      double tolerance = 1e-5;
      peri = 40.87890604759809;
      orbits::Elements truth;
      truth[Elements::E] = 0.06971140563997051;
      truth[Elements::I] = 2.182573356133807 * DEGREE;
      truth[Elements::LAN] = 359.3985763074872 * DEGREE;
      truth[Elements::AOP] = 0.7782031838266231 * DEGREE;
      truth[Elements::TOP] = (2448965.9524885667 - JD2000) * DAY;
      truth[Elements::A] = peri / (1-truth[Elements::E]);
      if (abs(qb[Elements::A]-truth[Elements::A]) > tolerance*peri ||
	  abs(qb[Elements::E]-truth[Elements::E]) > tolerance ||
	  abs(qb[Elements::I]-truth[Elements::I]) > tolerance ||
	  abs(qb[Elements::LAN]-truth[Elements::LAN]) > tolerance ||
	  abs(qb[Elements::AOP]-truth[Elements::AOP]) > tolerance ||
	  abs(qb[Elements::TOP]-truth[Elements::TOP]) > DAY) {
	fail = true;
	cout << "***FAILURE: QB1 elements mismatch" << endl;
	printElements(truth);
      }
    }  
    
    cout << endl;
    cout << "------ Rotate LAN by 10 degrees -----" << endl;
    spinIt(s,10.*DEGREE);
    double oldLAN = qb[Elements::LAN];
    qb = getElements(s, true);
    printElements(qb);
    double diff = qb[Elements::LAN] - oldLAN - 10*DEGREE;
    if (diff > PI) diff -= 2*PI; // Correct any wraps
    if (diff < PI) diff += 2*PI;
    if (abs(diff)>1e-5*DEGREE) {
      fail = true;
      cout << "****FAILURE: incorrent LAN rotation, " << diff/DEGREE << endl;
    }
    
    cout << endl;
    cout << "--------- d(elements) / d(state) in LAN quadrants ----" << endl;
    // Now get the state somewhere away from perihelion and check derivatives
    // against finite differences.
    tdb += 0.5*YEAR;
    s = getState(e, tdb);
    spinIt(s, 10.*DEGREE); // Move away from LAN=180 to avoid bad derivs.

    for (int ispin = 0; ispin < 4; ispin++) {
      auto el = getElements(s);
      printElements(el);
      testDerivs(s);
      if (ispin<3) spinIt(s, 90.*DEGREE);
      cout << endl;
    }

    cout << "--------- Rotate AOP through quadrants ------" << endl;
    // Now explore derivative errors with changing AOP quadrants
    auto el0 = orbits::getElements(s);
    el0[Elements::AOP] = 280*DEGREE;
    s = getState(el0, tdb);
    for (int ispin=0; ispin<4; ispin++) {
      auto el = getElements(s);
      printElements(el);
      testDerivs(s);
      el0[Elements::AOP] += 90*DEGREE;
      s = getState(el0,tdb);
    }      

  } catch (std::runtime_error& e) {
    cerr << e.what() << endl;
    exit(1);
  }
  
  exit(fail ? 1 : 0);
}
/***
For QB1, from Horizons:
Initial IAU76/J2000 heliocentric ecliptic osculating elements (au, days, deg.):
  EPOCH=  2449825.5 ! 1995-Apr-18.00 (TDB)         Residual RMS= .70646        
  Equivalent ICRF heliocentric equatorial cartesian coordinates (au, au/d):
   X= 4.080592008663736E+01  Y= 2.262269771879630E+00  Z= 1.103045097983797E+00
  VX=-1.607010770801691E-04 VY= 2.504643884337053E-03 VZ= 1.201121986118657E-03
Asteroid physical parameters (km, seconds, rotational period in hours):        
   GM= n.a.                RAD= n.a.               ROTPER= n.a.                
   H= 7.1                  G= .150                 B-V= n.a.                   
                           ALBEDO= n.a.            STYP= n.a.                  
*******************************************************************************
JDTDB
   X     Y     Z
   VX    VY    VZ
*******************************************************************************
$$SOE
2449825.500000000 = A.D. 1995-Apr-18 00:00:00.0000 TDB [del_T=     61.185610 s]
 X = 4.080421182418573E+01 Y = 2.269168822159551E+00 Z = 1.106062792338629E+00
 VX=-1.678826496378512E-04 VY= 2.505310726716511E-03 VZ= 1.201609718550474E-03
***/	    
