// Testing program for Ephemeris interface to JPL DE430 (or other)
// ephemeris via cspice.

#include "Ephemeris.h"
#include <iostream>
#include <iomanip>
using namespace std;

int
main(int argc,
     char *argv[])
{
  bool fail = false;  // Flag any failures

  // Test creation of ephemeris with bad filename
  bool caught = false;
  try {
    orbits::Ephemeris junk("nonexistent.file");
  } catch (orbits::SpiceError& e) {
    caught = true;
    cout << "Caught expected exception.  Message:" << endl;
    cout << e.what() << endl;
  }
  if (!caught) {
    cout << "***FAILURE: should have had exception for bad spice kernel" << endl;
    fail = true;
  }
    
  // Put rest of program into try block:
  try {
    orbits::Ephemeris eph;
  
    cout << "-------------" << endl;
    // Test conversion to TDB
    string utc = "2018 Feb 21 06:00:00";
    double tdb = eph.utc2tdb(utc);
    cout << "TDB for " + utc + ": " << std::setprecision(10) << tdb << endl;

    {
      // Check against horizons
      const double tolerance = 0.1 * TIMESEC;
      const double tdbMinusUT = 69.185262*TIMESEC; // From Horizons for this date
      const double answer = 2458170.75;  // Correct UTC JD for this date
      double error = tdb - tdbMinusUT - (answer - JD2000)*DAY;
      if (abs(error) > tolerance) {
	fail = true;
	cout << "***FAILURE: us - Horizons = " << error/TIMESEC << " seconds" << endl;
      }
    } 

    cout << "-------------" << endl;
    // Check conversion of our UT class into TDB
    astrometry::UT ut2(2018,2,21,6,0,0.);
    tdb = eph.utc2tdb(ut2);
    cout << "TDB from UT object: " << std::setprecision(10) << tdb << endl;

    {
      // Check against horizons
      const double tolerance = 0.1 * TIMESEC;
      const double tdbMinusUT = 69.185262*TIMESEC; // From Horizons for this date
      const double answer = 2458170.75;  // Correct UTC JD for this date
      bool thisFail = false;
      double error = tdb - tdbMinusUT - (answer - JD2000)*DAY;
      if (abs(error) > tolerance) thisFail = true;
      if (thisFail) {
	fail = true;
	cout << "***FAILURE: ";
      }
      cout << "us - Horizons = " << error/TIMESEC << " seconds" << endl;

      // Since Horizons won't do vectors with UTC, we'll be setting TDB equal to utc string.
      tdb -= tdbMinusUT;
    }


    cout << "-------------" << endl;
    // Check heliocentric position of Saturn against Horizons.
    auto s = eph.state(orbits::SATURN, tdb);
    cout << "Saturn position: " << s.x << endl;
    {
      // Check against horizons
      double tolerance = 1 * METER;
      astrometry::Vector3 answer;
      answer[0] = 3.180486880720022E-01;
      answer[1] =-9.287696036768155E+00;
      answer[2] =-3.849999488669508E+00;				      
      bool thisFail = false;
      astrometry::Vector3 error = s.x.getVector() - answer;
      for (int i=0; i<3; i++) 
	if (abs(error[i]) > tolerance) 
	  thisFail = true;
      if (thisFail) {
	fail = true;
	cout << "***FAILURE: ";
      }
      error /= METER;
      cout << "us-Horizons = " << std::fixed << std::setprecision(3)
	   << error[0] << " " << error[1] << " " << error[2] << " meters" << endl;

      // Build cache for Saturn
      eph.cachePositions(orbits::SATURN);

      astrometry::Vector3 posn = eph.position(orbits::SATURN, tdb).getVector();
      cout << endl;
      cout << "...and using cache: " << endl;
      error = posn - answer;
      tolerance = 1e4 * METER;
      for (int i=0; i<3; i++) 
	if (abs(error[i]) > tolerance) 
	  thisFail = true;
      if (thisFail) {
	fail = true;
	cout << "***FAILURE: ";
      }
      error /= METER;
      cout << "us-Horizons = " << std::fixed << std::setprecision(3)
	   << error[0] << " " << error[1] << " " << error[2] << " meters" << endl;

      double tdb2 = 16.838;
      astrometry::Vector3 cached = eph.position(orbits::SATURN,tdb2).getVector();
      astrometry::Vector3 spice = eph.position(orbits::SATURN,tdb2,false).getVector();
      cout << endl;
      cout << "...and using cache at " << tdb2 << ": " << endl;
      error = cached-spice;
      for (int i=0; i<3; i++) 
	if (abs(error[i]) > tolerance) 
	  thisFail = true;
      if (thisFail) {
	fail = true;
	cout << "***FAILURE: ";
      }
      error /= METER;
      cout << "cache-spice = " << std::fixed << std::setprecision(3)
	   << error[0] << " " << error[1] << " " << error[2] << " meters" << endl;

      tdb2 = 17.555;
      cached = eph.position(orbits::SATURN,tdb2).getVector();
      spice = eph.position(orbits::SATURN,tdb2,false).getVector();
      cout << endl;
      cout << "...and using cache at " << tdb2 << ": " << endl;
      error = cached-spice;
      for (int i=0; i<3; i++) 
	if (abs(error[i]) > tolerance) 
	  thisFail = true;
      if (thisFail) {
	fail = true;
	cout << "***FAILURE: ";
      }
      error /= METER;
      cout << "cache-spice = " << std::fixed << std::setprecision(3)
	   << error[0] << " " << error[1] << " " << error[2] << " meters" << endl;
      
    }

    cout << "-------------" << endl;
    cout << "Saturn velocity: " << s.v << endl;

    {
      // Check against horizons
      const double tolerance = 1 * METER / TIMESEC;
      astrometry::Vector3 answer;
      answer[0] = 5.269512929058580E-03/DAY;
      answer[1] = 2.298986960400073E-04/DAY;
      answer[2] =-1.319317410084509E-04/DAY;
      bool thisFail = false;
      astrometry::Vector3 error = s.v.getVector() - answer;
      for (int i=0; i<3; i++) 
	if (abs(error[i]) > tolerance) 
	  thisFail = true;
      if (thisFail) {
	fail = true;
	cout << "***FAILURE: ";
      }
      error /= (METER/TIMESEC);
      cout << "us-Horizons = " << std::fixed << std::setprecision(8)
	   << error[0] << " " << error[1] << " " << error[2] << " m/s" << endl;
    }

    cout << "-------------" << endl;
    // Now get observatory position of CTIO
    double lon = 289.194100 * 3.14159265358979 / 180.;
    double lat = -30.169117 *3.141592653588979 / 180.;
    double elev = 2.3888790;

    astrometry::CartesianICRS obs = eph.observatory(lon, lat, elev, tdb);
    cout << "CTIO Position: " << obs << endl;
    {
      // Compare to horizons
      const double tolerance = 1 * METER;
      astrometry::Vector3 answer;
      answer[0] = -8.742296494899074E-01;
      answer[1] = 4.276268427758197E-01;
      answer[2] = 1.852377240773099E-01;
      bool thisFail = false;
      astrometry::Vector3 error = obs.getVector() - answer;
      for (int i=0; i<3; i++) 
	if (abs(error[i]) > tolerance) 
	  thisFail = true;
      if (thisFail) {
	fail = true;
	cout << "***FAILURE: ";
      }
      error /= METER;
      cout << "us-Horizons = " << std::fixed << std::setprecision(3)
	   << error[0] << " " << error[1] << " " << error[2] << " meters" << endl;
    }

    cout << "-------------" << endl;
    // Now repeat using observatory id
    obs = eph.observatory(807, tdb);
    cout << "CTIO Position from 807: " << obs << endl;
    {
      // Compare to horizons
      const double tolerance = 10 * METER;
      astrometry::Vector3 answer;
      answer[0] = -8.742296494899074E-01;
      answer[1] = 4.276268427758197E-01;
      answer[2] = 1.852377240773099E-01;
      bool thisFail = false;
      astrometry::Vector3 error = obs.getVector() - answer;
      for (int i=0; i<3; i++) 
	if (abs(error[i]) > tolerance) 
	  thisFail = true;
      if (thisFail) {
	fail = true;
	cout << "***FAILURE: ";
      }
      error /= METER;
      cout << "us-Horizons = " << std::fixed << std::setprecision(3)
	   << error[0] << " " << error[1] << " " << error[2] << " meters" << endl;
    }
    
  } catch (std::runtime_error& e) {
    cerr << e.what() << endl;
    exit(1);
  }
  
  exit(fail ? 1 : 0);
}

