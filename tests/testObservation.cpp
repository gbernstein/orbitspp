// Test whether RA, Dec for 1992 QB1 agree with Horizons when
// given Horizons heliocentric elements.

#include "OrbitTypes.h"
#include "Ephemeris.h"
#include "Trajectory.h"
#include <iostream>
#include <iomanip>
#include "StringStuff.h"

using namespace std;

using namespace orbits;

string horizons_text =
  " 2000-Jan-01 00:00 2451544.500000000  7.9765082   3.8891337  88.9980 \n"
  " 2000-Jul-01 00:00 2451726.500000000 11.2049037   5.2895046  87.0941 \n"
  " 2001-Jan-01 00:00 2451910.500000000  9.2720952   4.5000016  89.6655 \n"
  " 2001-Jul-01 00:00 2452091.500000000 12.4985461   5.8908537  85.4304 \n"
  " 2002-Jan-01 00:00 2452275.500000000 10.5644964   5.1064611  91.3517 \n"
  " 2002-Jul-01 00:00 2452456.500000000 13.7915551   6.4876986  83.7664 \n"
  " 2003-Jan-01 00:00 2452640.500000000 11.8578625   5.7097622  93.0340 \n"
  " 2003-Jul-01 00:00 2452821.500000000 15.0841021   7.0797024  82.1147 \n"
  " 2004-Jan-01 00:00 2453005.500000000 13.1527667   6.3096961  94.7138 \n"
  " 2004-Jul-01 00:00 2453187.500000000 16.3833337   7.6697793  81.4047 \n"
  " 2005-Jan-01 00:00 2453371.500000000 14.4510949   6.9064467  95.3747 \n"
  " 2005-Jul-01 00:00 2453552.500000000 17.6770379   8.2516721  79.7509 \n"
  " 2006-Jan-01 00:00 2453736.500000000 15.7503089   7.4987881  97.0528 \n"
  " 2006-Jul-01 00:00 2453917.500000000 18.9722248   8.8282827  78.1122 \n"
  " 2007-Jan-01 00:00 2454101.500000000 17.0534560   8.0874505  98.7336 \n"
  " 2007-Jul-01 00:00 2454282.500000000 20.2696332   9.3995239  76.4647 \n"
  " 2008-Jan-01 00:00 2454466.500000000 18.3611048   8.6722118 100.4148 \n"
  " 2008-Jul-01 00:00 2454648.500000000 21.5792967   9.9694859  75.7605 \n"
  " 2009-Jan-01 00:00 2454832.500000000 19.6730518   9.2523543 101.0761 \n"
  " 2009-Jul-01 00:00 2455013.500000000 22.8838928  10.5299583  74.1198 \n"
  " 2010-Jan-01 00:00 2455197.500000000 20.9907797   9.8285436 102.7604 \n"
  " 2010-Jul-01 00:00 2455378.500000000 24.1921554  11.0845795  72.4745 \n"
  " 2011-Jan-01 00:00 2455562.500000000 22.3137822  10.3999484 104.4441 \n"
  " 2011-Jul-01 00:00 2455743.500000000 25.5038063  11.6328703  70.8373 \n"
  " 2012-Jan-01 00:00 2455927.500000000 23.6415602  10.9659754 106.1239 \n"
  " 2012-Jul-01 00:00 2456109.500000000 26.8289128  12.1787961  70.1374 \n"
  " 2013-Jan-01 00:00 2456293.500000000 24.9701015  11.5243355 106.7823 \n"
  " 2013-Jul-01 00:00 2456474.500000000 28.1450808  12.7125175  68.4962 \n"
  " 2014-Jan-01 00:00 2456658.500000000 26.3042133  12.0771263 108.4626 \n"
  " 2014-Jul-01 00:00 2456839.500000000 29.4625015  13.2380765  66.8733 \n"
  " 2015-Jan-01 00:00 2457023.500000000 27.6414326  12.6226467 110.1387 \n"
  " 2015-Jul-01 00:00 2457204.500000000 30.7808517  13.7550403  65.2460 \n"
  " 2016-Jan-01 00:00 2457388.500000000 28.9817593  13.1605210 111.8110 \n"
  " 2016-Jul-01 00:00 2457570.500000000 32.1140496  14.2684925  64.5602 \n"
  " 2017-Jan-01 00:00 2457754.500000000 30.3203024  13.6882189 112.4544 \n"
  " 2017-Jul-01 00:00 2457935.500000000 33.4363442  14.7679675  62.9479 \n"
  " 2018-Jan-01 00:00 2458119.500000000 31.6679794  14.2100742 114.1213 \n";

int
main(int argc,
     char *argv[])
{
  bool fail = false;  // Flag any failures

  try {
    orbits::Ephemeris eph;

    Elements qb;
    /** These are elements from Horizons:
	EPOCH=  2449825.5 ! 1995-Apr-18.00 (TDB)         Residual RMS= .70646  
        EC= .06971140563997051  QR= 40.87890604759809   TP= 2448965.9524885667      
	OM= 359.3985763074872   W=  .7782031838266231   IN= 2.182573356133807       
    **/
    double peri = 40.87890604759809;
    double epoch = eph.jd2tdb(2449825.5);
    qb[Elements::E] = 0.06971140563997051;
    qb[Elements::I] = 2.182573356133807 * DEGREE;
    qb[Elements::LAN] = 359.3985763074872 * DEGREE;
    qb[Elements::AOP] = 0.7782031838266231 * DEGREE;
    qb[Elements::TOP] = (2448965.9524885667 - JD2000) * DAY;
    qb[Elements::A] = peri / (1-qb[Elements::E]);

    // Convert heliocentric elements into an initial state of orbit
    State s0 = getState(qb, epoch, true, &eph);

    // Calculate the orbit:
    Trajectory orb = Trajectory(eph, s0, Gravity::GIANTS, 2*DAY);

    // Now check against horizons:
    std::istringstream iss(horizons_text);
    string line;

    double tolerance = 0.03; // arcsecond accuracy demanded
    cout << "# Comparison of 1992 QB1 predicted position to JPL Horizons" << endl;
    cout << "# Date      RA_predict      Dec_predict  RA_err  DEC_err" << endl;
    cout << "#                                          (arcseconds) " << endl;

    while (stringstuff::getlineNoComment(iss,line)) {
      auto fields = stringstuff::split(line);
      auto fptr = fields.begin();
      string date = (*fptr++);
      *fptr++;  // Skip time of day
      double jd = atof((fptr++)->c_str());
      double ra = atof((fptr++)->c_str()) * DEGREE;
      double dec = atof((fptr++)->c_str()) * DEGREE;

      double tdb = eph.jd2tdb(jd);

      // All observations are for CTIO:
      auto observer = eph.observatory(807, tdb);

      // Get expected position:
      astrometry::SphericalICRS predict = orb.observe(tdb,observer);
      double ra_p, dec_p;
      predict.getLonLat(ra_p, dec_p);
      double ra_err = (ra_p - ra) / ARCSEC / cos(dec_p);
      double dec_err = (dec_p - dec) / ARCSEC;
      cout << date << " " << predict << std::fixed << std::setprecision(4)
	   << " " << ra_err << " " << dec_err << endl;
      if ( abs(ra_err)>tolerance || abs(dec_err)>tolerance)
	fail = true;
    }

    if (fail)
      cout << "****FAILURE: Difference from Horizons too large****" << endl;

  } catch (std::runtime_error& e) {
    cerr << e.what() << endl;
    exit(1);
  }
  
  exit(fail ? 1 : 0);
}
