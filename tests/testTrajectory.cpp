// Testing orbit dynamics program against Horizons

#include "Trajectory.h"
#include <iostream>
#include <iomanip>
using namespace std;

int
main(int argc,
     char *argv[])
{
  bool fail = false;  // Flag any failures

  // Put rest of program into try block:
  try {
    orbits::Ephemeris eph;
  
    cout << "-------- Integrate QB1 from 2018 Apr ------" << endl;

    orbits::State s;
    s.tdb = (2458226.750000000 - JD2000) * DAY;
    s.x[0] = 3.341187430342263E+01;
    s.x[1] = 2.193374689734561E+01;
    s.x[2] = 1.053504201793683E+01;
    s.v[0] = -1.539912581639087E-03/DAY;
    s.v[1] = 2.060058504418397E-03/DAY;
    s.v[2] = 9.874446157718628E-04/DAY;

    // Create Trajectory
    orbits::Trajectory qb(eph, s, orbits::GIANTS, 1.18*DAY);
    
    // single future point
    auto x = qb.position((2459322.75-JD2000)*DAY);
    cout << x << endl;

    // Bunch of past/future points (from Horizons)
    vector<double> mjd;
    vector<double> xx, yy, zz;
    mjd.push_back(2452382.7500);
    xx.push_back( 3.980012947065434E+01);
    yy.push_back( 8.613596052420819E+00);
    zz.push_back( 4.148756167814337E+00);
    mjd.push_back(2452535.7500);
    xx.push_back( 3.970398283577915E+01);
    yy.push_back( 8.987228593512992E+00);
    zz.push_back( 4.327928663190391E+00);
    mjd.push_back(2452688.7500);
    xx.push_back( 3.960382249943922E+01);
    yy.push_back( 9.359952500963606E+00);
    zz.push_back( 4.506663593134572E+00);
    mjd.push_back(2452838.7500);
    xx.push_back( 3.950174104458336E+01);
    yy.push_back( 9.724450016057530E+00);
    zz.push_back( 4.681451876758996E+00);
    mjd.push_back(2452991.7500);
    xx.push_back( 3.939366726993092E+01);
    yy.push_back( 1.009526447518327E+01);
    zz.push_back( 4.859267512704392E+00);
    mjd.push_back(2453143.7500);
    xx.push_back( 3.928236157007868E+01);
    yy.push_back( 1.046264595732711E+01);
    zz.push_back( 5.035435112596856E+00);
    mjd.push_back(2453296.7500);
    xx.push_back( 3.916637251253273E+01);
    yy.push_back( 1.083139179062372E+01);
    zz.push_back( 5.212255101920475E+00);
    mjd.push_back(2453447.7500);
    xx.push_back( 3.904802630611830E+01);
    yy.push_back( 1.119424623569682E+01);
    zz.push_back( 5.386248246156004E+00);
    mjd.push_back(2453600.7500);
    xx.push_back( 3.892420192894588E+01);
    yy.push_back( 1.156078523899639E+01);
    zz.push_back( 5.562006328310526E+00);
    mjd.push_back(2453753.7500);
    xx.push_back( 3.879645569136225E+01);
    yy.push_back( 1.192615925453758E+01);
    zz.push_back( 5.737203919505711E+00);
    mjd.push_back(2453904.7500);
    xx.push_back( 3.866654875564041E+01);
    yy.push_back( 1.228557963082659E+01);
    zz.push_back( 5.909544888222892E+00);
    mjd.push_back(2454057.7500);
    xx.push_back( 3.853105524741635E+01);
    yy.push_back( 1.264853189094759E+01);
    zz.push_back( 6.083577514215363E+00);
    mjd.push_back(2454208.7500);
    xx.push_back( 3.839353302905756E+01);
    yy.push_back( 1.300549234140988E+01);
    zz.push_back( 6.254735265695357E+00);
    mjd.push_back(2454361.7500);
    xx.push_back( 3.825035539923472E+01);
    yy.push_back( 1.336588171572835E+01);
    zz.push_back( 6.427535259561065E+00);
    mjd.push_back(2454514.7500);
    xx.push_back( 3.810333544302270E+01);
    yy.push_back( 1.372492830292861E+01);
    zz.push_back( 6.599689513977605E+00);
    mjd.push_back(2454665.7500);
    xx.push_back( 3.795448691863407E+01);
    yy.push_back( 1.407793056425187E+01);
    zz.push_back( 6.768943792630128E+00);
    mjd.push_back(2454818.7500);
    xx.push_back( 3.779988438014895E+01);
    yy.push_back( 1.443420496231127E+01);
    zz.push_back( 6.939765069178111E+00);
    mjd.push_back(2454969.7500);
    xx.push_back( 3.764358738085761E+01);
    yy.push_back( 1.478440334432896E+01);
    zz.push_back( 7.107671228376171E+00);
    mjd.push_back(2455122.7500);
    xx.push_back( 3.748147397111526E+01);
    yy.push_back( 1.513776831290464E+01);
    zz.push_back( 7.277093730186759E+00);
    mjd.push_back(2455273.7500);
    xx.push_back( 3.731780084494582E+01);
    yy.push_back( 1.548502820594003E+01);
    zz.push_back( 7.443587232029667E+00);
    mjd.push_back(2455426.7500);
    xx.push_back( 3.714825144389821E+01);
    yy.push_back( 1.583534824461314E+01);
    zz.push_back( 7.611545990992531E+00);
    mjd.push_back(2455579.7500);
    xx.push_back( 3.697498890752036E+01);
    yy.push_back( 1.618408514718865E+01);
    zz.push_back( 7.778743774213709E+00);
    mjd.push_back(2455730.7500);
    xx.push_back( 3.680037029432700E+01);
    yy.push_back( 1.652667839702075E+01);
    zz.push_back( 7.942994131657692E+00);
    mjd.push_back(2455883.7500);
    xx.push_back( 3.661979032677896E+01);
    yy.push_back( 1.687217037360951E+01);
    zz.push_back( 8.108632272375296E+00);
    mjd.push_back(2456035.7500);
    xx.push_back( 3.643677471379729E+01);
    yy.push_back( 1.721373816077393E+01);
    zz.push_back( 8.272387086739061E+00);
    mjd.push_back(2456188.7500);
    xx.push_back( 3.624893670634914E+01);
    yy.push_back( 1.755584347807266E+01);
    zz.push_back( 8.436397618003094E+00);
    mjd.push_back(2456341.7500);
    xx.push_back( 3.605749033291544E+01);
    yy.push_back( 1.789620111716038E+01);
    zz.push_back( 8.599568284743496E+00);
    mjd.push_back(2456491.7500);
    xx.push_back( 3.586631637786746E+01);
    yy.push_back( 1.822815734391577E+01);
    zz.push_back( 8.758709285693811E+00);
    mjd.push_back(2456644.7500);
    xx.push_back( 3.566779019016800E+01);
    yy.push_back( 1.856495890878862E+01);
    zz.push_back( 8.920171163735869E+00);
    mjd.push_back(2456795.7500);
    xx.push_back( 3.546838734307538E+01);
    yy.push_back( 1.889555085069176E+01);
    zz.push_back( 9.078654180075350E+00);
    mjd.push_back(2456948.7500);
    xx.push_back( 3.526284878528223E+01);
    yy.push_back( 1.922865940876595E+01);
    zz.push_back( 9.238341623615826E+00);
    mjd.push_back(2457099.7500);
    xx.push_back( 3.505657115280673E+01);
    yy.push_back( 1.955554550441974E+01);
    zz.push_back( 9.395044107976117E+00);
    mjd.push_back(2457252.7500);
    xx.push_back( 3.484411396538270E+01);
    yy.push_back( 1.988483775374222E+01);
    zz.push_back( 9.552898015910751E+00);
    mjd.push_back(2457405.7500);
    xx.push_back( 3.462821102408589E+01);
    yy.push_back( 2.021216335412904E+01);
    zz.push_back( 9.709807103136431E+00);
    mjd.push_back(2457557.7500);
    xx.push_back( 3.441033164909292E+01);
    yy.push_back( 2.053537215623550E+01);
    zz.push_back( 9.864740679673003E+00);
    mjd.push_back(2457710.7500);
    xx.push_back( 3.418763399251348E+01);
    yy.push_back( 2.085868700736372E+01);
    zz.push_back( 1.001972301247134E+01);
    mjd.push_back(2457861.7500);
    xx.push_back( 3.396454255785577E+01);
    yy.push_back( 2.117575915889163E+01);
    zz.push_back( 1.017171082289694E+01);
    mjd.push_back(2458014.7500);
    xx.push_back( 3.373517316879448E+01);
    yy.push_back( 2.149495874946377E+01);
    zz.push_back( 1.032471631830249E+01);
    mjd.push_back(2458167.7500);
    xx.push_back( 3.350248485646207E+01);
    yy.push_back( 2.181204351357604E+01);
    zz.push_back( 1.047670596582468E+01);
    mjd.push_back(2458317.7500);
    xx.push_back( 3.327116266512767E+01);
    yy.push_back( 2.212083010860435E+01);
    zz.push_back( 1.062471594044462E+01);
    mjd.push_back(2458470.7500);
    xx.push_back( 3.303198010005568E+01);
    yy.push_back( 2.243364196570766E+01);
    zz.push_back( 1.077465320658009E+01);
    mjd.push_back(2458621.7500);
    xx.push_back( 3.279274825409331E+01);
    yy.push_back( 2.274020808554873E+01);
    zz.push_back( 1.092159464131932E+01);
    mjd.push_back(2458774.7500);
    xx.push_back( 3.254715704453595E+01);
    yy.push_back( 2.304862186126397E+01);
    zz.push_back( 1.106941952491420E+01);
    mjd.push_back(2458926.7500);
    xx.push_back( 3.230001805146345E+01);
    yy.push_back( 2.335278716492671E+01);
    zz.push_back( 1.121520592326293E+01);
    mjd.push_back(2459079.7500);
    xx.push_back( 3.204810739926869E+01);
    yy.push_back( 2.365667907943359E+01);
    zz.push_back( 1.136085909753957E+01);
    mjd.push_back(2459232.7500);
    xx.push_back( 3.179306902864494E+01);
    yy.push_back( 2.395826219011573E+01);
    zz.push_back( 1.150540347107284E+01);
    mjd.push_back(2458114.5);  // 2017-Dec-27 00:00:00.0000 TDB
    xx.push_back( 3.358384441110039E+01);
    yy.push_back( 2.170192732308413E+01);
    zz.push_back( 1.042392372599419E+01);
    astrometry::CartesianICRS vv(-1.524127150690516E-03/DAY,  // Horizons velocity for this date
				 2.070340276060543E-03/DAY,
				 9.923832708824847E-04/DAY);
    cout << endl;
    cout << "---Difference against Horizons (km) ---" << endl;
    cout << " DT (yrs)     x         y       z " << endl;

    // This integration with dt = 1 day should be good to about 1 km
    // over the 16-yr integration above.
    const double tolerance = 3e3 * METER; 
    linalg::Vector<double> future(mjd.size());
    for (int i=0; i<future.size(); i++)
      future[i] = (mjd[i] - JD2000)*DAY;

    linalg::DMatrix velocity;
    auto m = qb.position(future, &velocity);
    for (int i=0; i<future.size(); i++) {
      cout << std::fixed << std::setprecision(3) << std::setw(8) << future[i] - s.tdb
	   << " " << std::setw(8) << std::setprecision(2) <<  (m(0,i) - xx[i]) / (1000.*METER)
	   << " " << std::setw(8) << std::setprecision(2) <<  (m(1,i) - yy[i]) / (1000.*METER)
	   << " " << std::setw(8) << std::setprecision(2) <<  (m(2,i) - zz[i]) / (1000.*METER)
	   << endl;
      if ( abs(m(0,i) - xx[i]) > tolerance ||
	   abs(m(1,i) - yy[i]) > tolerance ||
	   abs(m(2,i) - zz[i]) > tolerance)
	fail = true;
    }
    if (fail)
      cerr << "***FAILURE: Difference from horizons > " << tolerance/(1000*METER) << " km" << endl;


    const double velocity_tolerance = 0.001 ; // v tolerance in m/s.
    {
      bool vfail = false;
      // Check last velocity of those extracted in bulk:
      int i = mjd.size() - 1;
      cout << "--- Velocity difference (m/s) ---" << endl;
      cout << " DT (yrs)     x         y       z " << endl;
      cout << future[i] << std::setprecision(4);
      for (int j=0; j<3; j++) {
	double diff = (velocity(j,i) - vv[j])*SECOND/METER;
	cout << "  " << diff;
	if (abs(diff) > velocity_tolerance) vfail = true;
      }
      cout << endl;

      if (vfail) {
	fail = true;
	cerr << "***FAILURE: Velocity difference from horizons***" << endl;
      }
    }
    
    // Get state at single time point, compare to horizons
    {
      double mjd1 = 2454466.5000; // = A.D. 2008-Jan-01 00:00:00.0000 TDB
      double tdb1 = (mjd1-JD2000)*DAY;
      astrometry::CartesianICRS x_horizons(3.815011333177522E+01,
			       1.361184665408451E+01,
			       6.545469772482830E+00);
      astrometry::CartesianICRS v_horizons(-9.655441061360712E-04/DAY,
			       2.345073633226195E-03/DAY,
			       1.124404442419986E-03/DAY);
      astrometry::CartesianICRS v_us;
      auto x_us = qb.position(tdb1, &v_us);

      x_us -= x_horizons;
      v_us -= v_horizons;
      v_us *= SECOND/METER;  // Convert to m/s
      {
	bool vfail = false;
	// Check last velocity of those extracted in bulk:
	cout << "--- Single-epoch position (km) and velocity (m/s) ---" << endl;
	cout << " DT (yrs)     x         y       z " << endl;
	cout << std::setprecision(2) << tdb1-s.tdb << std::setprecision(4) << " x:";
	for (int j=0; j<3; j++) {
	  double diff = x_us.getVector()[j];
	  cout << "  " << diff / (1e3*METER);
	  if (abs(diff) > tolerance) vfail = true;
	}
	cout << endl;
	cout << std::setprecision(2) << tdb1-s.tdb << std::setprecision(4) << " v:";
	for (int j=0; j<3; j++) {
	  double diff = v_us.getVector()[j];
	  cout << "  " << diff;
	  if (abs(diff) > velocity_tolerance) vfail = true;
	}
	cout << endl;

	if (vfail) {
	  fail = true;
	  cerr << "***FAILURE: Posn/velocity difference from horizons***" << endl;
	}
      }
    
    }
      

  } catch (std::runtime_error& e) {
    cerr << e.what() << endl;
    exit(1);
  }
  
  exit(fail ? 1 : 0);
}
