// Predict positions and uncertainties from an ABG file.
#include "Pset.h"
#include "Fitter.h"
#include "Elements.h"
#include "StringStuff.h"
#include "FTable.h"
#include "FitsTable.h"
#include "Astrometry.h"
#include "Exposures.h"

#include <iostream>
#include <fstream>

const string usage =
  "PrintEphemeris: produce predictions and uncertainties for previously fitted orbit\n"
  "at specified time intervals, or at times specified in a file.\n\n"
  "Usage:\n"
  "  PrintEphemeris [parameter file...] [-<key> <value>...]\n"
  "  where any parameter file(s) given will be scanned first, then parameter\n"
  "   key/value pairs on cmd line will be read and override file values.\n"
  "  Program options are listed below.\n"
  "\n"
  "Input orbit fitting results are provided as ASCII file from fitter.save().\n"
  "Starting/ending dates, or dates in file can be given either as TDB or MJD.\n"
  "Dates in file can also be yyyy-mm-dd.dddd format.\n"
  "Interval is in days.\n"
  "Output will be to stdout.";

using namespace std;
using namespace orbits;

int main(int argc,
	 char *argv[])
{
  try {
    string ephemerisPath;
    double startDate, endDate, stepDate;
    string dateFile;
    int obscode;
    Pset parameters;
   
    {
      const int def=PsetMember::hasDefault;
      const int low=PsetMember::hasLowerBound;
      const int up=PsetMember::hasUpperBound;
      const int lowopen = low | PsetMember::openLowerBound;
      const int upopen = up | PsetMember::openUpperBound;

      parameters.addMember("dateFile",&dateFile, def,
			   "File holding TDBs or MJDs", "");
      parameters.addMember("startDate",&startDate, def,
			   "Start of ephem, TDB or MJD", 2015.0);
      parameters.addMember("endDate",&endDate, def,
			   "End of ephem, TDB or MJD", 2023.0);
      parameters.addMember("stepDate",&stepDate, def,
			   "Ephem interval (days)", 30.);
      parameters.addMember("ephemerisFile",&ephemerisPath, def,
			   "SPICE file (null=>environment)", "");
      parameters.addMember("obscode",&obscode, def | low,
			   "Observatory code (default: CTIO)", 807, 0);
    }
    parameters.setDefault();
    
    if (argc > 1 && (string(argv[1])=="-h" || string(argv[1])=="--help")) {
      cout << usage << endl;
      cout << "\nParameters:" << endl;
      parameters.dump(cerr);
      exit(1);
    }
    
    {
      // Read any parameter files
      int nPositional=0;
      for (int iarg=1; iarg < argc && argv[iarg][0]!='-'; iarg++) {
	std::ifstream ifs(argv[iarg]);
	if (!ifs) {
	  cerr << "Can't open parameter file " << argv[iarg] << endl;
	  cerr << usage << endl;
	  exit(1);
	}
	try {
	  parameters.setStream(ifs);
	} catch (std::runtime_error &m) {
	  cerr << "In file " << argv[iarg] << ":" << endl;
	  quit(m,1);
	}
	nPositional++;
      }
      // And read arguments from the remaining command line entries
      parameters.setFromArguments(argc, argv);
    }

    // Read the ephemeris
    Ephemeris eph(ephemerisPath);

    // Load a Fitter from stdin
    Fitter fit(eph);
    fit.restore(std::cin);

    int nobs=0;
    std::vector<double> tdb_in;
    
    if (dateFile.empty()) {
      // Parse the start and end dates
      if (startDate > 3000.) {
	// It's an MJD, convert to TDB (yrs past 2000)
	startDate = eph.mjd2tdb(startDate);
      } else {
	startDate -= 2000.;
      }
      if (endDate > 3000.) {
	// It's an MJD, convert to TDB (yrs past 2000)
	endDate = eph.mjd2tdb(endDate);
      } else {
	endDate -= 2000.;
      }
      nobs = int( floor((endDate-startDate)/(stepDate*DAY))) + 1;

      for (int i=0; i<nobs; i++) {
	tdb_in.push_back(startDate + i*stepDate*DAY);
      }
      
    } else {
      // Read MJD's or TDB's from a file
      std::ifstream ifs(dateFile);
      string buffer;
      while (stringstuff::getlineNoComment(ifs,buffer)) {
	string word = stringstuff::split(buffer).front();
	auto loc = word.find('-');
	if (loc==std::string::npos || loc==0) {
	  // Input is a number (potentially starting with negative sign)
	  std::istringstream iss(word);
	  double d;
	  iss >> d;
	  if (d > 3000.) {
	    // It's an MJD
	    tdb_in.push_back(eph.mjd2tdb(d));
	  } else {
	    // It's a TDB
	    tdb_in.push_back(d-2000.);
	  }
	} else {
	  // Input should be YMD
	  auto ymd = stringstuff::split(word,'-');
	  int y = atoi(ymd.front().c_str());
	  ymd.pop_front();
	  int m = atoi(ymd.front().c_str());
	  ymd.pop_front();
	  double d = atof(ymd.front().c_str());
	  astrometry::UT t(y,m,d);
	  tdb_in.push_back(eph.mjd2tdb(t.getMJD()));
	}
      }
      nobs = tdb_in.size();
    }
    if (nobs <= 0) {
      cerr << "No observing dates." << endl;
      exit(1);
    }
    
    cerr << "Calculating " << nobs << " positions" << endl;
      
    DVector tobs(nobs);
    DMatrix earth(nobs,3);
    astrometry::CartesianICRS xyz;
    for (int i=0; i<nobs; i++) {
      tobs[i] = tdb_in[i];
      xyz = eph.observatory(obscode, tobs[i]);
      earth.row(i) = xyz.getVector().transpose();
    }
    // Put data into frame
    tobs.array() -= fit.getFrame().tdb0;
    DMatrix xE = fit.getFrame().fromICRS(earth);

    // Destinations for prediction data
    DVector x(nobs);
    DVector y(nobs);
    DVector covXX(nobs);
    DVector covYY(nobs);
    DVector covXY(nobs);
    fit.predict(tobs, xE, &x, &y, &covXX, &covYY, &covXY);
    
    // Get major axes of error ellipses
    DVector tr2 = 0.5*(covXX + covYY);
    DVector det = covXX.array()*covYY.array() - covXY.array()*covXY.array();
    DVector a = sqrt(tr2.array() + sqrt(tr2.array()*tr2.array() - det.array()));
    // Convert to arcsec
    a /= ARCSEC;

    // Print a header
    cout << "#      UTC         RA         Dec    Error   Range    " << endl;
    cout << "# Y  M   D         deg        deg    asec  bary   geo " << endl;
    
    // Now stuff result back into array or write to screen
    for (int i=0; i<nobs; i++) {
      double ra, dec, tdb;
      astrometry::Gnomonic gn(x[i],y[i],fit.getFrame().orient);
      astrometry::SphericalICRS icrs(gn);
      icrs.getLonLat(ra, dec);
      ra /= DEGREE;
      dec /= DEGREE;

      // Let's also calculate bary and geocentric range
      astrometry::Vector3 x = fit.predictState(tobs[i]).x.getVector();
      double barycentricRange = sqrt(x.dot(x));
      tdb = tobs[i] + fit.getFrame().tdb0;
      astrometry::Vector3 dx = x - eph.observatory(obscode, tdb).getVector();
      double geocentricRange = sqrt(dx.dot(dx));
	
      
      // Reformat time as a date
      astrometry::UT ut;
      double mjd = eph.tdb2mjd(tdb);
      ut.setMJD(mjd);
      cout << fixed << setprecision(3);
      ut.writeYMD(cout);
      cout << " " << fixed << noshowpos << setprecision(6) << setw(10) << (ra<0 ? ra+360. : ra)
	   << " " << showpos << setw(10) << dec
	   << " " << noshowpos << setprecision(3) << a[i]
	   << " " << setprecision(2) << barycentricRange
	   << " " << geocentricRange
	   << endl;
    }

  } catch (std::runtime_error& e) {
    quit(e);
  }
  exit(0);
}
