// Predict positions given orbital elements in a FITS file
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
  "BulkEphemeris: produce predicted positions on a sequence of dates\n"
  "for a bunch of objects with elements given in a FITS table.\n\n"
  "Usage:\n"
  "  BulkEphemeris [parameter file...] [-<key> <value>...]\n"
  "  where any parameter file(s) given will be scanned first, then parameter\n"
  "   key/value pairs on cmd line will be read and override file values.\n"
  "  Program options are listed below.\n"
  "\n"
  "Starting/ending dates can be given either as TDB or MJD.\n"
  "Dates in file can also be yyyy-mm-dd.dddd format.\n"
  "Interval is in days.\n"
  "Output will be to stdout, or in a FITS table if a name is given.";

using namespace std;
using namespace orbits;

int main(int argc,
	 char *argv[])
{
  try {
    string ephemerisPath;
    double startDate, endDate, stepDate;
    string dateFile;
    string elementFile, outFile;
    int obscode;
    Pset parameters;
   
    {
      const int def=PsetMember::hasDefault;
      const int low=PsetMember::hasLowerBound;
      const int up=PsetMember::hasUpperBound;
      const int lowopen = low | PsetMember::openLowerBound;
      const int upopen = up | PsetMember::openUpperBound;

      parameters.addMember("elementFile",&elementFile, 0,
			   "Orbital element table");
      parameters.addMember("dateFile",&dateFile, def,
			   "File holding TDBs or MJDs", "");
      parameters.addMember("outFile",&outFile, def,
			   "Output FITS filename for ephemeris", "");
      parameters.addMember("startDate",&startDate, def,
			   "Start of ephem, TDB or MJD", 2021.0);
      parameters.addMember("endDate",&endDate, def,
			   "End of ephem, TDB or MJD", 2033.0);
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
    Ephemeris ephem(ephemerisPath);

    // Open the table of input elements
    img::FTable eTable;
    {
      FITS::FitsTable ft(elementFile);
      eTable = ft.extract();
    }

    // Flag set if printing output to stdout instead of table
    bool textOut = outFile.empty();

    // Make an output arrays (in case writing a table)
    vector<double> tdb, ra, dec, helio, geo, phase, phasePA;
    vector<int> sourceID;

    int nobs=0;
    std::vector<double> tdb_in;
    
    if (dateFile.empty()) {
      // Parse the start and end dates
      if (startDate > 3000.) {
	// It's an MJD, convert to TDB (yrs past 2000)
	startDate = ephem.mjd2tdb(startDate);
      } else {
	startDate -= 2000.;
      }
      if (endDate > 3000.) {
	// It's an MJD, convert to TDB (yrs past 2000)
	endDate = ephem.mjd2tdb(endDate);
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
	    tdb_in.push_back(ephem.mjd2tdb(d));
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
	  tdb_in.push_back(ephem.mjd2tdb(t.getMJD()));
	}
      }
      nobs = tdb_in.size();
    }
    if (nobs <= 0) {
      cerr << "No observing dates." << endl;
      exit(1);
    }
    
    cerr << "Calculating " << nobs << " dates" << endl;
      
    DVector tobs(nobs);
    DMatrix earth(nobs,3);
    astrometry::CartesianICRS xyz;
    for (int i=0; i<nobs; i++) {
      tobs[i] = tdb_in[i];
      xyz = ephem.observatory(obscode, tobs[i]);
      earth.row(i) = xyz.getVector().transpose();
    }

    if (textOut) {
      // Print a header
      cout << "# Source    MJD    RA         Dec    bary  geo  phase  PA " << endl;
      cout << "# ------------------------------------------------------- " << endl;
    }
    // Begin loop through all objects in file
    int nobj = eTable.nrows();
    cerr << "Looping through " << nobj << " orbits" << endl;
    for (int source=0; source < nobj; source++) {
      vector<double> ve;
      double epoch;
      eTable.readCell(ve, "elements", source);
      eTable.readCell(epoch, "epoch", source);
      Elements el;
      // Put into standard form
      el[Elements::A] = ve[0];
      el[Elements::E] = ve[1];
      el[Elements::I] = ve[2]*DEGREE;
      el[Elements::LAN] = ve[3]*DEGREE;
      el[Elements::AOP] = ve[4]*DEGREE;
      el[Elements::TOP] = ephem.mjd2tdb(ve[5]);
      // Build trajectory - elements from MPC are heliocentric
      Trajectory traj(ephem, getState(el, ephem.mjd2tdb(epoch), true, &ephem),
		      GIANTS, 2*DAY);  // Integrate with giant planets at time step 2d
      
      // Loop over observations (not super efficient but good enough...)
      double ra0, dec0;
      for (int i=0; i<nobs; i++) {
	astrometry::CartesianICRS observer(earth.row(i));
	auto posn = traj.observe(tobs[i], observer);
	auto circ = traj.getCircumstances(tobs[i], observer.getVector());
	posn.getLonLat(ra0,dec0);
	if (textOut) {
	  cout << source << " "
	       << tobs[i] << " "
	       << ra0/DEGREE << " "
	       << dec0/DEGREE << " "
	       << circ.sunDistance << " "
	       << circ.obsDistance << " "
	       << circ.phaseAngle / DEGREE << " "
	       << circ.phasePA / DEGREE
	       << endl;
	} else {
	  // Append to arrays
	  sourceID.push_back(source);
	  tdb.push_back(tobs[i]);
	  ra.push_back(ra0/DEGREE);
	  dec.push_back(dec0/DEGREE);
	  helio.push_back(circ.sunDistance);
	  geo.push_back(circ.obsDistance);
	  phase.push_back(circ.phaseAngle/DEGREE);
	  phasePA.push_back(circ.phasePA/DEGREE);
	}
      } // End observation loop
    } // end object loop

    if (!textOut) {
      // Write a results table
      FITS::FitsTable ft(outFile,FITS::OverwriteFile);
      auto tab = ft.use();
      tab.addColumn(sourceID,"source");
      tab.addColumn(tdb,"tdb");
      tab.addColumn(ra,"ra");
      tab.addColumn(dec,"dec");
      tab.addColumn(helio,"helio");
      tab.addColumn(geo,"geo");
      tab.addColumn(phase,"phase");
      tab.addColumn(phasePA,"phasepa");
    }
    
    
  } catch (std::runtime_error& e) {
    quit(e);
  }
  exit(0);
}
