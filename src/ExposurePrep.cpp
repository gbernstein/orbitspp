// Program to merge information on exposures into one table.

#include <vector>
#include <iostream>

#include "Std.h"
#include "StringStuff.h"
#include "FitsTable.h"
#include "Astrometry.h"
#include "Ephemeris.h"

using namespace std;

/** Output table will include:
expnum
mjd_mid
ra, dec
astrometric covariance matrix
covariance warning flag
observatory position (ICRS)
ecliptic lat, lon of pointing
"elongation" = ecliptic longitude minus observatory lon
 (or maybe just include ecliptic lon of observatory)
**/

string usage =
  "Merge information from astrometric files and ephemeris to create\n"
  "a summary table of all exposures for use in orbit fitting.\n"
  "usage: ExposurePrep <exposure table> <ccd table> <output table>\n"
  "  <exposure table> is name of FITS file with astrometric data for exposures.\n"
  "  <ccd table> is name of FITS file with astrometric data for CCDs.\n"
  "  <output table> is name of FITS file for new merged table";

int
main(int argc, char *argv[]) {
  if (argc != 4) {
    cerr << usage << endl;
    exit(1);
  }

  string exposureFile = argv[1];
  img::FTable exposureTable;
  try {
    FITS::FitsTable ff(exposureFile);
    vector<string> columnsToKeep;
    columnsToKeep.push_back("expnum");
    columnsToKeep.push_back("ra");
    columnsToKeep.push_back("dec");
    columnsToKeep.push_back("filter");
    columnsToKeep.push_back("cov");
    columnsToKeep.push_back("covwarn");
    
    exposureTable = ff.extract(0,-1, columnsToKeep);
  } catch (std::runtime_error& m) {
    cerr << "Error reading exposure table" << endl;
    quit(m);
  }
  
  cerr << "Read " << exposureTable.nrows() << " exposures' data." << endl;
  
  // Get MJD-MID from the CCD table
  string ccdFile = argv[2];
  try {
    FITS::FitsTable ff(ccdFile);
    vector<string> columnsToKeep;
    columnsToKeep.push_back("expnum");
    columnsToKeep.push_back("mjd_mid");
    auto ccdTable = ff.extract(0,-1, columnsToKeep);

    cerr << "Read " << ccdTable.nrows() << " CCDs' data." << endl;

    // Read all of the mjd-mid values into a map
    map<int,double> mjdmap;
    for (int i=0; i<ccdTable.nrows(); i++) {
      int expnum;
      double mjd;
      ccdTable.readCell(expnum,"expnum",i);
      ccdTable.readCell(mjd,"mjd_mid",i);
      mjdmap[expnum] = mjd;
    }

    // And put them into the exposure table
    vector<double> mjd(exposureTable.nrows(), 0.);
    for (int i=0; i<exposureTable.nrows(); i++) {
      int expnum;
      exposureTable.readCell(expnum, "expnum",i);
      auto iter = mjdmap.find(expnum);
      if (iter!=mjdmap.end())
	mjd[i] = iter->second;
    }
    exposureTable.addColumn(mjd, "mjd_mid");
  } catch (std::runtime_error& m) {
    cerr << "Error reading ccd table" << endl;
    quit(m);
  }

  string outFile = argv[3];
  exit(0);
}

  