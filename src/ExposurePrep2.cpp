// Add a second extension to exposure table containing
// time and pointing info for exposures not in the first table.

#include <vector>
#include <iostream>

#include "Std.h"
#include "StringStuff.h"
#include "FitsTable.h"
#include "Astrometry.h"
#include "Ephemeris.h"

using namespace std;

/** New output table will include:
expnum
mjd_mid
ra, dec (from initial header, less accurate)
observatory position (ICRS)
ecliptic lat, lon of pointing
ecliptic lon of observatory
**/

string usage =
  "Add a table of times/points for exposures without astrometric solutions\n"
  "to the FITS file holding info on exposures with astrometry.\n"
  "usage: ExposurePrep2 <exposure list input> <exposure FITS table to augment>\n"
  "  <exposure list input> is FITS table of all exposures from DESDB with MJD, etc.\n"
  "  <exposure FITS...> is name of existing astrometric exposure list to give new extension\n";

int
main(int argc, char *argv[]) {
  if (argc != 3) {
    cerr << usage << endl;
    exit(1);
  }
  string allExposureFile = argv[1];
  string targetExposureFile = argv[2];

  // Get info on all exposures
  img::FTable allExposuresTable;
  try {
    FITS::FitsTable ff(allExposureFile);
    allExposuresTable = ff.extract();
  } catch (std::runtime_error& m) {
    cerr << "Error reading all-exposure table" << endl;
    quit(m);
  }

  try {

    orbits::Ephemeris ephem;
    
    // Record which exposures have astrometric information
    set<int> knownExposures;
    {
      FITS::FitsTable ff(targetExposureFile,FITS::ReadOnly,1);
      vector<string> useColumns(1,"EXPNUM");
      img::FTable tab = ff.extract(0,-1,useColumns);
      vector<int> expnum;
      tab.readCells(expnum,"EXPNUM");
      knownExposures.insert(expnum.begin(), expnum.end());
    }

    // Collect info on new exposures
    vector<int> expnum;
    vector<double> mjd_mid;
    vector<string> filter;
    vector<vector<double>> observatory;  //ICRS observatory position
    vector<double> ra;
    vector<double> dec;
    vector<double> ecl_lat;
    vector<double> ecl_lon;
    vector<double> obs_ecl_lon;

  
    for (int i=0; i<allExposuresTable.nrows(); i++) {
      LONGLONG expnumInLong;
      // Skip any exposure in the first table
      allExposuresTable.readCell(expnumInLong, "EXPNUM", i);
      int expnumIn = expnumInLong;
      if (knownExposures.count(expnumIn)) continue;

      // Collect other info for this exposure
      double mjdIn, raIn, decIn;
      string bandIn;
      allExposuresTable.readCell(mjdIn, "MJD_MID", i);
      allExposuresTable.readCell(raIn, "RADEG", i);
      allExposuresTable.readCell(decIn, "DECDEG", i);
      allExposuresTable.readCell(bandIn, "BAND", i);
    
      expnum.push_back(expnumIn);
      mjd_mid.push_back(mjdIn);
      ra.push_back(raIn);
      dec.push_back(decIn);
      filter.push_back(bandIn);
    
      // Get observatory position
      astrometry::CartesianICRS obspos = ephem.observatory(807, ephem.mjd2tdb(mjdIn));
      std::vector<double> v3 = {obspos[0],
				obspos[1],
				obspos[2]};
      observatory.push_back(v3);
  

      // Get ecliptic longitude of observatory
      {
	astrometry::CartesianEcliptic cc(obspos);
	astrometry::SphericalEcliptic ecl(cc);
	double lon,lat;
	ecl.getLonLat(lon, lat);
	obs_ecl_lon.push_back(lon/DEGREE);
      }
	
      // Get ecliptic coordinates of pointing
      {
	astrometry::SphericalICRS icrs(raIn*DEGREE, decIn*DEGREE);
	astrometry::SphericalEcliptic ecl(icrs);
	double lon,lat;
	ecl.getLonLat(lon, lat);
	ecl_lon.push_back(lon/DEGREE);
	ecl_lat.push_back(lat/DEGREE);
      }

    } // End loop over input exposures

    // Write a new table to the exposure table
    FITS::FitsTable ff(targetExposureFile,FITS::Create+FITS::OverwriteHDU,2);
    img::FTable out = ff.use();
    out.addColumn(expnum,"expnum");
    out.addColumn(mjd_mid,"mjd_mid");
    out.addColumn(ra,"ra");
    out.addColumn(dec,"dec");
    out.addColumn(filter,"filter");
    out.addColumn(observatory,"observatory",3);
    out.addColumn(ecl_lon,"ecl_lon");
    out.addColumn(ecl_lat,"ecl_lat");
    out.addColumn(obs_ecl_lon,"obs_ecl_lon");
  } catch (std::runtime_error& m) {
    quit(m);
  }

  exit(0);
}

  
