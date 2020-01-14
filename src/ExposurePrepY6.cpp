// Program to merge information on exposures into one table
// Y6A1 version.

#include <vector>
#include <iostream>
#include <utility>
#include <algorithm>

#include "Std.h"
#include "StringStuff.h"
#include "FitsTable.h"
#include "Astrometry.h"
#include "Ephemeris.h"

using namespace std;

template <class T1, class T2>
void indirect_sort(vector<T1>& v, const vector<T2>& key) {
  vector<std::pair<T2,T1>> vp;
  vp.reserve(v.size());
  for (int i=0; i<v.size(); i++)
    vp.push_back(std::make_pair(key[i],v[i]));
  std::sort(vp.begin(), vp.end());
  for (int i=0; i<v.size(); i++)
    v[i] = vp[i].second;
  return;
}
     
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
  "usage: ExposurePrep <astrometric exposures> <all exposures> <output table>\n"
  "  <astrometric exposures> is name of FITS file with astrometric data for exposures.\n"
  "  <all exposures> is name of FITS file with data for all exposures.\n"
  "  <output table> is name of FITS file for new merged table";

int
main(int argc, char *argv[]) {
  if (argc != 4) {
    cerr << usage << endl;
    exit(1);
  }

  string exposureFile = argv[1];
  string allExposureFile = argv[2];
  string outFile = argv[3];

  orbits::Ephemeris ephem;

  img::FTable exposureTable;
  try {
    FITS::FitsTable ff(exposureFile);
    vector<string> columnsToKeep;
    columnsToKeep.push_back("expnum");
    columnsToKeep.push_back("ra");
    columnsToKeep.push_back("dec");
    columnsToKeep.push_back("band");
    columnsToKeep.push_back("cov");
    columnsToKeep.push_back("mjd");
    
    exposureTable = ff.extract(0,-1, columnsToKeep);
  } catch (std::runtime_error& m) {
    cerr << "Error reading exposure table" << endl;
    quit(m);
  }

  cerr << "Read " << exposureTable.nrows() << " exposures' data." << endl;

  // Convert expnum from FITS "K" to "J" type
  std::vector<int> expnum(exposureTable.nrows(),0);
  for (int i=0; i<exposureTable.nrows(); i++) {
    LONGLONG expnumInLong;
    exposureTable.readCell(expnumInLong, "expnum", i);
    expnum[i] = expnumInLong;
  }
  exposureTable.eraseColumn("expnum");
  exposureTable.addColumn(expnum,"expnum");

  // Set up covwarn=true where cov is set to zero
  std::vector<bool> covwarn(exposureTable.nrows(),false);
  vector<double> cov;
  for (int i=0; i<exposureTable.nrows(); i++) {
    exposureTable.readCell(cov,"cov",i);
    if (cov[0]<=0.)
      covwarn[i] = true;
  }
  exposureTable.addColumn(covwarn,"covwarn");
  
  // Change band to FILTER
  vector<string> filter;
  exposureTable.readCells(filter, "band");
  exposureTable.eraseColumn("band");
  exposureTable.addColumn(filter,"filter",1,1);
  
  // Rename MJD to MJD_MID
  vector<double> mjd;
  exposureTable.readCells(mjd,"mjd");
  exposureTable.eraseColumn("mjd");
  exposureTable.addColumn(mjd,"mjd_mid");

  try {
    // Create new information about exposures
    std::vector<std::vector<double>> obspos; // Observatory ICRS coords
    std::vector<double> ecl_lat;
    std::vector<double> ecl_lon;
    std::vector<double> obs_ecl_lon;  // Ecliptic longitude of observatory

    // Covariance to insert where warning is set.
    const std::vector<double> DEFAULT_COV = {pow(15.,2.),pow(15.,2.),0.};

    double mjd, ra, dec;
    for (int i=0; i<exposureTable.nrows(); i++) {
      double mjd;
      exposureTable.readCell(mjd, "mjd_mid",i);
      // Get observatory position
      astrometry::CartesianICRS observatory = ephem.observatory(807, ephem.mjd2tdb(mjd));
      std::vector<double> v3 = {observatory[0],
				observatory[1],
				observatory[2]};
      obspos.push_back(v3);

      {
	// Get ecliptic longitude of observatory
	astrometry::CartesianEcliptic cc(observatory);
	astrometry::SphericalEcliptic ecl(cc);
	double lon,lat;
	ecl.getLonLat(lon, lat);
	obs_ecl_lon.push_back(lon/DEGREE);
      }
	
      {
	// Get ecliptic coordinates of pointing
	exposureTable.readCell(ra,"ra",i);
	exposureTable.readCell(dec,"dec",i);
	astrometry::SphericalICRS icrs(ra*DEGREE, dec*DEGREE);
	astrometry::SphericalEcliptic ecl(icrs);
	double lon,lat;
	ecl.getLonLat(lon, lat);
	ecl_lon.push_back(lon/DEGREE);
	ecl_lat.push_back(lat/DEGREE);
      }

      // Replace cov with default if warned
      bool warn;
      exposureTable.readCell(warn,"covwarn",i);
      if (warn)
	exposureTable.writeCell(DEFAULT_COV, "cov", i);
      
    } 
    exposureTable.addColumn(obspos,"observatory",3);
    exposureTable.addColumn(ecl_lon,"ecl_lon");
    exposureTable.addColumn(ecl_lat,"ecl_lat");
    exposureTable.addColumn(obs_ecl_lon,"obs_ecl_lon");
  } catch (std::runtime_error& m) {
    quit(m);
  }

  
  // Write the exposure table to disk
  try {
    FITS::FitsTable out(outFile,FITS::Create+FITS::OverwriteFile);
    out.copy(exposureTable);
  } catch (FITS::FITSError& m) {
    cerr << "Error writing first output table to " << outFile << endl;
    quit(m);
  }

  //////////////////////////////////////////////////////////////////
  // Now make secondary extension for exposures not in the above,
  // a.k.a. "non-astrometric" (which have no covariance matrix columns)
  //////////////////////////////////////////////////////////////////

  img::FTable allExposuresTable;
  try {
    FITS::FitsTable ff(allExposureFile);
    allExposuresTable = ff.extract();
  } catch (std::runtime_error& m) {
    cerr << "Error reading all-exposure table" << endl;
    quit(m);
  }

  // Record which exposures have astrometric information
  set<int> knownExposures;
  expnum.clear();
  exposureTable.readCells(expnum,"EXPNUM");
  knownExposures.insert(expnum.begin(), expnum.end());

  try {
    // Collect info on new exposures
    expnum.clear();
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
      allExposuresTable.readCell(raIn, "RA", i);
      allExposuresTable.readCell(decIn, "DEC", i);
      allExposuresTable.readCell(bandIn, "BAND", i);

      // Skip exposures for which filter name is >1 char
      stringstuff::stripWhite(bandIn);
      if (bandIn.size()>1) continue;
    
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

    // Sort all fields by expnum
    indirect_sort(mjd_mid,expnum);
    indirect_sort(ra,expnum);
    indirect_sort(dec,expnum);
    indirect_sort(filter,expnum);
    indirect_sort(observatory,expnum);
    indirect_sort(ecl_lon,expnum);
    indirect_sort(ecl_lat,expnum);
    indirect_sort(obs_ecl_lon,expnum);
    std::sort(expnum.begin(), expnum.end());

    // Write a new table to the exposure table
    FITS::FitsTable ff(outFile,FITS::Create+FITS::OverwriteHDU,2);
    img::FTable out = ff.use();
    out.addColumn(expnum,"expnum");
    out.addColumn(mjd_mid,"mjd_mid");
    out.addColumn(ra,"ra");
    out.addColumn(dec,"dec");
    out.addColumn(filter,"filter",1,1);
    out.addColumn(observatory,"observatory",3);
    out.addColumn(ecl_lon,"ecl_lon");
    out.addColumn(ecl_lat,"ecl_lat");
    out.addColumn(obs_ecl_lon,"obs_ecl_lon");
  } catch (std::runtime_error& m) {
    quit(m);
  }
  
  
  exit(0);
}
