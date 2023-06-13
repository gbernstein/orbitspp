// Program to merge information on exposures into one table.

#include <vector>
#include <iostream>

#include "Std.h"
#include "StringStuff.h"
#include "FitsTable.h"
#include "Astrometry.h"
#include "Ephemeris.h"

using namespace std;

/** 
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
  "Add information on observatory position and ecliptic coords\n"
  "to a table of exposure data.\n"
  "  Input table should have columns for ra, dec, mjd_mid.\n"
  "  Output table will have additional columns for cov, covwarn\n"
  "if they were not there before, and will gain/overwrite columns for:\n"
  "  observatory = ICRS cartesian coordinates (AU) of observatory\n"
  "  ecl_lon, ecl_lat = ecliptic coords (degrees) of sightline\n"
  "  obs_ecl_lat = ecliptic longitude of observatory.\n"
  "usage: AddObservatory <exposure table> <obsCode>\n"
  "  <exposure table> is name of FITS file with exposure data to be augmented.\n"
  "  <obsCode> is the observatory code (see observatories.dat)";


int
main(int argc, char *argv[]) {
  if (argc < 2) {
    cerr << usage << endl;
    exit(1);
  }
  int obsCode = 807;

  string exposureFile = argv[1];

  if (argc > 2) {
      obsCode = atoi(argv[2]); 

    }
  else {
      cerr << "## Assuming observatory is CTIO/DECam" << endl;
    }

  try {
    orbits::Ephemeris ephem;

    FITS::FitsTable ff(exposureFile,FITS::ReadWrite);
    img::FTable table=ff.use();
    //make sure it has the needed columns
    if (!(table.hasColumn("ra") &&
	  table.hasColumn("dec") &&
	  (table.hasColumn("mjd_mid") ||
	   table.hasColumn("mjd-mid") ||
	   table.hasColumn("mjd")))) {
      cerr << "Input table is missing require column ra, dec, or mjd_mid" << endl;
      exit(1);
    }

      cerr << obsCode << endl;

    // See if it has cov, covwarn - if not, add them as we go
    bool addCov = !table.hasColumn("cov");
    if (addCov) {
      vector<vector<double>> vv;
      table.addColumn(vv, "cov", 3);  // fixed 3-element vectors
    }
    bool addCovWarn = !table.hasColumn("covwarn");
    if (addCovWarn) {
      vector<bool> v;
      table.addColumn(v,"covwarn");
    }
    
    vector<double> v3(3);
    {
      // Add column for observatory position
      vector<vector<double>> vv;
      if (table.hasColumn("observatory"))
	  table.eraseColumn("observatory");
      table.addColumn(vv, "observatory", 3);  // fixed 3-element vectors
    }
    {
      vector<double> v;
      if (table.hasColumn("ecl_lon"))
	  table.eraseColumn("ecl_lon");
      if (table.hasColumn("ecl_lat"))
	  table.eraseColumn("ecl_lat");
      if (table.hasColumn("obs_ecl_lon"))
	  table.eraseColumn("obs_ecl_lon");
      table.addColumn(v, "ecl_lon");
      table.addColumn(v, "ecl_lat");
      table.addColumn(v, "obs_ecl_lon");
    }

    // calculate new info
    const std::vector<double> covDefault(3,0.);
    for (int row=0; row<table.nrows(); row++) {
      double mjd, ra, dec;
      if (table.hasColumn("mjd_mid")) {
	table.readCell(mjd,"mjd_mid",row);
      } else if (table.hasColumn("mjd-mid")) {
	table.readCell(mjd,"mjd-mid",row);
      } else if (table.hasColumn("mjd")) {
	table.readCell(mjd,"mjd",row);
      }
      table.readCell(ra,"ra",row);
      table.readCell(dec,"dec",row);

      if (addCov)
	table.writeCell(covDefault,"cov",row);
      if (addCovWarn)
	table.writeCell(true,"covwarn",row);

      // Get observatory position
      astrometry::CartesianICRS observatory = ephem.observatory(obsCode, ephem.mjd2tdb(mjd));
      std::vector<double> v3 = {observatory[0],
				observatory[1],
				observatory[2]};
      table.writeCell(v3,"observatory",row);

      {
	// Get ecliptic longitude of observatory
	astrometry::CartesianEcliptic cc(observatory);
	astrometry::SphericalEcliptic ecl(cc);
	double lon,lat;
	ecl.getLonLat(lon, lat);
	table.writeCell(lon/DEGREE, "obs_ecl_lon",row);
      }
	
      {
	// Get ecliptic coordinates of pointing
	table.readCell(ra,"ra",row);
	table.readCell(dec,"dec",row);
	astrometry::SphericalICRS icrs(ra*DEGREE, dec*DEGREE);
	astrometry::SphericalEcliptic ecl(icrs);
	double lon,lat;
	ecl.getLonLat(lon, lat);
	table.writeCell(lon/DEGREE,"ecl_lon",row);
	table.writeCell(lat/DEGREE,"ecl_lat",row);
      }
    } // End row loop
    // Table will write back to FITS now.
  } catch (std::runtime_error& m) {
    quit(m);
  }
  
  exit(0);
}

  
