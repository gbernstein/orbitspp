// Run orbit fits given sets of observations stored in FITS tables
#include "Pset.h"
#include "OrbitTypes.h"
#include "Fitter.h"
#include "Elements.h"
#include "FTable.h"
#include "FitsTable.h"
#include "Astrometry.h"

#include <iostream>

const string usage =
  "BulkElements: calculate orbital elements and their covariance from ABG\n"
  "values fit to data by BulkFit.\n"
  "Usage:\n"
  "  BulkElements [parameter file...] [-<key> <value>...]\n"
  "  where any parameter file(s) given will be scanned first, then parameter key/value\n"
  "  pairs on cmd line will be read and override file values.\n"
  "  Program options are listed below.\n"
  "\n"
  "The input orbit table (in FITS file) will on exit have additional\n"
  "columns for ELEMENTS and ELCOV.  Stored elements use radians\n"
  "for angles.";

using namespace std;
using namespace orbits;

int main(int argc,
	 char *argv[])
{
  try {
    string ephemerisPath;
    string orbitPath;
    Pset parameters;
   
    {
      const int def=PsetMember::hasDefault;
      const int low=PsetMember::hasLowerBound;
      const int up=PsetMember::hasUpperBound;
      const int lowopen = low | PsetMember::openLowerBound;
      const int upopen = up | PsetMember::openUpperBound;

      parameters.addMember("orbitFile",&orbitPath, def,
			   "FITS file for input/output orbital parameters", "");
      parameters.addMember("ephemerisFile",&ephemerisPath, def,
			   "SPICE file (null=>environment", "");
    }
    parameters.setDefault();

    if (argc<2 || string(argv[1])=="-h" || string(argv[1])=="--help") {
      cout << usage << endl;
      parameters.dump(cerr);
      exit(1);
    }
    
    {
      // Read any parameter files
      int nPositional=0;
      for (int iarg=1; iarg < argc && argv[iarg][0]!='-'; iarg++) {
	ifstream ifs(argv[iarg]);
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

    FITS::FitsTable ft(orbitPath,FITS::ReadWrite,1);
    auto orbitTable = ft.use();

    // Set reference frame from orbit table header info
    Frame frame;
    {
      double ra0, dec0, pa0, tdb0, x0, y0, z0;
      orbitTable.header()->getValue("RA0",ra0);
      orbitTable.header()->getValue("DEC0",dec0);
      orbitTable.header()->getValue("TDB0",tdb0);
      orbitTable.header()->getValue("PA0",pa0);
      orbitTable.header()->getValue("X0",x0);
      orbitTable.header()->getValue("Y0",y0);
      orbitTable.header()->getValue("Z0",z0);
      astrometry::CartesianICRS origin(x0,y0,z0);
      astrometry::SphericalICRS pole(ra0*DEGREE, dec0*DEGREE);
      astrometry::Orientation orient(pole, pa0*DEGREE);
      frame = Frame(origin, orient, tdb0);
    }

    // Make new columns
    {
      vector<vector<double>> vv;
      vv.push_back(vector<double>(6,0.));
      orbitTable.addColumn(vv,"ELEMENTS");
      vv.clear();
      vv.push_back(vector<double>(36,0.));
      orbitTable.addColumn(vv,"ELCOV");
    }

    // Begin reading input data.  One Fitter should do for everything.
    Fitter fit(eph, Gravity::GIANTS);
    fit.setFrame(frame);

    for (int row=0; row<orbitTable.nrows(); row++) {
      int flag;
      orbitTable.readCell(flag, "FLAGS", row);
      if (flag>0) {
	orbitTable.writeCell(vector<double>(6,0.),"ELEMENTS",row);
	orbitTable.writeCell(vector<double>(36,0.),"ELCOV",row);
	continue; // Do not calculate elements for bad orbit
      }
      vector<double> v;
      orbitTable.readCell(v,"ABG",row);
      ABG abg;
      for (int i=0; i<6; i++) abg[i] = v[i];
      orbitTable.readCell(v,"ABGCOV",row);
      ABGCovariance cov;
      for (int i=0; i<6; i++)
	for (int j=0; j<6; j++)
	  cov(i,j) = v[6*i+j];
      fit.setABG(abg,cov);
      Elements e= fit.getElements();
      v.resize(6);
      for (int i=0; i<6; i++) v[i]=e[i];
      orbitTable.writeCell(v,"ELEMENTS",row);
      ElementCovariance ecov = fit.getElementCovariance();
      v.resize(36);
      for (int i=0; i<6; i++)
	for (int j=0; j<6; j++)
	  v[6*i+j] = ecov(i,j);
      orbitTable.writeCell(v,"ELCOV",row);
    }
  } catch (std::runtime_error& e) {
    quit(e);
  }
  exit(0);
}
