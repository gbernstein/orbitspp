// Convert abg or aei to ICRS state vectors
#include "Pset.h"
#include "OrbitTypes.h"
#include "Fitter.h"
#include "Elements.h"
#include "FTable.h"
#include "FitsTable.h"
#include "Astrometry.h"

#include <iostream>

const string usage =
  "BulkState: calculate state vectors and their uncertainties from the\n"
  "           ABG values or elements in a table.  Use ABG if both\n"
  "           are present.  State is given for the Frame's reference time.\n"
  "Usage:\n"
  "  BulkState [parameter file...] [-<key> <value>...]\n"
  "  where any parameter file(s) given will be scanned first, then parameter key/value\n"
  "  pairs on cmd line will be read and override file values.\n"
  "  Program options are listed below.\n"
  "\n"
  "The input orbit table (in FITS file) will on exit have additional\n"
  "columns for XV and XVCOV.  State vectors are barycentric ICRS\n"
  "in units of AU and Julian years";
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
			   "FITS file for input/output orbital info", "");
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

    // See what kind of input we're going to use
    bool useABG = false;
    if ( orbitTable.hasColumn("ABG") && orbitTable.hasColumn("ABGCOV"))
      useABG = true;
    else if ( orbitTable.hasColumn("ELEMENTS") && orbitTable.hasColumn("ELEMENTCOV"))
      useABG = false;  // Use elements
    else {
      cerr << "Can find neither columns for ABG nor for elements in table"
	   << endl;
      exit(1);
    }

    // Make new columns
    {
      vector<vector<double>> vv;
      vv.push_back(vector<double>(6,0.));
      orbitTable.addColumn(vv,"XV");
      vv.clear();
      vv.push_back(vector<double>(36,0.));
      orbitTable.addColumn(vv,"XVCOV");
    }

    // Begin reading input data.  One Fitter should do for everything.
    Fitter fit(eph, Gravity::GIANTS);
    fit.setFrame(frame);

    for (int row=0; row<orbitTable.nrows(); row++) {
      int flag = 0;
      if (orbitTable.hasColumn("FLAGS")) {
      // int flag;
      orbitTable.readCell(flag, "FLAGS", row);
    }
      if (flag>0) {
	orbitTable.writeCell(vector<double>(6,0.),"XV",row);
	orbitTable.writeCell(vector<double>(36,0.),"XVCOV",row);
	continue; // Do not calculate elements for bad orbit
      }

      State s;
      Matrix66 sCov;

      if (useABG) {
	vector<double> v;
	orbitTable.readCell(v,"ABG",row);
	ABG abg;
	for (int i=0; i<6; i++) abg[i] = v[i];
	orbitTable.readCell(v,"ABGCOV",row);
	ABGCovariance cov;
	for (int i=0; i<6; i++)
	  for (int j=0; j<6; j++)
	    cov(i,j) = v[6*i+j];

	Vector3 xFrame;
	Vector3 vFrame;
	// Get state vector and covariance in Frame
	// at reference epoch of frame
	abg.getState(0., xFrame, vFrame); 
	// Derivative of state vector in reference frame wrt ABG:
	Matrix66 dSdABG_frame = abg.getStateDerivatives();
	 // Rotate x and v derivatives into ICRS
	Matrix66 dSdABG;
	// Note that Frame is working with Nx3 arrays so we need
	// to put xyz on the 2nd index.
	DMatrix tmp = dSdABG_frame.subMatrix(0,3,0,6).transpose();
	dSdABG.subMatrix(0,3,0,6) = frame.toICRS(tmp , true).transpose();
	tmp = dSdABG_frame.subMatrix(3,6,0,6).transpose();
	dSdABG.subMatrix(3,6,0,6) = frame.toICRS(tmp, true).transpose();

	// And now transform to ICRS
	s.x = frame.toICRS(xFrame);
	s.v = frame.toICRS(vFrame,true);
	s.tdb = frame.tdb0;
	sCov = dSdABG * cov * dSdABG.transpose();
	
      } else {
	// Start with elements
	vector<double> v(6);
	orbitTable.readCell(v,"ELEMENTS",row);
	Elements e;
	for (int i=0; i<6; i++) e[i] = v[i];
	v.resize(36);
	orbitTable.readCell(v,"ELEMENTCOV",row);
	ElementCovariance eCov;
	for (int i=0; i<6; i++)
	  for (int j=0; j<6; j++)
	    eCov(i,j) = v[6*i+j];

	// Get state vector from elements
	s = getState(e,frame.tdb0);
	// And the derivatives - need to invert
	auto dEdS = getElementDerivatives(s);
	// Invert to get dSdE
	// ?? Catch a degenerate matrix?
	Matrix66 dSdE = dEdS.inverse();
	sCov = dSdE * eCov * dSdE.transpose();
      }

      vector<double> v(6);
      for (int i=0; i<3; i++) {
	v[i]=s.x[i];
	v[i+3] = s.v[i];
      }
      orbitTable.writeCell(v,"XV",row);
      v.resize(36);
      for (int i=0; i<6; i++)
	for (int j=0; j<6; j++)
	  v[6*i+j] = sCov(i,j);
      orbitTable.writeCell(v,"XVCOV",row);
    }
  } catch (std::runtime_error& e) {
    quit(e);
  }
  exit(0);
}
