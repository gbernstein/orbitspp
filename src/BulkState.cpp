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
  "columns for XV and XVINVCOV.  State vectors are barycentric ICRS\n"
  "in units of AU and Julian years.\n"
  "\n"
  "If the input ABG or ELEMENTS have degeneracies in the conversion to state vector,\n"
  "or if ELEMENTCOV is singular, then the output XVINVCOV will be zeros.  If the\n"
  "input orbit is unbound, then both XV and XVINVCOV will be zeros.";
using namespace std;
using namespace orbits;

int main(int argc,
	 char *argv[])
{
  try {
    string ephemerisPath;
    string orbitPath;
    bool fixOrient;
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
			   "SPICE file (null=>environment)", "");
      parameters.addMember("fixOrient",&fixOrient, def,
			   "Apply fix for MergeOrbits FRAME bug", false);
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

    // See what kind of input we're going to use
    bool useABG = false;
    if ( orbitTable.hasColumn("ABG") && orbitTable.hasColumn("ABGINVCOV"))
      useABG = true;
    else if ( orbitTable.hasColumn("ABG") && orbitTable.hasColumn("ABGCOV"))
    	useABG = true;
    else if ( orbitTable.hasColumn("ELEMENTS") && orbitTable.hasColumn("ELEMENTCOV"))
      useABG = false;  // Use elements
    else {
      cerr << "Can find neither columns for ABG nor for elements in table"
	   << endl;
      exit(1);
    }

    // Set reference frame from orbit table header info
    Frame frame;
    // This variable will be true if we are to get the
    // reference frame for each row from a column in the table,
    // vs having one common one in the header of the ORBIT table
    bool needIndividualFrames = false; 
    {
      // See if there is a reference frame in the header
      double ra0, dec0, pa0, tdb0, x0, y0, z0;
      if (orbitTable.header()->getValue("RA0",ra0) &&
	  orbitTable.header()->getValue("DEC0",dec0) &&
	  orbitTable.header()->getValue("TDB0",tdb0) &&
	  orbitTable.header()->getValue("PA0",pa0) &&
	  orbitTable.header()->getValue("X0",x0) && 
	  orbitTable.header()->getValue("Y0",y0) &&
	  orbitTable.header()->getValue("Z0",z0)) {
	// The Frame is in the header
	astrometry::CartesianICRS origin(x0,y0,z0);
	astrometry::SphericalICRS pole(ra0*DEGREE, dec0*DEGREE);
	astrometry::Orientation orient(pole, pa0*DEGREE);
	frame = Frame(origin, orient, tdb0);
      } else {
	// Need to get Frame from a column in the table
	needIndividualFrames = true;
      }
    }

    if (needIndividualFrames && !orbitTable.hasColumn("FRAME")) {
      cerr << "ERROR: Reference frame found neither in header nor in a column" << endl;
      exit(1);
    }


    // Make new columns
    {
      // Get rid of old columns if they are there
      if (orbitTable.hasColumn("XV"))
	orbitTable.eraseColumn("XV");
      if (orbitTable.hasColumn("XVINVCOV"))
	orbitTable.eraseColumn("XVINVCOV");
      
      vector<vector<double>> vv;
      vv.push_back(vector<double>(6,0.));
      orbitTable.addColumn(vv,"XV", 6);
      vv.clear();
      vv.push_back(vector<double>(36,0.));
      orbitTable.addColumn(vv,"XVINVCOV", 36);
    }

    for (int row=0; row<orbitTable.nrows(); row++) {
      int flag = 0;
      if (orbitTable.hasColumn("FLAGS")) {
	orbitTable.readCell(flag, "FLAGS", row);
      }
      if (flag>0) {
	orbitTable.writeCell(vector<double>(6,0.),"XV",row);
	orbitTable.writeCell(vector<double>(36,0.),"XVINVCOV",row);
	continue; // Do not calculate elements for bad orbit
      }

      
      State s;
      Matrix66 sInvCov;
      bool singularCov = false;  // Set this for bad inversion

      if (useABG) {
	vector<double> v;
	orbitTable.readCell(v,"ABG",row);
	ABG abg;
	for (int i=0; i<6; i++) abg[i] = v[i];
	orbitTable.readCell(v,"ABGINVCOV",row);
	Matrix66 aInvCov;
	for (int i=0; i<6; i++)
	  for (int j=0; j<6; j++)
	    aInvCov(i,j) = v[6*i+j];

	Vector3 xFrame;
	Vector3 vFrame;
	// Get state vector and covariance in Frame
	// at reference epoch of frame
	abg.getState(0., xFrame, vFrame);
	  
	// Derivative of state vector in reference frame wrt ABG:
	Matrix66 dSdABG_frame = abg.getStateDerivatives();
	 // Rotate x and v derivatives into ICRS
	Matrix66 dSdABG;

	if (needIndividualFrames) {
	  // Get reference frame from a column in the table
	  vector<double> v(7);
	  orbitTable.readCell(v,"FRAME",row);
	  astrometry::CartesianICRS origin(v[3],v[4],v[5]);
	  astrometry::SphericalICRS pole(v[0],v[1]);
	  astrometry::Orientation orient(pole, v[2]);
	  /* Do the following if input file has messed-up frames*/
	  if (fixOrient) orient.alignToEcliptic();
	  frame = Frame(origin, orient, v[6]);
	}

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

	/**/ if (false) { //** row==778) {
	  auto el = getElements(s);
	  cerr << "State from ABG:      " << s << endl;
	  cerr << "Elements from ABG: " << el << endl;
	  vector<double> ee(6);
	  orbitTable.readCell(ee,"ELEMENTS",row);
	  for (int i=0; i<6; i++) el[i] = ee[i];
	  cerr << "Elements from Table: " <<  el << endl;
	  cerr << "State from ELEMENTS: " << getState(el,frame.tdb0) << endl;
	  Fitter fit(eph);
	  fit.setFrame(frame);
	  fit.setABG(abg, aInvCov.inverse());
	  el = fit.getElements();
	  cerr << "Elements from Fitter: " <<  el << endl;
	}
	// We need to invert the derivatives to transform
	// the inverse covariance
	Eigen::FullPivLU<Eigen::Matrix<double,6,6>> lu(dSdABG);
	if (lu.isInvertible()) {
	  Matrix66 dABGdS = lu.inverse();
	  sInvCov = dABGdS.transpose() * aInvCov * dABGdS;
	} else {
	  singularCov = true;
	}
	
      } else {
	// Start with elements
	vector<double> v(6);
	orbitTable.readCell(v,"ELEMENTS",row);
	Elements e;
	for (int i=0; i<6; i++) e[i] = v[i];
	try {
	  // Get state vector from elements
	  if (needIndividualFrames) {
	    // Get reference frame from a column in the table
	    vector<double> v(7);
	    orbitTable.readCell(v,"FRAME",row);
	    s = getState(e,v[6]);  // Last element is TDB0
	  } else {
	    s = getState(e,frame.tdb0);
	  }
	  // Get covariance
	  v.resize(36);
	  orbitTable.readCell(v,"ELEMENTCOV",row);
	  Matrix66 eCov;
	  for (int i=0; i<6; i++)
	    for (int j=0; j<6; j++)
	      eCov(i,j) = v[6*i+j];
	  Eigen::LDLT<Eigen::Matrix<double,6,6>> chol(eCov);
	  if (!singularCov && chol.info()==Eigen::Success && chol.isPositive()) {
	    // Get state derivatives
	    auto dEdS = getElementDerivatives(s);
	    Matrix66 tmp = chol.solve(dEdS);
	    sInvCov = dEdS.transpose() * tmp;
	    // Above does this: sInvCov = dEdS.transpose() * eInvCov * dEdS;
	  } else {
	    // Bad cov matrix from elements
	    singularCov = true;
	  }
	} catch (std::runtime_error& e) {
	  // Catch non-elliptical orbits
	  s.x = Vector3(0.);
	  s.v = Vector3(0.);
	  singularCov = true;
	}
      }  // End ABG/ELEMENT choice

      vector<double> v(6);
      for (int i=0; i<3; i++) {
	v[i]=s.x[i];
	v[i+3] = s.v[i];
      }
      orbitTable.writeCell(v,"XV",row);
      v.resize(36);
      for (int i=0; i<6; i++)
	for (int j=0; j<6; j++)
	  v[6*i+j] = singularCov ? 0. : sInvCov(i,j);
      orbitTable.writeCell(v,"XVINVCOV",row);
    } // end row loop
  } catch (std::runtime_error& e) {
    quit(e);
  }
  exit(0);
}
