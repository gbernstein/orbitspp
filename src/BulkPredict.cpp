// Run orbit fits given sets of observations stored in FITS tables
#include "Pset.h"
#include "Fitter.h"
#include "Elements.h"
#include "StringStuff.h"
#include "FTable.h"
#include "FitsTable.h"
#include "Astrometry.h"
#include "Exposures.h"

#include <iostream>

const string usage =
  "BulkPredict: produce predictions and uncertainties for previously fitted orbits\n"
  "at a list of times or expnums.\n\n"
  "Usage:\n"
  "  FitsFitter [parameter file...] [-<key> <value>...]\n"
  "  where any parameter file(s) given will be scanned first, then parameter\n"
  "   key/value pairs on cmd line will be read and override file values.\n"
  "  Program options are listed below.\n"
  "\n"
  "Input orbit fitting results are provided as binary FITS table from BulkFitter.cpp\n"
  "The input (orbit ID, expnum) or (orbit ID, MJD) pairs can be provided either\n"
  "from the input observation file, at at stdin. More efficient if all observations\n"
  "for a given orbit are contiguous in the input. FITS input file should have columns\n"
  "for ORBIT_ID, and EXPNUM (=DECam exposure number, in which case exposure file must\n"
  "be provided) or MJD.  Input at stdin should have one pair of numbers per line.\n"
  "\n"
  "Output will be either to FITS file or stdout, with input info followed by a\n"
  "FLAG column giving status of input orbit (0=valid, others are invalid),\n"
  "then RA, DEC of predicted position in decimal degrees, and ERROR_A, the\n"
  "semimajor axis of error ellipse, in arcseconds.";

using namespace std;
using namespace orbits;

// Possible flags for failures:
const int NO_ORBIT = 32;

int main(int argc,
	 char *argv[])
{
  try {
    string ephemerisPath;
    string exposurePath;
    string observationPath;
    string orbitPath;
    string predictPath;
    double bindingFactor;
    double gamma0;
    double dGamma;
    int obscode;
    Pset parameters;
   
    {
      const int def=PsetMember::hasDefault;
      const int low=PsetMember::hasLowerBound;
      const int up=PsetMember::hasUpperBound;
      const int lowopen = low | PsetMember::openLowerBound;
      const int upopen = up | PsetMember::openUpperBound;

      parameters.addMember("observationFile",&observationPath, def,
			   "FITS file holding observation list (null=>stdin)", "");
      parameters.addMember("orbitFile",&orbitPath, def,
			   "FITS file for input orbit fits", "");
      parameters.addMember("predictFile",&predictPath, def,
			   "FITS file for output orbit predictions (null=>stdout)", "");
      parameters.addMember("ephemerisFile",&ephemerisPath, def,
			   "SPICE file (null=>environment", "");
      parameters.addMember("exposureFile",&exposurePath, def,
			   "FITS file holding DECam exposure info", "");
      parameters.addMember("obscode",&obscode, def | low,
			   "Observatory code (default: CTIO)", 807, 0);
    }
    parameters.setDefault();
    if (argc<2 || string(argv[1])=="-h" || string(argv[1])=="--help") {
      cout << usage << endl;
      cout << "\nParameters:" << endl;
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

    // Read the orbit data and build an index
    if (orbitPath.empty()) {
      cerr << "Must specify orbitFile" << endl;
      exit(1);
    }
    img::FTable orbitTable;
    try {
      FITS::FitsTable ft(orbitPath);
      orbitTable = ft.extract();
    } catch (std::runtime_error& m) {
      cerr << "Failure opening orbit file " << orbitPath << endl;
      quit(m);
    }

    map<int,int> orbitIndex;
    {
      vector<int> id;
      orbitTable.readCells(id, "ORBITID");
      for (int i=0; i<id.size(); i++)
	orbitIndex[id[i]]=i;
    }
    
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
    
    // Open the input file if we have it, set up output arrays for table if 
    // if we are using FITS output.
    bool observationsFromFile = !observationPath.empty();
    bool predictToFile = !predictPath.empty();
    bool useExpnum; // Set true if input will have expnum instead of MJD.

    // Input that is also output:
    vector<int> orbitID;
    vector<int> expnum;
    vector<double> mjd;  // Might use this instead of expnum
    
    if (observationsFromFile) {
      FITS::FitsTable ft(observationPath);
      img::FTable obsTable = ft.extract();
      obsTable.readCells(orbitID, "ORBITID");
      // Are we using expnum or MJD?
      auto colnames = obsTable.listColumns();
      useExpnum = (std::find(colnames.begin(),colnames.end(),"EXPNUM")!=colnames.end());
      if (useExpnum) {
	obsTable.readCells(expnum, "EXPNUM");
      } else {
	// We'll use MJD's
	obsTable.readCells(mjd,"MJD");
      }
    } else {
      // Reading observation info from stdin.
      // Let's read in all the input at the start to simplify code
      string buffer;
      int idThis, expnumThis;
      double mjdThis;
      // Let's read one line to decide whether we are expnum or mjd people:
      stringstuff::getlineNoComment(cin, buffer);
      {
	std::istringstream iss(buffer);
	iss >> idThis >> mjdThis;
      }
      if (mjdThis > 100000.) {
	useExpnum = true;
	expnum.push_back( static_cast<int> (round(mjdThis)));
      } else {
	useExpnum = false;
	mjd.push_back(mjdThis);
      }
      orbitID.push_back(idThis);

      // Now read everything
      while (stringstuff::getlineNoComment(cin, buffer)) {
	std::istringstream iss(buffer);
	if (useExpnum) {
	  iss >> idThis >> expnumThis;
	  orbitID.push_back(idThis);
	  expnum.push_back(expnumThis);
	} else {
	  iss >> idThis >> mjdThis;
	  orbitID.push_back(idThis);
	  mjd.push_back(mjdThis);
	}
      }
    }  // Done getting inputs

    // Load the exposure table if needed
    ExposureTable* exposureTable = useExpnum? new ExposureTable(exposurePath) : nullptr;

    // Create vectors we'll add to output
    int nobs = predictToFile ? orbitID.size() : 0; // No need for space if using stdout
    vector<int> flags(nobs,0.);
    vector<double> ra(nobs,0.);
    vector<double> dec(nobs,0.);
    vector<double> semimajor(nobs,0.);
      
    // Begin reading input data.  One Fitter should do for everything.
    Fitter fit(eph, Gravity::GIANTS);
    fit.setFrame(frame);
    
    int inRow = 0;
    while (inRow < orbitID.size()) {

      // Do all of the observations of a given orbitID at once
      int beginRow = inRow;
      int idThis = orbitID[beginRow];
      int endRow = beginRow;
      for ( ; endRow < orbitID.size() && orbitID[endRow]==idThis; endRow++) {}

      // Set up Fitter with ABG for this orbitID
      int errorThis=0;
      int index;
      if (!orbitIndex.count(idThis)) {
	// No orbit found.  Set flags, move on.
	errorThis = NO_ORBIT;
      } else {
	index = orbitIndex[idThis];
	orbitTable.readCell(errorThis, "FLAGS", index);
      }
      if (errorThis > 0) {
	// No prediction to make.
	for (int i = beginRow; i<endRow; i++) {
	  if (predictToFile) {
	    flags[i] = errorThis;
	  } else if (useExpnum) {
	    // Stdout, expnum format
	    cout << orbitID[i]
		 << " " << expnum[i]
		 << " " << errorThis
		 << " " << 0. << " " << 0. << " " << 0. << endl;
	  } else {
	    // Stdout, MJD format
	    cout << orbitID[i]
		 << " " << mjd[i]
		 << " " << errorThis
		 << " " << 0. << " " << 0. << " " << 0. << endl;
	  }
	} // end loop over bad orbit's observations
	inRow = endRow;
	continue;
      }

      // We have a legit orbit.  Load it into Fitter.
      vector<double> v;
      orbitTable.readCell(v,"ABG",index);
      ABG abg;
      for (int i=0; i<6; i++) abg[i] = v[i];
      orbitTable.readCell(v,"ABGCOV",index);
      ABGCovariance cov;
      for (int i=0; i<6; i++)
	for (int j=0; j<6; j++)
	  cov(i,j) = v[6*i+j];
      fit.setABG(abg,cov);
      
      // Make arrays for observing time and location in our frame
      DVector tobs(endRow-beginRow);
      DMatrix earth(3,endRow-beginRow);
      double mjdThis;
      astrometry::CartesianICRS xyz;
      for (int i=beginRow; i<endRow; i++) {
	if (useExpnum) {
	  if (!exposureTable->observingInfo(expnum[i],mjdThis,xyz)) {
	    cerr << "ERROR: Missing information for exposure " << expnum[i] << endl;
	    exit(1);
	  }
	  tobs[i-beginRow] = eph.mjd2tdb(mjdThis);
	} else {
	  double tdb = eph.mjd2tdb(mjd[i]);
	  xyz = eph.observatory(obscode, tdb);
	  tobs[i - beginRow] = tdb;
	}
	earth.col(i-beginRow) = xyz.getVector();
      }
      // Put data into frame
      tobs.array() -= frame.tdb0;
      DMatrix xE = frame.fromICRS(earth).transpose(); // Fitter wants Nx3, Frame is Nx3
      // Destinations for prediction data
      DVector x(endRow-beginRow);
      DVector y(endRow-beginRow);
      DVector covXX(endRow-beginRow);
      DVector covYY(endRow-beginRow);
      DVector covXY(endRow-beginRow);
      fit.predict(tobs, xE, &x, &y, &covXX, &covYY, &covXY);
      
      // Get major axes
      DVector tr2 = 0.5*(covXX + covYY);
      DVector det = covXX.array()*covYY.array() - covXY.array()*covXY.array();
      DVector a = sqrt(tr2.array() + sqrt(tr2.array()*tr2.array() - det.array()));
      // Convert to arcsec
      a /= ARCSEC;

      // Now stuff result back into array or write to screen
      for (int i=beginRow; i<endRow; i++) {
	int ii = i - beginRow;
	double raThis, decThis;
	astrometry::Gnomonic gn(x[ii],y[ii],frame.orient);
	astrometry::SphericalICRS icrs(gn);
	icrs.getLonLat(raThis, decThis);
	raThis /= DEGREE;
	decThis /= DEGREE;

	if (predictToFile) {
	  flags[i] = errorThis;
	  ra[i] = raThis;
	  dec[i] = decThis;
	  semimajor[i] = a[ii];
	} else {
	  cout << setw(6) << idThis;
	  if (useExpnum)
	    cout << " " << setw(6) << expnum[i];
	  else
	    cout << " " << fixed << setprecision(5) << setw(11) <<mjd[i];
	  cout << " " << setw(2) << errorThis
	       << " " << fixed << noshowpos << setprecision(6) << setw(10) << (raThis<0 ? raThis+360. : raThis)
	       << " " << showpos << setw(10) << decThis
	       << " " << noshowpos << setprecision(3) << a[ii]
	       << endl;
	}
      }

      // Advance row pointers and move to next orbit
      inRow = endRow;
    }
    
    // Write the FITS output file if needed
    if (predictToFile) {
      // Write the output table
      FITS::FitsTable ft(predictPath, FITS::Create + FITS::OverwriteFile);
      img::FTable out = ft.use();
      out.addColumn(orbitID,"ORBITID");
      if (useExpnum)
	out.addColumn(expnum,"EXPNUM");
      else
	out.addColumn(mjd,"MJD");
      out.addColumn(flags,"FLAGS");
      out.addColumn(ra,"RA");
      out.addColumn(dec,"DEC");
      out.addColumn(semimajor,"ERROR_A");
    }
  } catch (std::runtime_error& e) {
    quit(e);
  }
  exit(0);
}
