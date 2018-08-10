// Program to 
// ?? Check corners to get CCDNUM too?

#include <iostream>

#include "StringStuff.h"
#include "Astrometry.h"
#include "Fitter.h"
#include "Elements.h"
#include "Exposures.h"
#include "Pset.h"
#include "FitsTable.h"

using namespace std;
using namespace orbits;

const string usage =
  "DESTracks: gives positions of a fake orbit on every astrometric\n"
  "DES exposure that it crosses, and which CCD it lands on.\n"
  "The ra0, dec0, radius parameters give the sky regions from\n"
  "which exposures will be drawn (in degrees).  Restricted to\n"
  "radius<85 degrees to avoid singularities in gnomonic projection.\n"
  "Defaults will cover the full DES footprint.\n"
  "\n"
  "Usage: DESTracks [parameter file] [-parameter value]\n"
  "   [parameter file] is optional name(s) of file(s) holding parameter\n"
  "                name/value pairs, one per line\n"
  "   [-parameter value] are additional parameter name/value pairs.\n"
  "   Parameters are listed below.\n"
  "\n"
  "Inputs: stdin specifies one orbit per line, which is either a state\n"
  "        vector or orbital elements.  Blank lines or #comments are\n"
  "        ignored.  If using state vectors (readState=true, the default)\n"
  "        then each line is\n"
  "        <id> <x> <y> <z> <vx> <vy> <vz>\n"
  "        in ICRS barycentric coordinates, units of AU and Julian years.\n"
  "        <id> is an integer orbit ID, and the tdb0 parameter gives\n"
  "        epoch of the state.  If using orbital elements, each line gives\n"
  "        osculating barycentric elements at tdb0:\n"
  "        <id> <a> <e> <i> <LAN> <AOP> <TOP>\n"
  "        with angles in degrees and TOP given as TDB since J2000.\n"
  "Outputs: For each orbit given, one row per DES exposure crossed.\n"
  "        <id> <expnum> <tdb> <ra> <dec> <ccdnum>\n"
  "        which will be written to a FITS table if positionFile is\n"
  "        specified, otherwise written to stdout.  CCDNUM=0 if\n"
  "        orbit position is not in detectable region of a CCD with\n"
  "        valid astrometry - an exposure is listed if target is within\n"
  "        1.1 degrees of its optic axis.";

int main(int argc,
	 char *argv[])
{
  // Controlling parameters / constants
  string positionPath; // output file
  string ephemerisPath;
  string exposurePath;
  string cornerPath;  // File with CCD corners per exposure

  double ra0;    // Center of field (ignore exposures >60 degrees away)
  double dec0;
  double searchRadius;
  double tdb0;   // Reference time (=epoch of orbits or state vectors)
  bool readState;  // Input will be ICRS state vectors
  // Else orbital elements

  const int obscode=807;
  const double fieldRadius = 1.1*DEGREE;
  
  Pset parameters;
   
  try {
    {
      const int def=PsetMember::hasDefault;
      const int low=PsetMember::hasLowerBound;
      const int up=PsetMember::hasUpperBound;
      const int lowopen = low | PsetMember::openLowerBound;
      const int upopen = up | PsetMember::openUpperBound;

      parameters.addMember("positionFile",&positionPath, def,
			   "FITS file for output positions of observations (null=>ASCII/stdout)", "");
      parameters.addMember("ephemerisFile",&ephemerisPath, def,
			   "SPICE file (null=>environment", "");
      parameters.addMember("exposureFile",&exposurePath, def,
			   "DECam exposure info file (null=>environment)", "");
      parameters.addMember("cornerFile",&cornerPath, def,
			   "CCD corners table file (null=>environment)", "");
      parameters.addMember("ra0",&ra0, def,
			   "RA of field center (deg)", 10.);
      parameters.addMember("dec0",&dec0, def,
			   "Dec of field center (deg)", -20.);
      parameters.addMember("radius",&searchRadius, def || lowopen || up,
			   "Radius of search area (deg)", 85.,0.,85.);
      parameters.addMember("TDB0",&tdb0, def,
			   "TDB of reference time (=orbit epoch), yrs since J2000", 0.);
      parameters.addMember("readState",&readState, def,
			   "Input is ICRS state vector (T, default) or elements (F)", true);
    }
    parameters.setDefault();

    if (argc>1 && (string(argv[1])=="-h" || string(argv[1])=="--help")) {
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

    // Are we using a file or stdout for output?
    bool useStdout = positionPath.empty();

    // Read the ephemeris
    Ephemeris ephem(ephemerisPath);

    // Establish the reference frame
    Frame frame;
    {
      astrometry::SphericalICRS pole(ra0*DEGREE, dec0*DEGREE);
      frame.orient = astrometry::Orientation(pole);
      frame.orient.alignToEcliptic();
      frame.tdb0 = tdb0;
      frame.origin = ephem.observatory(obscode, tdb0);
    }

    // Read the exposure table
    ExposureTable et(exposurePath);

    // Get all exposures of relevance and build big arrays
    vector<Exposure*> exposures = et.getPool(frame, ephem, searchRadius*DEGREE);
    int nExposures = exposures.size();
    DVector tdb(nExposures);
    DMatrix earth(nExposures,3);  // ICRS observatory position
    DMatrix axis(nExposures,3);   // ICRS direction cosines of pointings
    for (int i=0; i<nExposures; i++) {
      tdb[i] = exposures[i]->tdb;
      earth.row(i) = exposures[i]->earthICRS.getVector().transpose();
      axis.row(i) = exposures[i]->axisICRS.getUnitVector().transpose();
    }

    cerr << "# Read " << exposures.size() << " exposures" << endl;
    
    // Load CCD corners
    {
      CornerTable tt(cornerPath);
      for (auto& eptr : exposures)
	tt.fillExposure(eptr);
    }

    // convert field radius to a maximum chord length between
    // unit vectors
    const double maxChordSq = pow(2. * sin(fieldRadius/2.),2.);

    // Set up output vectors (used only for table output)
    vector<LONGLONG> idOut;
    vector<int> expnumOut;
    vector<int> ccdnumOut;
    vector<double> tdbOut;
    vector<double> raOut;
    vector<double> decOut;
    

    // Begin reading input lines
    string buffer;
    while (stringstuff::getlineNoComment(cin, buffer)) {
      State xv;
      LONGLONG orbitID;
      xv.tdb = tdb0;
      istringstream iss(buffer);
      if (readState) {
	// Read state and orbitID
	iss >> orbitID
	    >> xv.x[0] >> xv.x[1] >> xv.x[2]
	    >> xv.v[0] >> xv.v[1] >> xv.v[2];
      } else {
	// Read elements and orbitID
	Elements el;
	iss >> orbitID >> el;
	// convert to state
	xv = getState(el, tdb0);
      }

      // Predict for all relevant exposures - get ICRS direction cosines
      Trajectory orbit(ephem, xv, Gravity::GIANTS);
      DMatrix target = orbit.observe(tdb, earth);
      // Find those within radius
      BVector hits = (target - axis).rowwise().squaredNorm().array() < maxChordSq;
      
      for (int i=0; i<nExposures; i++) {
	if (!hits[i]) continue;

	astrometry::SphericalICRS radec;
	radec.setUnitVector(target.row(i).transpose());
	double ra,dec;
	radec.getLonLat(ra,dec);

	if (useStdout) {
	  cout << setw(8) << orbitID
	       << " " << setw(6) << exposures[i]->expnum
	       << " " << fixed << setprecision(8) << setw(11) << exposures[i]->tdb
	       << " " << setprecision(7) << setw(11) << ra/DEGREE
	       << " " << showpos << setprecision(7) << setw(11) << dec/DEGREE
	       << " " << noshowpos << setw(2) << exposures[i]->whichCCD(radec)
	       << endl;
	} else {
	  idOut.push_back(orbitID);
	  expnumOut.push_back(exposures[i]->expnum);
	  ccdnumOut.push_back(exposures[i]->whichCCD(radec));
	  tdbOut.push_back(exposures[i]->tdb);
	  raOut.push_back(ra/DEGREE);
	  decOut.push_back(dec/DEGREE);
	}
      }
    } // End the input reading loop.

    if (!useStdout) {
      // Write the results to a FITS file
      FITS::FitsTable ft(positionPath, FITS::Create + FITS::OverwriteFile);
      img::FTable out = ft.use();
      out.addColumn(idOut,"ORBITID");
      out.addColumn(expnumOut,"EXPNUM");
      out.addColumn(ccdnumOut,"CCDNUM");
      out.addColumn(tdbOut,"TDB");
      out.addColumn(raOut,"RA");
      out.addColumn(decOut,"DEC");
    }
  } catch (std::runtime_error& e) {
    quit(e);
  }
  exit(0);
}

    

