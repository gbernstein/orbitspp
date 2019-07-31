// Run orbit fits given sets of observations stored in FITS tables
#include "Pset.h"
#include "OrbitTypes.h"
#include "Elements.h"
#include "Fitter.h"
#include "FTable.h"
#include "FitsTable.h"
#include "Astrometry.h"
#include "Legacy.h"

#include <iostream>

const string usage =
  "FitsToABG: Read state vectors and uncertainties in ABG basis from a FITS\n"
  "           file and output for each object a file in the old .abg ASCII\n"
  "           format plus some heliocentric orbital elements in the old aei format.\n"
  "Usage:\n"
  "  FitsToABG <input FITS file> [-param value ...]\n"
  "  Program options are listed below.\n"
  "  The input epoch will be defaulted to reference time.  Values can either be a year\n"
  "  or an MJD, or a JD.\n"
  " ";

using namespace std;
using namespace orbits;
using namespace astrometry;


int main(int argc,
	 char *argv[])
{
  try {
    string ephemerisPath;
    double outEpoch;
    Pset parameters;

    {
      const int def=PsetMember::hasDefault;
      const int low=PsetMember::hasLowerBound;
      const int up=PsetMember::hasUpperBound;
      const int lowopen = low | PsetMember::openLowerBound;
      const int upopen = up | PsetMember::openUpperBound;

      parameters.addMember("outEpoch",&outEpoch, def,
			   "epoch for output elements", 0.);
      parameters.addMember("ephemerisFile",&ephemerisPath, def,
			   "SPICE file (null=>environment", "");
    }
    parameters.setDefault();

    if (argc<2 || string(argv[1])=="-h" || string(argv[1])=="--help") {
      cout << usage << endl;
      parameters.dump(cerr);
      exit(1);
    }
    
    parameters.setFromArguments(argc, argv);
    
    string inFileName = argv[1];
    
    // And read arguments from the remaining command line entries
    parameters.setFromArguments(argc, argv);

    // Read the ephemeris
    Ephemeris eph(ephemerisPath);

    // If the 
    bool useReferenceEpoch = false;
    if (outEpoch <= 0.) {
      useReferenceEpoch = true;
    } else if (outEpoch < 3000.) {
      // This is an epoch; subtract 2000 to give our tdb
      outEpoch -= 2000.;
    } else if (outEpoch<2400000.) {
      // This looks like MJD, convert to tdb
      outEpoch = eph.mjd2tdb(outEpoch);
    } else {
      // Looks like a JD
      outEpoch = eph.jd2tdb(outEpoch);
    }
   
    // Acquire input data
    FITS::FitsTable ft(inFileName,FITS::ReadOnly,1);
    auto inputTable = ft.use();

    for (int row=0; row < inputTable.nrows(); row++) {
      vector<double> vframe;
      inputTable.readCell(vframe, "frame", row);
      vector<double> vabg;
      inputTable.readCell(vabg, "ABG", row);
      vector<double> vabgcov;
      inputTable.readCell(vabgcov, "ABGINVCOV", row);
      
      // Reconstruct structures from the tabular input:
      
      // Convert frame to ReferenceFrame - unpack what MergeOrbits does
      double ra = vframe[0];
      double dec = vframe[1];
      double pa = vframe[2];
      double originX = vframe[3];
      double originY = vframe[4];
      double originZ = vframe[5];
      double tdb0 = vframe[6];

      if (useReferenceEpoch)
	outEpoch = tdb0;
      Orientation orient(SphericalICRS(ra,dec), pa);
      CartesianICRS origin(originX,originY,originZ);
      Frame frame(origin, orient, tdb0);
      ABG abg;
      for (int i=0; i<6; i++)
	abg[i] = vabg[i];
      ABGCovariance abgcov;
      {
	Matrix66 invcov;
	int k=0;
	for (int i=0; i<6; i++)
	  for (int j=0; j<6; j++, k++)
	    invcov(i,j) = vabgcov[k];
	abgcov = invcov.inverse();
      }
      
      // Names of output files
      ostringstream oss;
      oss << "obj" << setfill('0') << setw(4) << row;
      string aeiName = oss.str() + ".aei";
      string abgName = oss.str() + ".abg";
      
      // Open and write ABG file - do this using same C I/O as in the orbits software
      writeOldABG(abgName, abg, abgcov, frame, eph);

      // The Fitter class knows how to propagate ABG and uncertainties
      Fitter fit(eph);
      fit.setFrame(frame);
      fit.setABG(abg,abgcov);

      // Get Heliocentric elements and uncertainty
      if (useReferenceEpoch) {
	// Do not need to give a tdb, using reference epoch
	auto el = fit.getElements(true);
	auto elCov = fit.getElementCovariance(true);
	/* Print out the results, with comments */
	writeOldAEI(aeiName, el, elCov, outEpoch, eph);
      } else {
	auto el = fit.getElements(outEpoch,true);
	auto elCov = fit.getElementCovariance(outEpoch,true);
	/* Print out the results, with comments */
	writeOldAEI(aeiName, el, elCov, outEpoch, eph);
      }
    }

  } catch (std::runtime_error& e) {
    quit(e);
  }
  exit(0);
}
