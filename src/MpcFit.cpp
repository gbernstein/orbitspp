// Quick fit of orbit to MPC-format observations
#include "Pset.h"
#include "OrbitTypes.h"
#include "Elements.h"
#include "Fitter.h"
#include "Astrometry.h"
#include "Legacy.h"

#include <iostream>

const string usage =
  "MpcFit: Do an orbit fit from MPC-format file\n"
  "Usage:\n"
  "  FitsToABG <name> [-param value ...]\n"
  "    name - output files will be old-style <name>.abg, <name>.aei\n"
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
    double epoch;
    bool helio;  // True for heliocentric element outputs
    Pset parameters;

    {
      const int def=PsetMember::hasDefault;
      const int low=PsetMember::hasLowerBound;
      const int up=PsetMember::hasUpperBound;
      const int lowopen = low | PsetMember::openLowerBound;
      const int upopen = up | PsetMember::openUpperBound;

      parameters.addMember("epoch",&epoch, def,
			   "epoch for output elements", 0.);
      parameters.addMember("helio",&helio, def,
			   "true for heliocentric", false);
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
    
    string outFileName = argv[1];
    
    // And read arguments from the remaining command line entries
    parameters.setFromArguments(argc, argv);

    // Read the ephemeris
    Ephemeris eph(ephemerisPath);

    bool useReferenceEpoch = false;
    if (epoch <= 0.) {
      useReferenceEpoch = true;
    } else if (epoch < 3000.) {
      // This is an epoch; subtract 2000 to give our tdb
      epoch -= 2000.;
    } else if (epoch<2400000.) {
      // This looks like MJD, convert to tdb
      epoch = eph.mjd2tdb(epoch);
    } else {
      // Looks like a JD
      epoch = eph.jd2tdb(epoch);
    }

    Fitter fit(eph);
    fit.readMPCObservations(std::cin);

    fit.chooseFrame(-1);
    //**fit.chooseFrame(0);
    fit.setLinearOrbit();
    fit.newtonFit();
    fit.printResiduals(cerr);

    string aeiName = outFileName + ".aei";
    string abgName = outFileName + ".abg";
    
    // Open and write ABG file - do this using same C I/O as in the orbits software
    writeOldABG(abgName, fit.getABG(), fit.getInvCovarABG().inverse(), fit.getFrame(), eph);

    // Get elements and uncertainty
    if (useReferenceEpoch) {
      // Do not need to give a tdb, using reference epoch
      auto el = fit.getElements(helio);
      auto elCov = fit.getElementCovariance(helio);
      /* Print out the results, with comments */
      writeOldAEI(aeiName, el, elCov, epoch, eph);
    } else {
      auto el = fit.getElements(epoch,helio);
      auto elCov = fit.getElementCovariance(epoch,helio);
      /* Print out the results, with comments */
      writeOldAEI(aeiName, el, elCov, epoch, eph);      
    }
    Matrix66 sCov;
    State s = fit.predictState(epoch - fit.getFrame().tdb0, &sCov);
    cout << s << endl;
    cout << sCov << endl;
    
  } catch (std::runtime_error& e) {
    quit(e);
  }
  exit(0);
}

