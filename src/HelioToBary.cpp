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
  "HelioToBary: convert orbital elements from MPC heliocentric to barycentric\n"
  "             in same format.\n"
  "Usage:\n"
  "  BulkElements <input file> <output file> [-param value ...]\n"
  "  Program options are listed below.\n"
  "\n"
  "Input and output files have degrees for angles.";

using namespace std;
using namespace orbits;

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
			   "Epoch for output elements (e.g. 2000., default is retain input)",-1.);
      parameters.addMember("ephemerisFile",&ephemerisPath, def,
			   "SPICE file (null=>environment", "");
    }
    parameters.setDefault();

    if (argc<3 || string(argv[1])=="-h" || string(argv[1])=="--help") {
      cout << usage << endl;
      parameters.dump(cerr);
      exit(1);
    }
    
    parameters.setFromArguments(argc, argv);
    
    string inFileName = argv[1];
    string outFileName = argv[2];
    
    // And read arguments from the remaining command line entries
    parameters.setFromArguments(argc, argv);

    // Put output epoch into years since J2000
    double outTDB = outEpoch > 0. ? outEpoch - 2000. : 0.;

    // Read the ephemeris
    Ephemeris eph(ephemerisPath);

    // Acquire input data
    FITS::FitsTable ft(inFileName,FITS::ReadOnly,1);
    auto inputTable = ft.use();

    vector<double> a;
    inputTable.readCells(a, "a");
    vector<double> e;
    inputTable.readCells(e, "e");
    vector<double> incl;
    inputTable.readCells(incl, "i");
    vector<double> tPeri;
    inputTable.readCells(tPeri, "Tp");
    vector<double> epoch;
    inputTable.readCells(epoch, "Epoch");
    vector<double> lan;
    inputTable.readCells(lan, "Node");
    vector<double> aop;
    inputTable.readCells(aop, "Peri");
    vector<string> name;
    inputTable.readCells(name, "Name");
    vector<string> number;
    inputTable.readCells(number, "Number");

    // Row-by-row, change elements to barycentric and overwrite
    for (int row=0; row<a.size(); row++) {
      Elements el;
      el[Elements::A] = a[row];
      el[Elements::E] = e[row];
      el[Elements::I] = incl[row] * DEGREE;
      el[Elements::LAN] = lan[row] * DEGREE;
      el[Elements::AOP] = aop[row] * DEGREE;
      el[Elements::TOP] = eph.jd2tdb(tPeri[row]);

      double inEpoch = eph.jd2tdb(epoch[row]);
      
      // Get state vector, assuming helio orbits
      auto s = getState(el, inEpoch, true, &eph);

      // If we want the output epoch to differ, we must integrate the
      // orbit to a new state.
      if (outEpoch > 0.) {
	Trajectory t(eph, s, orbits::GIANTS, 0.3*DAY);
	astrometry::CartesianICRS v;
	astrometry::CartesianICRS x = t.position(outTDB, &v);
	s.x = x;
	s.v = v;
	s.tdb = outTDB;
      }

      el = getElements(s);
	  
      // Now get back Bary elements
      a[row] = 	   el[Elements::A];
      e[row] = 	   el[Elements::E];
      incl[row] =  el[Elements::I] / DEGREE;
      lan[row]  =  el[Elements::LAN] / DEGREE;
      aop[row]  =  el[Elements::AOP] / DEGREE;
      tPeri[row]=  eph.tdb2jd(el[Elements::TOP]);

      if (outEpoch > 0.) {
	epoch[row] = eph.tdb2jd(outTDB);
      }
    } // end row loop

    // Make output table
    {
      FITS::FitsTable ff(outFileName, FITS::Create + FITS::OverwriteFile);
      img::FTable outTab = ff.use();
      outTab.addColumn(name,"Name",-1,14);
      outTab.addColumn(number,"Number",-1,8);
      outTab.addColumn(a,"a");
      outTab.addColumn(e,"e");
      outTab.addColumn(incl,"i");
      outTab.addColumn(lan,"Node");
      outTab.addColumn(aop,"Peri");
      outTab.addColumn(tPeri,"Tp");
      outTab.addColumn(epoch,"Epoch");
    }

  } catch (std::runtime_error& e) {
    quit(e);
  }
  exit(0);
}
