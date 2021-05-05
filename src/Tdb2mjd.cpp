// Change T_p column (in tdb years) to UT JDB in
// objects table.

#include <vector>
#include <iostream>

#include "Std.h"
#include "FitsTable.h"
#include "Ephemeris.h"

using namespace std;

string usage =
  "Change columns T_p and sigma_T in table from being TDB years\n"
  "to MJD.\n"
  "usage: ExposurePrep <filename>";

int
main(int argc, char *argv[]) {
  if (argc != 2) {
    cerr << usage << endl;
    exit(1);
  }

  string filename = argv[1];
  string dataColumn = "T_p";
  string errorColumn = "sigma_T";

  try {
    orbits::Ephemeris ephem;
    FITS::FitsTable ff(filename,FITS::ReadWrite,1);
    auto table = ff.use();
    vector<double> v;
    table.readCells(v,dataColumn);
    for (int i=0; i<v.size(); i++) {
      // Routine expects tdb is yrs past 2000
      // v[i] = ephem.tdb2mjd(v[i]-2000.);
      v[i] = ephem.tdb2mjd(v[i]); // But auxmrt is already this
    }
    table.writeCells(v,dataColumn);
    v.clear();
    table.readCells(v,errorColumn);
    for (int i=0; i<v.size(); i++) {
      // 
      v[i] = v[i]/DAY;
    }
    table.writeCells(v,errorColumn);
    
  } catch (std::runtime_error& m) {
    quit(m);
  }

  exit(0);
}
