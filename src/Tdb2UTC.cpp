#include <vector>
#include <iostream>

#include "Std.h"
#include "FitsTable.h"
#include "Ephemeris.h"

using namespace std;

string usage =
  "Print UTC time values corresponding to a column of TDB's\n"
  "in a FITS table\n"
  "usage: Tdb2UTC <filename> <column>";

int
main(int argc, char *argv[]) {
  if (argc != 3) {
    cerr << usage << endl;
    exit(1);
  }

  string filename = argv[1];
  string dataColumn = argv[2];

  try {
    orbits::Ephemeris ephem;
    astrometry::UT ut;
    FITS::FitsTable ff(filename,FITS::ReadWrite,1);
    auto table = ff.use();
    vector<double> v;
    vector<LONGLONG> id;
    int y, m;
    double d;
    table.readCells(v,dataColumn);
    table.readCells(id,"ORBITID");
    cout.precision(6);
    cout.fill('0');
    cout.setf( std::ios::fixed, std:: ios::floatfield );
    for (int i=0; i<v.size(); i++) {
      ut.setMJD(ephem.tdb2mjd(v[i]));
      ut.getYMD(y,m,d);
      cout << setw(3) << std::right << id[i] << " "
	   << setw(4) << std::right << y
	   << " " << setw(2) << std::right << m
	   << " " << setw(cout.precision()+3) << d
	   << endl;
    }
    table.writeCells(v,dataColumn);
    
  } catch (std::runtime_error& m) {
    quit(m);
  }

  exit(0);
}
