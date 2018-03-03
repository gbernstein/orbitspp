// Test MPCObservation
#include "OrbitTypes.h"
#include "StringStuff.h"
#include <iostream>
#include <iomanip>

using namespace std;
using namespace orbits;

int main(int argc,
	 char *argv[])
{
  string line;
  while (stringstuff::getlineNoComment(cin, line)) {
    try {
      MPCObservation obs(line);
      cout << "Position " << obs.radec
	   << " at jd " << fixed << setprecision(7) << obs.mjd
	   << " obsid " << obs.obscode
	   << endl;
    } catch (runtime_error& m) {
      cout << "Caught: " << m.what() << endl;
    }
  }
  exit(0);
}

   
