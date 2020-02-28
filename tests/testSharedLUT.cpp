/* See if the SharedLUT does what we want */

#include <iostream>
#include "SharedLUT.h"

using namespace std;

int main(int argc, char *argv[]) {
  // Argument is size of vector to build

  int exitcode=0;
  SharedLUT<int> lut;
  
#pragma omp parallel for
  for (int i=0; i<atoi(argv[1]); i++) {
    // Each thread will fill the table 10 entries past its
    // number if the table does not already reach its number
    if (lut.size() <=i) {
      vector<int> v;
      int s=lut.size();
      for (int j=s; j<=i+10; j++)
	v.push_back(2*j);
#pragma omp critical(io)
      cerr << "Appending " << s << " to " << i+10 << endl;
      lut.append(s, v.begin(), v.end());
    }
    // Every loop step reads its value
    if (lut[i]!=2*i) {
      cerr << "Error in value at i=" << i << " " << lut[i] << endl;
      exitcode=1;
    }
  }

  exit(exitcode);
}

    
      
