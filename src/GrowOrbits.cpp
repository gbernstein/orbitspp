// Find all additional detections that are consistent with orbit fitting input detections.
#include <list>
#include <vector>
#include "Astrometry.h"
#include "OrbitTypes.h"
#include "PlaneGeometry.h"  // ???
#include "Ephemeris.h"
#include "Fitter.h"
#include "Exposures.h"

using namespace std;
using namespace orbits;

/**

Parameters:
Full detection catalog
Initial orbit detection catalog

exposure catalog
ephemeris

output file

some parameters that will matter:
* binding factor
* chisq thresholds for refits
* false positive rate or # detections to keep for output
* max error region for N+1 exposure?

**/

struct GlobalResources {
  Ephemeris* ephem;
  vector<Exposure*> allExposures;
  vector<double> reducedChisqToKeep; // Definition of valid orbit, indexed by # detections
} globals;

class FitStep {
  // Class that holds an orbit fit to N objects, plus a vector of pointers
  // to exposures that could hold the (N+1)th.
public:
  FitStep(Fitter* fitptr_,
	  vector<Exposure*> possibleExposures_,
	  int nIndependent_)
    fitptr(fitptr_), possibleExposures(possibleExposures_), nIndependent(nIndependent_) {}
  ~FitStep() {delete fp;}

  const Fitter* getFitter() const {return fp;}

  list<FitStep*> search(); // return all possible N+1 -object FitSteps.
private:
  Fitter* fitptr;  // An orbit fit to the N objects

  // Exposures potentially holding the (N+1)th.
  // This list should not include exposures already having a detection, but might
  // include some v close in time.
  vector<Exposure*> possibleExposures;
  
  // Number of detections (so far) with independent asteroids
  // (i.e. don't count consecutive exposures as >1)
  int nIndependent;  

  double fpr;  // False positive rate
};
    
// This method of the FitStep is the key ingredient:
list<FitStep*>
FitStep::search() {
  list<FitStep*> out;

  // Get predictions from orbit for all potential exposures.
  int nExpo = possibleExposures.size();
  DVector tobs(nExpo);
  DMatrix xE(nExpo, 3); // Fitter wants Nx3
  DMatrix axes(nExpo, 2); // Exposure pointings
  // Destinations for prediction data
  DVector x(nExpo);
  DVector y(nExpo);
  DVector covXX(nExpo);
  DVector covYY(nExpo);
  DVector covXY(nExpo);

  for (int i=0, auto eptr = possibleExposures.begin();
       eptr!=possibleExposures.end();
       ++i, ++eptr) {
    tobs[i] = eptr->tobs;
    xE.row(i) = eptr->earth.transpose();
    axes.row(i) = eptr->axis.transpose();
  }

  fit.predict(tobs, xE, &x, &y, &covXX, &covYY, &covXY);
  // Get major axes
  DVector tr2 = 0.5*(covXX + covYY);
  DVector det = covXX.array()*covYY.array() - covXY.array()*covXY.array();
  DVector a = sqrt(tr2.array() + sqrt(tr2.array()*tr2.array() - det.array()));
  
}

int
main(int argc, char **argv) {

  // Read parameters

  // Open inputs

  // Create output arrays

  // Read ephemeris

  // set up reference frame, pick gamma range


  // Pluck out all relevant exposure data.  Then close
  // transient catalog and exposure catalog.


  /* Loop over all initial detection sets:

     Fit an orbit
     quit if not good enough.
     Assign an input ID

     Create an object which is a combination of fitted orbit to N objects
     and a catalog of exposures with possible (N+1)th detection.  
     And the input ID.  Add it to a queue?

     Loop through the queue.  For each object in the queue:
        Predict object position on all candidate exposures.
	Pop the impossible exposures out of the candidate list.
	For each other exposure:
	   Calculate chisq for all detections
	   Make a list of the ones that could be N+1.
	   Get FPR by counting number within chisq corresponding to CCD-sized area,
	     and rescaling to real search ellipse area.
	   (skip exposures with too-large error regions somehow?)
	   (treat new exposures within time ~2 hrs as not independent)
	Get total FPR.
	Now loop through potential N+1's:
	   Spawn a new Fitter.
	   Check its chisq - if no good, discard.
	   If good, and FPR/Ndetect threshold is met - save as output
	   Then make a new exposure/fit object, put on queue for next search.
	Destroy this exposure/fit

      May wish to cull multiple fits that correspond to same input.
      ?? Is there a way to remove well-fit detections from the detection pool??
      Have a "claimed" bit for each detection in the master Exposure list, which is altered when an exposure/fit reaches "terminal" stage (no children) with an adequately good fit?  If multithreaded this would starve any current /queued / future fitters from grabbing the likely additional detections.  Hopefully one of them would find everything it needs first.

  */

   


}
