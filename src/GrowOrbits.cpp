// Find all additional detections that are consistent with orbit fitting input detections.
#include <list>
#include <vector>
#include "AstronomicalConstants.h"
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

const double FIELD_RADIUS = 1.1*DEGREE;  // Generous DECam field radius
// Chisq small enough to consider old orbit matched to new detection:
const double SINGLE_CHISQ_THRESHOLD = 16.;  
const double NOMINAL_POSITION_ERROR = 0.1*ARCSEC; // Position error for typical detections
// Area to be used for assessing false positive rate:
const double FALSE_POSITIVE_COUNTING_AREA = 0.5*DEGREE*DEGREE;
// We won't attempt to link to exposures that have too many
// expected detections, defined as being above this number
const double MAXIMUM_FPR_PER_EXPOSURE=10.;
// and at least this many independent detections so far:
const int START_IGNORING_HIGH_FPR=4;
// Time that must pass between exposures to be considered independent detections
// (e.g. when asteroids or defects would have moved out of linking range)
const double INDEPENDENT_TIME_INTERVAL = 0.1*DAY;

struct GlobalResources {
  Ephemeris* ephem;
  vector<Exposure*> allExposures;
  vector<double> reducedChisqToKeep; // Definition of valid orbit, indexed by # detections
} globals;

struct Detection {
  Exposure* eptr;
  int objectRow;
  Detection(Exposure* eptr_, int objectRow_): eptr(eptr_), objectRow(objectRow_) {}
};
  
class FitStep {
  // Class that holds an orbit fit to N objects, plus a vector of pointers
  // to exposures that could hold the (N+1)th.
public:
  FitStep(Fitter* fitptr_,
	  vector<Exposure*> possibleExposures_,
	  vector<Detection> members_,
	  int nIndependent_,
	  int orbitID_,
	  double fpr_):
    fitptr(fitptr_), possibleExposures(possibleExposures_),
    members(members_), nIndependent(nIndependent_),
    orbitID(orbitID), fpr(fpr_) {}
  ~FitStep() {delete fitptr;}

  const Fitter* getFitter() const {return fitptr;}

  list<FitStep*> search(); // return all possible N+1 -object FitSteps.
private:
  Fitter* fitptr;  // An orbit fit to the N objects

  // Exposures potentially holding the (N+1)th.
  // This list should not include exposures already having a detection, but might
  // include some v close in time.
  vector<Exposure*> possibleExposures;
  
  vector<Detection> members;  // Detections that comprise the orbit
  
  // Number of detections (so far) with independent asteroids
  // (i.e. don't count consecutive exposures as >1)
  int nIndependent;  

  // A unique ID for each starting orbit
  int orbitID;
  
  double fpr;  // False positive rate for having found this guy.
};
    
// This method of the FitStep is the key ingredient:
list<FitStep*>
FitStep::search() {
  list<FitStep*> out;

  struct ExposureInfo {
    Exposure* eptr;
    double x,y;  // Position of this fit at exposure's time
    double covXX, covYY, covXY;  // and covariance
  };

  // Narrow the list of exposures with potential links
  list<ExposureInfo> infoList;

  {
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

    auto eptr = possibleExposures.begin();
    for (int i=0;
	 eptr!=possibleExposures.end();
	 ++i, ++eptr) {
      tobs[i] = (*eptr)->tobs;
      xE.row(i) = (*eptr)->earth.transpose();
      axes.row(i) = (*eptr)->axis.transpose();
    }

    fitptr->predict(tobs, xE, &x, &y, &covXX, &covYY, &covXY);
    // Get major axes
    DVector tr2 = 0.5*(covXX + covYY);
    DVector det = covXX.array()*covYY.array() - covXY.array()*covXY.array();
    DVector a = sqrt(tr2.array() + sqrt(tr2.array()*tr2.array() - det.array()));
    // Add field radius to look for matches
    a.array() += FIELD_RADIUS;
    axes.col(0) -= x;
    axes.col(1) -= y;
    BVector useExposure = (axes.array()*axes.array()).rowwise().sum() <= a.array()*a.array();
    
    // Collect info for possible exposures, and
    // pop the ones that don't cross orbit
    auto eiter=possibleExposures.begin();
    for (int i=0;
	 eiter!=possibleExposures.end();
	 i++) {
      if (useExposure[i]) {
	ExposureInfo info;
	info.eptr = *eiter;
	info.x = x[i];  info.y = y[i];
	info.covXX=covXX[i]; info.covYY=covYY[i]; info.covXY=covXY[i];
	infoList.push_back(info);
	++eiter;
      } else {
	eiter = possibleExposures.erase(eiter);
      }
    }
  } // Done making exposure list.
  
  // Now walk through all viable exposures, collecting
  // information on the actual matches and the
  // false positive expectation.
  
  double newFPR = 0.; // False positive rate that will be handed to children
  list<Detection> matches;  // List of N+1's.

  for (auto& info : infoList) {
    // Determine the search ellipse area on this exposure
    // with some rough allowance for measurement errors
    double searchArea = PI * sqrt( (info.covXX*info.covYY - info.covXY*info.covYY)
				   + NOMINAL_POSITION_ERROR*NOMINAL_POSITION_ERROR);
    // False positive rate will be calculated by scaling error ellipse
    // to have this size, and counting detections in it.
    double chisqForFalsePositive = FALSE_POSITIVE_COUNTING_AREA / searchArea;
    searchArea *= SINGLE_CHISQ_THRESHOLD;
    
    // Get chisq of all points
    DVector chisq = info.eptr->chisq(info.x,info.y,info.covXX,info.covYY,info.covXY);

    // ?? This would be where to mask detections for those already claimed.

    double thisFPR = (chisq.array() <= chisqForFalsePositive).count()
      * searchArea / FALSE_POSITIVE_COUNTING_AREA;
    if ((nIndependent>=START_IGNORING_HIGH_FPR)
	&& (thisFPR > MAXIMUM_FPR_PER_EXPOSURE)) {
      // This exposure is too busy for this stage of linking, do not link to it.
      continue;
    } else {
      // Add this exposure to FPR
      newFPR += thisFPR;
    }

    // See if this exposure is close in time to any of the detections
    // already in the orbit
    
    bool closeInTime = false;
    for (auto& m : members) {
      if (abs(info.eptr->tobs-m.eptr->tobs) < INDEPENDENT_TIME_INTERVAL) {
	closeInTime = true;
	break;
      }
    }
    
    // Now refit orbit to all potential candidates
    // Refit for all good points
    for (int i=0; i<chisq.size(); i++) {
      if (chisq[i] > SINGLE_CHISQ_THRESHOLD)
	continue;
      auto nextFit = fitptr->augmentObservation(info.eptr->tobs,
						info.eptr->xy(i,0), info.eptr->xy(i,1),
						info.eptr->covXX[i], info.eptr->covYY[i], info.eptr->covXY[i],
						info.eptr->earth,
						true);  // recalculate gravity ???
      // Check chisq to see if we keep it
      double chisqPerDetection = nextFit->getChisq() / (nextFit->nObservations());
      double threshold = (nextFit->nObservations()<globals.reducedChisqToKeep.size()) ?
	globals.reducedChisqToKeep[nextFit->nObservations()] :
	globals.reducedChisqToKeep.back();
      if (chisqPerDetection > threshold) {
	delete nextFit;
	continue;
      }

      // Make a new FitStep that has this detection
      out.push_back(new FitStep(nextFit, possibleExposures, members,
			       (closeInTime ? nIndependent : nIndependent+1),
			       orbitID, newFPR));
      // Add this detection to the new FitStep
      out.back()->members.push_back(Detection(info.eptr, i));
    }
  } // End of N+1 exposure loop
  return out;
}

  /**
      May wish to cull multiple fits that correspond to same input.
      ?? Is there a way to remove well-fit detections from the detection pool??
      Have a "claimed" bit for each detection in the master Exposure list, which is altered when an exposure/fit reaches "terminal" stage (no children) with an adequately good fit?  If multithreaded this would starve any current /queued / future fitters from grabbing the likely additional detections.  Hopefully one of them would find everything it needs first.
  **/

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
  */

   


}
