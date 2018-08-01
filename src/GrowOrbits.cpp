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

const string usage =
  "GrowOrbits: Given lists of transient detections potentially corresponding to the same\n"
  "object, looks through the transient catalog added possible additional detections\n"
  "one by one until nothing new fits the orbit.";

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
  vector<Exposure*> allExposures;
  vector<double> reducedChisqToKeep; // Definition of valid orbit, indexed by # detections
  bool goodFit(Fitter* fitptr) {
    // Criterion we'll use for whether an orbit fit is good enough to pursue:
    double chisqPerDetection = fitptr->getChisq() / (fitptr->nObservations());
    double threshold = (fitptr->nObservations()< reducedChisqToKeep.size()) ?
      reducedChisqToKeep[fitptr->nObservations()] :
      reducedChisqToKeep.back();
    return chisqPerDetection <= threshold;
  }
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
      // Exclude points that are tagged as invalid or bad matches to orbit.
      if (chisq[i] > SINGLE_CHISQ_THRESHOLD || !info.eptr->valid[i])
	continue;
      auto nextFit = fitptr->augmentObservation(info.eptr->tobs,
						info.eptr->xy(i,0), info.eptr->xy(i,1),
						info.eptr->covXX[i], info.eptr->covYY[i], info.eptr->covXY[i],
						info.eptr->earth,
						true);  // recalculate gravity ???
      // Check chisq to see if we keep it
      if (!globals.goodFit(nextFit))
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
  string tripletPath;
  string ephemerisPath;
  string transientPath;
  string exposurePath;
  double bindingFactor;
  double gamma0;
  double dGamma;
  double searchRadius;
  try {
    Pset parameters;
   
    {
      const int def=PsetMember::hasDefault;
      const int low=PsetMember::hasLowerBound;
      const int up=PsetMember::hasUpperBound;
      const int lowopen = low | PsetMember::openLowerBound;
      const int upopen = up | PsetMember::openUpperBound;

      parameters.addMember("tripletFile",&tripletPath, def,
			   "FITS file holding initial triplets (null->ASCII stdin)", "");
      parameters.addMember("ephemerisFile",&ephemerisPath, def,
			   "SPICE file (null=>environment", "");
      parameters.addMember("exposureFile",&exposurePath, def,
			   "FITS file holding DECam exposure info", "");
      parameters.addMember("transientFile",&transientPath, def,
			   "FITS file holding transient detections", "");
      parameters.addMember("bindingFactor",&bindingFactor, def | low,
			   "Chisq penalty at unbinding", 4., 0.);
      parameters.addMember("gamma0",&gamma0, def | lowopen,
			   "Center of gamma search region", 0.03, 0.);
      parameters.addMember("dGamma",&dGamma, def | lowopen,
			   "Half-width of gamma search region", 0.01, 0.);
      parameters.addMember("searchRadius",&searchRadius, def | lowopen,
			   "Radius of alpha/beta values at reference time", 1., 0.);
    }
    parameters.setDefault();

    if (argc<2 || string(argv[1])=="-h" || string(argv[1])=="--help") {
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

    // Read the ephemeris
    Ephemeris eph(ephemerisPath);

    // Open inputs
    // ??? read triplets

    // ??? Set up the chisq thresholds in globals
    
    // set up reference frame, pick gamma range
    Frame frame;
    {
      // Hardwire for zone029: ???
      double ra0 = 26.;
      double dec0 = -5.;
      frame.tdb0 = 15.9;  // Late Nov, Y3
      astrometry::SphericalICRS pole(ra0*DEGREE, dec0*DEGREE);
      astrometry::Orientation orient(pole);
      orient.alignToEcliptic();  // We'll do this by default.
      frame.orient = orient;
      // Put origin at position of observatory at MJD0
      frame.origin = eph.observatory(807, frame.tdb0);
    }    

    // Pluck out all relevant exposure data.  Then close
    // transient catalog and exposure catalog.
    auto exposurePool = selectExposures(frame,
					eph,
					gamma0,
					dgamma,
					searchRadius,
					transientPath,
					exposurePath);

    // Make an initial exposure list of all of these for the input orbits
    list<Exposure*> allExposures;
    allExposures.insert(allExposures.end(), exposurePool.begin(), exposurePool.end());
    
    // Make a lookup table from object ID back to the exposure it
    // was taken in.  This keeps us from having to keep the
    // transient file open.
    std::map<int, Detection> detectionIndex;
    for (auto& eptr : exposurePool) {
      for (int i=0; i<eptr->id.size(); i++)
	detectionIndex[eptr->id[i]] = Detection(eptr,i);
    }

    int orbitID = 0;

    // Now enter a loop of growing
    list<FitStep*> fitQueue;

    do {
      if (fitQueue.empty()) {
	// Stock the queue with fresh orbit from input
	string buffer;
	/* Loop over all initial detection sets: ?? simple stdin format to start */
	while (getlineNoComment(cin, buffer)) {
	  istringstream iss(buffer);
	  // Read object id's on a line
	  vector<Detection> members;
	  int objectID;
	  while (iss >> objectID) {
	    if (!detectionIndex.count(objectID)) {
	      cerr << "Object with ID " << objectID
		   << " is not in one of our search exposures"
		   << endl;
	      continue;
	    }
	    members.push_back(detectionIndex[objectID]);
	  }
      
	  if (members.size() < 3)
	    continue; // don't try these...

	  // Make an orbit fit
	  auto fitptr = new Fitter(eph, Gravity::BARY);  // ?? Need giants' gravity??
	  fitptr->setFrame(frame);
	  if (bindingFactor > 0.)
	    fitptr->setBindingConstraint(bindingFactor);
	  fitptr->setGammaPrior(gamma0, dGamma/2.); // Give gamma prior so search bound is 2 sigma
	  // Pass in the detections' data
	  {
	    int nobs = members.size();
	    DVector tobs(nobs);
	    DVector thetaX(nobs);
	    DVector thetaY(nobs);
	    DVector covxx(nobs);
	    DVector covyy(nobs);
	    DVector covxy(nobs);
	    DMatrix xE(nobs,3);
	    for (int i=0; i<members.size(); i++) {
	      Detection& det = members[i];
	      tobs[i] = det.eptr->tobs;
	      thetaX[i] = det.eptr->xy(i,0);
	      thetaY[i] = det.eptr->xy(i,1);
	      covxx[i] = det.eptr->covXX[i];
	      covyy[i] = det.eptr->covYY[i];
	      covxy[i] = det.eptr->covXY[i];
	      xE.row(i) = det.eptr->earth.transpose();
	    }
	    fitptr->setObservationsInFrame(tobs, thetaX, thetaY, covxx, covyy, covxy, xE);
	  }
	  bool goodFit = true;
	  try {
	    fitptr->setLinearOrbit();
	    fit.newtonFit();  // ?? any more elaborate fitting needed ??
	  } catch (std::runtime_error& e) {
	    goodFit = false;
	  }

	  // Skip the orbit if fit failed or is poor
	  if (!goodFit || !globals.goodFit(fitptr)) {
	    delete fitptr;
	    continue;
	  }

	  // Fit is good - queue up a FitStep
	  // First make a list of possible exposures that excludes the
	  // ones that the original detections come from.
	  list<Exposure*> possibleExposures(allPossibleExposures);
	  for (auto& m : members) {
	    possibleExposures.remove(m.eptr);
	  }
	  // Now queue up
	  fitQueue.push_back(new FitStep(fitptr, possibleExposures, members,
					 members.size(), members, orbitID,
					 1000.)); // Set a high FPR so we don't trigger on this
	  ++orbitID; // Increment orbit id
	} // End of input-reading loop.

	if (fitQueue.empty()) {
	  // Nothing left to fit!
	  exit(0);
	}

      } // End of searching for new inputs for queue.

	// Search the next item on the queue:
      auto newFits = fitQueue.front()->search();
      if (newFits.empty()) {
	// No new points can be added to this fit.  So it's a completed search.
	// ??? Have some criterion for output ???
	const FitStep& fs = *fitQueue.front();
	cout << "Orbit " << fs.orbitID
	     << " " << fs.nIndependent
	     << " " << fs.fptr->getChisq() << " / " << fs.fptr->getDOF()
	     << " FPR " << fs.fpr
	     << endl;
	fs.fptr->getABG().write(cout);

	// ?? Better criterion for when a fit is secure
	if (fs.nIndependent > 5) {
	  // Take this fit's detections out of circulation for future fits
	  for (auto& m : fs.members) {
	    m.eptr->valid[m.objectRow] = false;
	  }
	}
      } else {
	// Not a terminal fit - add new searches to queue
	fitQueue.insert(fitQueue.end(), newFits.begin(), newFits.end());
      }
      // Done with this fit
      delete fitQueue.front();
      fitQueue.pop_front();
    } // End of the queue loop.

    cerr << "ERROR: Should not get here." << endl;
  } catch (std::runtime_error& e) {
    quit(e,1);
  }
}  
