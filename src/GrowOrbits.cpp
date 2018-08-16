// Find all additional detections that are consistent with orbit fitting input detections.
#include <list>
#include <vector>
#include <map>
#include <algorithm>
#include "AstronomicalConstants.h"
#include "Astrometry.h"
#include "OrbitTypes.h"
#include "Ephemeris.h"
#include "Fitter.h"
#include "Exposures.h"
#include "Pset.h"

using namespace std;
using namespace orbits;

////////////////////////////////////////////////////////////////////////////////////
// Constants / magic numbers 
////////////////////////////////////////////////////////////////////////////////////
const string usage =
  "GrowOrbits: Given lists of transient detections potentially corresponding to the same\n"
  "object, looks through the transient catalog added possible additional detections\n"
  "one by one until nothing new fits the orbit.";

const bool DEBUG=false;

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

// What maximum FPR to retain?
const int MIN_DETECTIONS_TO_OUTPUT = 4; // ???

// What combination of FPR and number of independent exposures are sufficient
// to consider this a completed search and flush all other FitSteps from the same
// starting orbit and pull the detections out of circulation.
const double MAX_FPR_EXCLUSIVE = 0.001;
const int MIN_DETECTIONS_EXCLUSIVE = 7; 

////////////////////////////////////////////////////////////////////////////////////
// Helper classes
////////////////////////////////////////////////////////////////////////////////////

struct GlobalResources {
  vector<Exposure*> allExposures;
  vector<double> chisqPerDetectionToKeep; // Definition of valid orbit, indexed by # detections
  bool goodFit(Fitter* fitptr) {
    // Criterion we'll use for whether an orbit fit is good enough to pursue:
    double chisqPerDetection = fitptr->getChisq() / (fitptr->nObservations());
    double threshold = (fitptr->nObservations()< chisqPerDetectionToKeep.size()) ?
      chisqPerDetectionToKeep[fitptr->nObservations()] :
      chisqPerDetectionToKeep.back();
    return chisqPerDetection <= threshold;
  }
} globals;

struct Detection {
  Exposure* eptr;
  int objectRow;
  Detection(Exposure* eptr_, int objectRow_): eptr(eptr_), objectRow(objectRow_) {}
  Detection(): eptr(nullptr), objectRow(0) {}
  int id() const {return eptr->id[objectRow];}
  bool operator<(const Detection & rhs) const {return id() < rhs.id();}
};
  
class FitStep {
  // Class that holds an orbit fit to N objects, plus a vector of pointers
  // to exposures that could hold the (N+1)th.
public:
  FitStep(Fitter* fitptr_,
	  list<Exposure*> possibleExposures_,
	  vector<Detection> members_,
	  int nIndependent_,
	  int orbitID_,
	  double fpr_):
    fitptr(fitptr_), possibleExposures(possibleExposures_),
    members(members_), nIndependent(nIndependent_),
    orbitID(orbitID_), fpr(fpr_),
    skippedExposures(false)  {}
  ~FitStep() {delete fitptr;}

  Fitter* fitptr;  // An orbit fit to the N objects

  // Exposures potentially holding the (N+1)th.
  // This list should not include exposures already having a detection, but might
  // include some v close in time.
  list<Exposure*> possibleExposures;
  
  vector<Detection> members;  // Detections that comprise the orbit

  // Flag to set if the search skipped some exposures for
  // exceeding MAXIMUM_FPR_PER_EXPOSURE
  bool skippedExposures;
  
  // Number of detections (so far) with independent asteroids
  // (i.e. don't count consecutive exposures as >1)
  int nIndependent;  

  // A unique ID for each starting orbit
  int orbitID;
  
  // False positive rate for having found last detection,
  // summed over all images searched for it.
  double fpr;

  // Return constituent transient IDs, ascending order
  vector<int> sortedIDs() const {
    vector<int> out;
    out.reserve(members.size());
    for (auto& m : members) 
      out.push_back(m.id());
    sort(out.begin(), out.end());
    return out;
  }
  
  // Return time span of detections
  double arc() const {
    double t0 = members.front().eptr->tdb;
    double t1 = t0;
    for (auto& m : members) {
      double t = m.eptr->tdb;
      if (t<t0) t0=t;
      if (t>t1) t1=t;
    }
    return t1-t0;
  }

  // return all possible N+1 -object FitSteps.
  // They are sorted in order of increasing false positive rate,
  // so the most secure one is first.
  // The total expected number of false positives is returned.
  list<FitStep*> search(double& totalFPR); 
};

struct FitResult {
  // What we'll save for a potentially real orbit
  int orbitID;      // Starting orbit
  int nIndependent; // Different nights in fit
  double fpr;       // False positive rate
  vector<int> detectionIDs;  // IDs of fitted transients, *ascending*
  double chisq;     // Chisq of best fit
  ABG abg;          // orbit
  ABGCovariance invcov; // and covariance
  double arc;       // Time span from first to last detection
  bool skippedExposures; // True if some exposures were skipped for large errors
  bool hasOverlap;  // True if this shares detections with another FitResult
  bool isSecure;    // No doubt that this is real, monopolize detections.
  EIGEN_NEW
  
  FitResult(const FitStep& fs): orbitID(fs.orbitID),
				nIndependent(fs.nIndependent),
				fpr(fs.fpr),
				detectionIDs(fs.sortedIDs()),
				chisq(fs.fitptr->getChisq()),
				abg(fs.fitptr->getABG()),
				invcov(fs.fitptr->getInvCovarABG()),
				arc(fs.arc()),
				skippedExposures(fs.skippedExposures),
				hasOverlap(false),
				isSecure(false)  {}
    
  // Comparison operators are looking at which detections were used to make orbits.
  bool operator==(const FitResult& rhs) const {return detectionIDs==rhs.detectionIDs;}
  bool operator<(const FitResult& rhs) const {return detectionIDs<rhs.detectionIDs;}
  bool includes(const FitResult& rhs) const {return std::includes(detectionIDs.begin(), detectionIDs.end(),
						       rhs.detectionIDs.begin(), rhs.detectionIDs.end());}
  // Do detection lists intersect?
  bool intersects(const FitResult& rhs) const {
    if (detectionIDs.empty() || rhs.detectionIDs.empty()) return false;
    auto ptr1 = detectionIDs.begin();
    auto ptr2 = rhs.detectionIDs.begin();
    while (true) {
      if (*ptr1==*ptr2) {
	return true;
      } else if (*ptr1<*ptr2) {
	if (++ptr1 == detectionIDs.end()) return false;
      } else {
	if (++ptr2 == rhs.detectionIDs.end()) return false;
      }
    }
  }

  void write(std::ostream& os) const {
    stringstuff::StreamSaver ss(os);
    os << "Orbit"
       << (hasOverlap ? "*" : "-")
       << (skippedExposures ? "?" : "-")
       << " " << setw(5) << orbitID
       << " detections " << setw(2) << detectionIDs.size() << " / " << setw(2) << nIndependent
       << " arc " << fixed << setprecision(3) << setw(5) << arc
       << " chisq/DOF " << fixed << setprecision(2) << setw(5) << chisq
       << " / " << setw(2) << 2*detectionIDs.size() - 6
       << " FPR " << defaultfloat << setprecision(3) << fpr
       << endl;
    abg.write(os);
    os << endl;
    // Line giving all detections
    os << "ids:";
    for (auto& m : detectionIDs)
      os << " " << m;
      os << endl;
  }
};


bool
fitStepCompare(const FitStep* lhs, const FitStep* rhs) {
  // A function used for sorting (pointers to)
  // FitStep objects in increasing false positive order.
  return lhs->fpr < rhs->fpr;
}

class FPAccumulator: public DMatrix {
  /** Class to keep track of total false positive expectations.
      This is a 2d array whose rows indicate number of independent
      detections in the orbit, and columns are for different upper
      limits on false positive rate per orbit.
      Elements of the array are the total expected false positive
      counts, i.e. the sum over all contributed orbits with the
      specified number of detections and upper limit on FPR.
  **/
public:
  FPAccumulator(): DMatrix(nMax-nMin+1, nLogFPR, 0.),
		   maxLogFPR(-2.), dLogFPR(-0.5) {}
  static const int nMin=3; // Minimum number of independent detections
  static const int nMax=10; // Max number of independent detections
  const double maxLogFPR; // upper limit on FPR in first column
  const double dLogFPR; // increment to bound per column
  static const int nLogFPR=6; // Number of FPR threshold columns
  void addOrbit(int nIndependent, double fpr) {
    double logFPR = log10(fpr);
    int i = nIndependent-nMin;
    if (i<0) return;
    if (i>=rows()) i=rows()-1; // Last n bin includes all higher n's too
    int jmax = static_cast<int> (floor((log10(fpr)-maxLogFPR)/dLogFPR));
    for (int j=0 ; j<cols() && j<jmax; j++)
      (*this)(i,j) += fpr;
  }
  void write(std::ostream& os) const {
    stringstuff::StreamSaver ss(os);
    os << "# Expected false positives" << endl;
    os << "# N=\\FPR<";
    for (int j=0; j<nLogFPR; j++)
      os << " " << setprecision(2) << setw(7) << std::pow(10., maxLogFPR+j*dLogFPR);
    os << endl;
    os << fixed << setprecision(3);
    for (int i=0; i<this->rows(); i++) {
      os << setw(3) << i+nMin+1 << "      ";
      for (int j=0; j<nLogFPR; j++)
	os << " " << setw(7) << (*this)(i,j);
      os << endl;
    }
  }
};

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

// This method of the FitStep is the key ingredient:
list<FitStep*>
FitStep::search(double& totalFPR) {
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
  
  if (DEBUG) cerr << "Search from nobs " << fitptr->nObservations()
		  << " independent " << nIndependent
		  << " from " << possibleExposures.size() << " exposures"
		  << " last " << members.back().id()
		  << endl;
  // Now walk through all viable exposures, collecting
  // information on the actual matches and the
  // false positive expectation.
  
  totalFPR = 0.; // False positive rate that will be handed to children
  list<Detection> matches;  // List of N+1's.

  for (auto& info : infoList) {
    // Determine the search ellipse area on this exposure
    // with some rough allowance for measurement errors
    double searchArea = PI * sqrt( (info.covXX + NOMINAL_POSITION_ERROR*NOMINAL_POSITION_ERROR)*
				   (info.covYY + NOMINAL_POSITION_ERROR*NOMINAL_POSITION_ERROR) - info.covXY*info.covXY);
				   
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
      skippedExposures = true;
      continue;
    } else {
      // Add this exposure to FPR
      totalFPR += thisFPR;
    }

    if (DEBUG)
      cerr << "Exposure " << info.eptr->expnum
	   << " area " << searchArea * pow(180.*60./3.14159,2)
	   << " fpr " << thisFPR
	   << " minchi " << chisq.minCoeff()
	   << endl;
    if (DEBUG && searchArea!=searchArea)
      cerr << "cov: " << info.covXX << " " << info.covYY << " " << info.covXY
	   << " det " << info.covXX*info.covYY - info.covXY*info.covXY
	   << endl;
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
      if (DEBUG) cerr << "....match exposure " << info.eptr->expnum << " id " << i << " chisq " << chisq[i] << endl;
      auto nextFit = fitptr->augmentObservation(info.eptr->tobs,
						info.eptr->xy(i,0), info.eptr->xy(i,1),
						info.eptr->covXX[i], info.eptr->covYY[i], info.eptr->covXY[i],
						info.eptr->earth,
						true);  // recalculate gravity ???

      // Check chisq to see if we keep it
      if (!globals.goodFit(nextFit)) {
	if (DEBUG) cerr << "....bad fit, deleting" << endl;
	delete nextFit;
	continue;
      }

      // Make a new FitStep that has this detection.
      // Temporarily install the FPR just for this exposure.
      out.push_back(new FitStep(nextFit, possibleExposures, members,
			       (closeInTime ? nIndependent : nIndependent+1),
			       orbitID, thisFPR));
      // Add this detection to the new FitStep
      out.back()->members.push_back(Detection(info.eptr, i));
      // And remove its exposure from the possibleExposures so we don't
      // keep finding ourselves...
      out.back()->possibleExposures.remove(info.eptr);
    } // End of N+1 detection loop
  } // End of N+1 exposure loop

  if (DEBUG) cerr << "->total FPR for search " << totalFPR << endl;

  // Sort the new FitSteps to have lowest FPR first, searched first
  out.sort(&fitStepCompare);
  // Then replace FPR of each with the total for all searches.
  // But if this was not a new independent detection, then
  // we should have it keep the previous FPR.
  for (auto& fsptr : out) {
    if (fsptr->nIndependent > nIndependent)
      fsptr->fpr = totalFPR;
    else
      fsptr->fpr = fpr;
  }
  return out;
}

////////////////////////////////////////////////////////////////////////////////////
// Main begins here
////////////////////////////////////////////////////////////////////////////////////

int
main(int argc, char **argv) {
  string tripletPath;
  string orbitPath;
  string ephemerisPath;
  string transientPath;
  string exposurePath;
  double bindingFactor;
  double searchRadius;
  double maxFPR;  
  bool cullDuplicates;

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
      parameters.addMember("orbitFile",&orbitPath, def,
			   "FITS file for output orbits/links (null->ASCII stdout)", "");
      parameters.addMember("ephemerisFile",&ephemerisPath, def,
			   "SPICE file (null=>environment)", "");
      parameters.addMember("exposureFile",&exposurePath, def,
			   "FITS file holding DECam exposure info", "");
      parameters.addMember("transientFile",&transientPath, def,
			   "FITS file holding transient detections", "");
      parameters.addMember("bindingFactor",&bindingFactor, def | low,
			   "Chisq penalty at unbinding", 4., 0.);
      parameters.addMember("searchRadius",&searchRadius, def | lowopen,
			   "Radius enclosing TNO positions at reference time (degrees)", 4., 0.);
      parameters.addMember("maxFPR",&maxFPR, def | lowopen,
			   "maximum false positive rate to keep", 0.01, 0.);
      parameters.addMember("cull",&cullDuplicates, def,
			   "Use memory to prune duplicate searches?", true);
    }
    parameters.setDefault();

    if (argc>=2 && (string(argv[1])=="-h" || string(argv[1])=="--help")) {
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

    cerr << "# Loading data" << endl;
    // Read the ephemeris
    Ephemeris eph(ephemerisPath);

    // Grab input info
    vector<LONGLONG> orbitIn;
    vector<LONGLONG> objectIn;
    Frame frame;
    double gamma0;
    double dGamma;
    
    bool useStdin = tripletPath.empty();
    {
      double ra0, dec0, tdb0;
      if (useStdin) {
	// Read configuration info from stdin
	string buffer;
	getlineNoComment(cin, buffer);
	{
	  std::istringstream iss(buffer);
	  iss >> ra0 >> dec0 >> tdb0;
	}
	// Next line should be gamma0, dgamma
	getlineNoComment(cin, buffer);
	{
	  std::istringstream iss(buffer);
	  iss >> gamma0 >> dGamma >> tdb0;
	}
      } else {
	FITS::FitsTable ft(tripletPath);
	auto tab = ft.use();
	tab.readCells(orbitIn, "ORBITID");
	tab.readCells(objectIn, "TRANSIENTID");

	// Reference frame defined in header: we will need
	// a central RA, Dec, and reference epoch.
	tab.header()->getValue("RA0",ra0);
	tab.header()->getValue("DEC0",dec0);
	tab.header()->getValue("TDB0",tdb0);

	// Read the central gamma and half-width from header too
	tab.header()->getValue("GAMMA0",gamma0);
	tab.header()->getValue("DGAMMA0",dGamma);
      }
      astrometry::Orientation orient(astrometry::SphericalICRS(ra0*DEGREE, dec0*DEGREE));
      orient.alignToEcliptic();
      frame = Frame(eph.observatory(807,tdb0), orient, tdb0);
    }      
    if (gamma0<=0. || dGamma<=0. || dGamma > gamma0) {
      cerr << "Invalid gamma0, dGamma" << gamma0 << ", " << dGamma << endl;
      exit(1);
    }
    //??transientPath = "/Users/garyb/DES/TNO/zone029.transients.fits";
    
    {
      // Choose chisq cutoffs for numbers of detections
      DVector pctile(8,0.);
      // Start with the values for 99th pctile of
      // chisq distribution for 2*n-5 DOF,
      // assuming that gammadot isn't really helping fit.
      pctile << 0., 0., 0., 6.6, 11.3, 15.1, 18.5, 21.7;
      // Multiply by a factor for error underestimates
      pctile *= 1.4;
      // Add a few to allow for gamma / binding priors
      pctile.array() += bindingFactor + 4.;
      // Save, dividing by number of points
      globals.chisqPerDetectionToKeep.resize(pctile.size());
      for (int i=0; i<3; i++)
	globals.chisqPerDetectionToKeep[i] = 0.; // Need >=3 pts
      for (int i=3; i<pctile.size(); i++)
	globals.chisqPerDetectionToKeep[i] = pctile[i]/i;
    }

    // We're going to keep track of the total expected false positive count,
    // as a function of the number of independent epochs.
    FPAccumulator accumulator;
    
    // Pluck out all relevant exposure data.  Then close
    // transient catalog and exposure catalog.
    vector<Exposure*> exposurePool;
    {
      ExposureTable et(exposurePath);
      exposurePool = et.getPool(frame,
				eph,
				gamma0,
				dGamma,
				searchRadius,
				true); // astrometric exposures only
    }
    // Get transient catalogs for all these exposures, deleting pool members with none
    // And build a map from object ID back to the exposure it
    // was taken in. 
    std::map<int, Detection> detectionIndex;
    list<Exposure*> allExposureList;  // List of valid exposures
    {
      TransientTable tt(transientPath);
      for (auto iter = exposurePool.begin();
	   iter != exposurePool.end();
	   ++iter) {
	if (tt.fillExposure(frame, *iter)) {
	  // Add all these transients to the index
	  for (int i=0; i<(*iter)->id.size(); i++)
	    detectionIndex[(*iter)->id[i]] = Detection(*iter,i);
	  // and add exposure to list of useful ones
	  allExposureList.push_back(*iter);
	} else {
	  delete *iter;
	  *iter = nullptr;
	}
      }
    }
	  
    int orbitID = -1;  // Starting point for an orbitID counter, if using stdin

    // These are orbitID's that have already been "cleaned out" from an
    // excellent orbit and should be henceforth ignored.
    set<int> settledOrbitIDs;
    
    // Now enter a loop of growing fits.
    list<FitStep*> fitQueue;

    // And build a list of successful linkages
    list<FitResult, Eigen::aligned_allocator<FitResult>> savedResults;
    
    // Here we'll keep a list of detection combinations that have already been fit
    // and have no need of repeating.  This saves time but takes memory.
    // Set cull=false at command line to *not* do this.
    set<vector<int>> alreadySearched;
      
    int nextInputRow = 0;  // Counter used if we are getting input from FITS file

    cerr << "# Starting linking" << endl;
    while (true) {
      if (fitQueue.empty()) {
	// Stock the queue with fresh orbit from input
	// Continue reading new orbits and trying them out
	// until we are out of input or find one with decent initial fit
	while (true) {
	  // Gather object ID's for next orbit trial
	  vector<int> memberIDs;
	  if (useStdin) {
	    string buffer;
	    /* Loop over all initial detection sets: ?? simple stdin format to start */
	    if (!getlineNoComment(cin, buffer)) {
	      break;  // No more input.
	    }
	    istringstream iss(buffer);
	    // Read object id's on a line
	    int objectID;
	    while (iss >> objectID) {
	      memberIDs.push_back(objectID);
	    }
	    orbitID++;
	  } else {
	    // Read members from table
	    if (nextInputRow >= orbitIn.size())
	      break;  // No more input available.
	    orbitID = orbitIn[nextInputRow];
	    for ( ;nextInputRow < orbitIn.size() && orbitIn[nextInputRow]==orbitID ; nextInputRow++)
	      memberIDs.push_back(objectIn[nextInputRow]);
	  }
	  
	  // Test object ID's for validity
	  vector<Detection> members;
	  for (auto id : memberIDs) {
	    if (!detectionIndex.count(id)) {
	      // No detection with this ID - probably something
	      // we're rejecting from large errors.  Just skip
	      continue;
	    }
	    members.push_back(detectionIndex[id]);
	  }

	  // Test the initial detections for sane orbit
	  if (members.size() < 3)
	    continue; // don't try these...

	  if (DEBUG) cerr << "Making initial with " << members.size() << " detections" << endl;
	  // Make an orbit fit
	  auto fitptr = new Fitter(eph, Gravity::BARY);  // ?? Need giants' gravity??
	  fitptr->setFrame(frame);
	  if (bindingFactor > 0.)
	    fitptr->setBindingConstraint(bindingFactor);
	  fitptr->setGammaConstraint(gamma0, dGamma/2.); // Give gamma prior so search bound is 2 sigma
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
	      int j = det.objectRow;
	      tobs[i] = det.eptr->tobs;
	      thetaX[i] = det.eptr->xy(j,0);
	      thetaY[i] = det.eptr->xy(j,1);
	      covxx[i] = det.eptr->covXX[j];
	      covyy[i] = det.eptr->covYY[j];
	      covxy[i] = det.eptr->covXY[j];
	      xE.row(i) = det.eptr->earth.transpose();
	    }
	    fitptr->setObservationsInFrame(tobs, thetaX, thetaY, covxx, covyy, covxy, xE);
	  }
	  bool goodFit = true;
	  try {
	    fitptr->setLinearOrbit();
	    fitptr->setLinearOrbit();
	    if (DEBUG) cerr << "Initial linear fit chisq " << fitptr->getChisq()
		     << " " << fitptr->getABG(true)[ABG::G] << endl;
	    fitptr->newtonFit();  // ?? any more elaborate fitting needed ??
	    if (DEBUG) cerr << "Newton fit chisq " << fitptr->getChisq()
		     << " " << fitptr->getABG()[ABG::G] << endl;
	  } catch (std::runtime_error& e) {
	    goodFit = false;
	  }

	  // Skip the orbit if fit failed or is poor
	  if (!goodFit || !globals.goodFit(fitptr)) {
	    if (DEBUG) cerr << "Not a good fit " << goodFit << endl;
	    if (DEBUG) fitptr->printResiduals(cerr);
	    delete fitptr;
	    continue;
	  }

	  // Fit is good - queue up a FitStep
	  // First make a list of possible exposures that excludes the
	  // ones that the original detections come from.
	  list<Exposure*> possibleExposures(allExposureList);
	  for (auto& m : members) {
	    possibleExposures.remove(m.eptr);
	  }
	  // Now queue up
	  if (DEBUG) cerr << "Queuing using " << possibleExposures.size() << " exposures" << endl;
	  fitQueue.push_back(new FitStep(fitptr, possibleExposures, members,
					 members.size(), orbitID,
					 1000.)); // Set a high FPR so we don't trigger output on this
	} // End of input-reading loop.

      } // End of searching for new inputs for queue.

      if (fitQueue.empty())
	break; 	  // Nothing left to fit!

      // Search the next item on the queue:
      if (DEBUG) cerr << "-->New search, queue size " << fitQueue.size() << endl;
      auto thisFit = fitQueue.front();
      fitQueue.pop_front();
      if (settledOrbitIDs.count(thisFit->orbitID)) {
	// We do not need to pursue this fit, it's finished elsewhere
	delete thisFit;
	continue;
      }
      if (cullDuplicates) {
	vector<int> ids = thisFit->sortedIDs();
	if (alreadySearched.count(ids)) {
	  // Duplicate.  Skip.
	  delete thisFit;
	  continue;
	} else {
	  alreadySearched.insert(ids);
	}
      }
      double totalFPR;
      auto newFits = thisFit->search(totalFPR);

      // Add the FPR of this search into the total for
      // searches with the same number of independent detections
      accumulator.addOrbit(thisFit->nIndependent, totalFPR);

      if (DEBUG) cerr << "search returns " << newFits.size() << endl;
      if (newFits.empty()) {
	// No new points can be added to this fit.  So it's a completed search.

	// Is this a secure orbit, likely to have all relevant detections?
	bool isSecure = thisFit->nIndependent >= MIN_DETECTIONS_EXCLUSIVE
	  && thisFit->fpr <= MAX_FPR_EXCLUSIVE;
	// Is it worthy of output?
	if (thisFit->nIndependent >= MIN_DETECTIONS_TO_OUTPUT
	    && thisFit->fpr <= maxFPR) {
	  savedResults.push_back(FitResult(*thisFit));
	  if (isSecure) savedResults.back().isSecure = true;
	  if (DEBUG) savedResults.back().write(cout);
	}
	// If this is a secure fit, shut down further use of its detections
	if (isSecure) {
	  for (auto& m : thisFit->members) {
	    m.eptr->valid[m.objectRow] = false;
	  }
	  settledOrbitIDs.insert(thisFit->orbitID);
	}
      } else {
	// Not a terminal fit - add new searches to queue.
	// Put them at the front so we are searching depth-first
	fitQueue.insert(fitQueue.begin(), newFits.begin(), newFits.end());
      }
      // Done with this fit
      delete thisFit;
    } // End of the queue loop.

    // Done searching!  Clean up the results by
    // deleting anything whose detections are a subset of another,
    // and mark those that share detections with another.

    cerr << "# Starting the purge" << endl;
    for (auto ptr1 = savedResults.begin();
	 ptr1 != savedResults.end(); ) {
      bool dup = false; // Set if this is contained in another
      auto ptr2 = ptr1;
      ++ptr2;
      while (ptr2!=savedResults.end()) {
	if (ptr2->includes(*ptr1)) {
	  // ptr1 is a subset, delete it
	  ptr1 = savedResults.erase(ptr1);
	  dup = true;
	  break;
	} else if (ptr1->includes(*ptr2)) {
	  // ptr2 is a subset, delete it
	  ptr2 = savedResults.erase(ptr2);
	} else if (ptr1->intersects(*ptr2)) {
	  if (ptr1->isSecure && ptr2->isSecure) {
	    // Weird to have distinct secure orbits overlap.  Keep both
	    cout << "# WEIRD - two secure orbits overlap" << endl;
	    ptr1->hasOverlap = true;
	    ptr2->hasOverlap = true;
	    ++ptr2;
	  } else if (ptr1->isSecure) {
	    // Eliminate #2
	    ptr2 = savedResults.erase(ptr2);
	  } else if (ptr2->isSecure) {
	    // Eliminate #1
	    ptr1 = savedResults.erase(ptr1);
	    dup = true;
	    break;
	  } else {
	    // Neither is secure, just mark both
	    ptr1->hasOverlap = true;
	    ptr2->hasOverlap = true;
	    ++ptr2;
	  }
	} else {
	  // Nothing interesting
	  ++ptr2;
	}
      }
      if (!dup) ++ptr1;
    }
    
    // Print all results
    for (auto& r : savedResults)
      r.write(cout);
	    
    // Report total false positive rates, then quit
    accumulator.write(cout);
    exit(0);

  } catch (std::runtime_error& e) {
    quit(e,1);
  }
}  
