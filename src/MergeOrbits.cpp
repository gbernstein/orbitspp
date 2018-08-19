// Program to do cleanup of orbit sets


#include <iostream>

#include "StringStuff.h"
#include "Astrometry.h"
#include "Fitter.h"
#include "Elements.h"
#include "Exposures.h"
#include "Pset.h"
#include "FitsTable.h"
#include "FitsImage.h"
#include "Eigen/StdVector"

using namespace std;
using namespace orbits;

const int DEBUGLEVEL=2;

const string usage =
  "MergeOrbits:\n"
  "\n"
  "Usage: MergeOrbits <orbitfile1> [orbitfile2] ... [-parameter value]\n"
  "   <orbitfileN> are any number >=1 of outputs from GrowOrbits\n"
  "   [-parameter value] are additional parameter name/value pairs.\n"
  "   Parameters are listed below.\n"
  "\n"
  "Outputs: New FITS file with orbit and detection tables that are\n"
  "        merged and updated, including orbital elements,\n"
  "        a custom reference frame for each orbit,\n"
  "        and information on predicted positions on exposures";
// ....

struct FitResult {
  // Everything we need to know about an orbit
  string inputFile; // Orbit file it came from
  LONGLONG inputID;      // Starting orbit ID
  int inputIndependent; // Different nights in original fit
  double fpr;       // False positive rate in original search

  vector<int> detectionIDs;  // IDs of fitted detections, *ascending*
  // The observations in ICRS coords:
  vector<Observation, Eigen::aligned_allocator<Observation>> obsICRS;  
  double arc;       // Time span from first to last detection

  class Opportunity {
  public:
    Opportunity(): eptr(nullptr), orbitCov(0.), hasDetection(false) {}
    // Class representing an exposure that could have seen this
    Exposure* eptr;
    astrometry::Gnomonic orbitPred;  // Orbit prediction (with frame Orientation)
    Matrix22 orbitCov;   // Model covariance in frame
    vector<int> ccdnums;
    bool hasDetection;   // set if there is a detection on the exposure
  };
  vector<Opportunity, Eigen::aligned_allocator<Opportunity>> opportunities;

  double chisq;     // Chisq of best fit
  Frame frame;       // Frame for its ABG
  ABG abg;          // orbit
  ABGCovariance invcov; // and covariance
  Elements el;
  ElementCovariance elCov;

  int friendGroup;  // Number of its overlap group (-1=loner)
  bool skippedExposures; // True if some exposures were skipped for large errors
  bool hasOverlap;  // True if this shares detections with another FitResult
  bool isSecure;    // No doubt that this is real, monopolize detections.
  bool changedDetectionList; // set if detections were added / dropped
  EIGEN_NEW

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
  // This is the main function of the program:
  bool fitAndFind(const Ephemeris& ephem,
		  double tdb0,
		  TransientTable& transients,
		  ExposureTable& exposures,
		  vector<Exposure*>& pool);

};

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
		   maxLogFPR(-1.5), dLogFPR(-0.5) {}
  static const int nMin=3; // Minimum number of independent detections
  static const int nMax=10; // Max number of independent detections
  const double maxLogFPR; // upper limit on FPR in first column
  const double dLogFPR; // increment to bound per column
  static const int nLogFPR=7; // Number of FPR threshold columns

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
  img::Image<float> getImage() const {
    img::Image<float> out(this->rows(), this->cols());
    out = 0.;
    for (int i=out.xMin(); i<=out.xMax(); i++)
      for (int j=out.yMin(); j<=out.yMax(); j++) {
	out(i,j) = (*this)(i-1,j-1);  // FITS is 1-indexed
      }
    int tmp = nMin;
    out.setHdrValue("NMIN",tmp);
    out.setHdrValue("MAXLOGFP",maxLogFPR);
    out.setHdrValue("DLOGFP",dLogFPR);
    return out;
  }
};

bool
FitResult::fitAndFind(const Ephemeris& ephem,
		      double tdb0,
		      TransientTable& transients,
		      ExposureTable& exposureTable,
		      vector<Exposure*>& pool) {
  // First fit an orbit to the detectionID's currently in the vector.
  // Then search all exposures for any additional transients.
  // Return true if any detections changed.

  const double SEARCH_CHISQ=9.; // Maximum chisq to consider a match in 2d
  // Acquire Observation for each detected transient
  obsICRS.clear();
  if (DEBUGLEVEL>1)
    cerr << "# fitAndFind: Reading detection info" << endl;
  for (auto id : detectionIDs) {
    obsICRS.push_back(transients.getObservation(id, ephem, exposureTable));
  }

  if (DEBUGLEVEL>1)
    cerr << "# ...Calculating arc length" << endl;
  // Calculate arc
  {
    double t1 = obsICRS.front().tdb;
    double t2 = t1;
    for (auto& obs : obsICRS) {
      if (obs.tdb < t1) t1 = obs.tdb;
      if (obs.tdb > t2) t2 = obs.tdb;
    }
    arc = t2 - t1;
  }

  if (DEBUGLEVEL>1)
    cerr << "# ...reference frame" << endl;

  // Choose a reference frame for this orbit.  Use provided TDB,
  // and mean RA/Dec.  Get this mean by summing 3d coords
  {
    Vector3 sum(0.);
    for (auto& obs : obsICRS)
      sum += obs.radec.getUnitVector();
    astrometry::CartesianICRS ss(sum);
    astrometry::SphericalICRS pole(ss);
    frame.orient = astrometry::Orientation(pole);
    frame.orient.alignToEcliptic();
    frame.tdb0 = tdb0;
    frame.origin = ephem.observatory(807, frame.tdb0);
  }

  if (DEBUGLEVEL>1)
    cerr << "# ...fitting in frame at "
	 << frame.orient.getPole() << endl;

  // Execute fit
  Fitter fit(ephem, Gravity::GIANTS);
  for (auto& obs : obsICRS)
    fit.addObservation(obs);
  fit.setFrame(frame);
  // ! No priors on the orbit.
  fit.setLinearOrbit();
  fit.newtonFit();

  if (DEBUGLEVEL>1)
    cerr << "# ...Fit complete" << endl;

  // Save fit results
  chisq = fit.getChisq();
  abg = fit.getABG();
  invcov = fit.getInvCovarABG();
  el = fit.getElements();
  elCov = fit.getElementCovariance();

  // Now go fishing for exposures this orbit crosses
  int nExposures = pool.size();
  DVector tdbAll(nExposures);
  DMatrix earthAll(nExposures,3);  // ICRS observatory position
  DMatrix axisAll(nExposures,3);   // ICRS direction cosines of pointings
  for (int i=0; i<nExposures; i++) {
    tdbAll[i] = pool[i]->tdb;
    earthAll.row(i) = pool[i]->earthICRS.getVector().transpose();
    axisAll.row(i) = pool[i]->localICRS.getPole().getUnitVector().transpose();
  }

  // convert field radius to a maximum chord length between
  // unit vectors
  const double fieldRadius = 1.15*DEGREE; // A bit larger to be complete
  const double maxChordSq = pow(2. * sin(fieldRadius/2.),2.);

  DMatrix target = fit.getTrajectory().observe(tdbAll, earthAll);
  // Find those within radius
  BVector hits = (target - axisAll).rowwise().squaredNorm().array() < maxChordSq;
  int nHits = hits.array().count();
  
  if (DEBUGLEVEL>1)
    cerr << "# ...exposure search has " << nHits << " hits" << endl;

  // Collect information on all candidate exposures (??? repetitive reading...)
  opportunities.clear();
  tdbAll.resize(nHits);
  earthAll.resize(nHits,3);
  for (int i=0; i<hits.size(); i++) {
    if (hits[i]) {
      Opportunity opp;
      opp.eptr = pool[i];
      tdbAll[opportunities.size()] = pool[i]->tdb;
      earthAll.row(opportunities.size()) = pool[i]->earthICRS.getVector().transpose();
      opp.hasDetection = false;
      opportunities.push_back(opp);
    }
  }
  // ??? Did we get all the exposures from original orbit?
  
  if (DEBUGLEVEL>1)
    cerr << "# ...getting predictions" << endl;

  // Predict orbit position at each exposure
  // Convert into current frame:
  tdbAll.array() -= frame.tdb0;
  earthAll = frame.fromICRS(earthAll);
  DVector xPred(nHits);
  DVector yPred(nHits);
  DVector covxxPred(nHits);
  DVector covyyPred(nHits);
  DVector covxyPred(nHits);
  fit.predict(tdbAll, earthAll, &xPred, &yPred,
	      &covxxPred, &covyyPred, &covxyPred);

  // Save info on every exposure, find new detections
  vector<int> newDetectionIDs;
  

  if (DEBUGLEVEL>1)
    cerr << "# ...finding detections" << endl;

  for (int i=0; i<nHits; i++) {
    Opportunity& opp = opportunities[i];
    opp.orbitPred = astrometry::Gnomonic(xPred[i], yPred[i], frame.orient);
    opp.orbitCov(0,0) = covxxPred[i];
    opp.orbitCov(0,1) = covxyPred[i];
    opp.orbitCov(1,0) = covxyPred[i];
    opp.orbitCov(1,1) = covyyPred[i];

    /**/cerr << "call whichCCDs on " << opp.eptr->expnum << endl;
    // Check in detail for error ellipse crossing a CCD
    opp.ccdnums = opp.eptr->whichCCDs(opp.orbitPred, opp.orbitCov*SEARCH_CHISQ);

    /**/cerr << "filling transients" << endl;
    // Load transient and corner information into each exposure, using frame
    transients.fillExposure(frame, opp.eptr);

    /**/cerr << "calculate chisq" << endl;
    // Find detections - did any change?? - do CCDNUM match prev?
    DVector allChi = opp.eptr->chisq(xPred[i], yPred[i],
				     covxxPred[i], covyyPred[i], covxyPred[i]);
    /**/if (allChi.size()>0) cerr << "min chisq " << allChi.array().minCoeff() << endl;
    for (int iTrans=0; iTrans<allChi.size(); iTrans++) {
      if (allChi[iTrans]<SEARCH_CHISQ) {
	/**/cerr << "Found detection!" << endl;
	opp.hasDetection = true;
	newDetectionIDs.push_back(opp.eptr->id[iTrans]);
	// Save residuals for each detection ?? (have the info to do this later)
      }
    }
    /**/cerr << "...done " << i << endl;
  }

    /**/cerr << "sorting" << endl;
  // Compare new and old detection lists
  std::sort(newDetectionIDs.begin(), newDetectionIDs.end());
  std::sort(detectionIDs.begin(), detectionIDs.end());

  bool changed = (newDetectionIDs==detectionIDs);

  if (DEBUGLEVEL>1)
    cerr << "# ...changed detections? " << changed << endl;

  if (changed) changedDetectionList = true;

  // ?? Check for detections lost from not on an candidate exposure??

  // Replace detection list with new one.
  detectionIDs = newDetectionIDs;
  return changed;
}
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

int main(int argc,
	 char *argv[])
{
  // Controlling parameters / constants
  string ephemerisPath;
  string exposurePath;
  string cornerPath;  // File with CCD corners per exposure

  string transientPath;

  double ra0;    // Center of field (ignore exposures >60 degrees away)
  double dec0;
  double radius;
  double tdb0;   // Reference time (=epoch of orbits or state vectors)
  // Else orbital elements

  Pset parameters;
   
  try {
    {
      const int def=PsetMember::hasDefault;
      const int low=PsetMember::hasLowerBound;
      const int up=PsetMember::hasUpperBound;
      const int lowopen = low | PsetMember::openLowerBound;
      const int upopen = up | PsetMember::openUpperBound;

      parameters.addMember("transientFile",&transientPath, def,
			   "transient catalog (null=>environment", "");
      parameters.addMember("ephemerisFile",&ephemerisPath, def,
			   "SPICE file (null=>environment", "");
      parameters.addMember("exposureFile",&exposurePath, def,
			   "DECam exposure info file (null=>environment)", "");
      parameters.addMember("cornerFile",&cornerPath, def,
			   "CCD corners table file (null=>environment)", "");
      parameters.addMember("ra0",&ra0, def,
			   "RA of field center (deg)", 10.);
      parameters.addMember("dec0",&dec0, def,
			   "Dec of field center (deg)", -20.);
      parameters.addMember("radius",&radius, def || lowopen || up,
			   "Radius of search area (deg)", 85.,0.,85.);
      parameters.addMember("TDB0",&tdb0, def,
			   "TDB of reference time (=orbit epoch), yrs since J2000", 0.);
    }
    parameters.setDefault();

    if (argc<2 || (string(argv[1])=="-h" || string(argv[1])=="--help")) {
      cout << usage << endl;
      parameters.dump(cerr);
      exit(1);
    }

    vector<string> orbitFiles;
    // Get names of orbitFiles
    for (int iarg=1; iarg < argc && argv[iarg][0]!='-'; iarg++)
      orbitFiles.push_back(argv[iarg]);
    
    // And read arguments from the remaining command line entries
    parameters.setFromArguments(argc, argv);

    // Read the ephemeris
    Ephemeris ephem(ephemerisPath);

    // Establish a reference frame
    Frame frame;
    {
      astrometry::SphericalICRS pole(ra0*DEGREE, dec0*DEGREE);
      frame.orient = astrometry::Orientation(pole);
      frame.orient.alignToEcliptic();
      frame.tdb0 = tdb0;
      frame.origin = ephem.observatory(807, tdb0);
    }

    // Read the exposure table
    ExposureTable et(exposurePath);

    // Get all exposures of relevance and build big arrays
    vector<Exposure*> pool = et.getPool(frame, ephem, radius*DEGREE);

    cerr << "# Read " << pool.size() << " exposures" << endl;
    cerr << "# reading corners" << endl;
    
    // Load CCD corners
    {
      CornerTable tt(cornerPath);
      for (auto& eptr : pool)
	tt.fillExposure(eptr);
    }

    cerr << "# reading transients" << endl;

    // Open the transients table
    TransientTable transients(transientPath);
    
    // Start reading input orbits, giving each one its own
    // reference frame
    list<FitResult*> orbits;

    for (auto orbitPath : orbitFiles) {
      cerr << "# Reading orbits from " << orbitPath << endl;
      // Open the tables/files
      img::FTable inputOrbitTable;
      img::FTable inputObjectTable;
      {
        FITS::FitsTable ft(orbitPath, FITS::ReadOnly,"ORBITS");
	inputOrbitTable = ft.extract();
      }
      {
        FITS::FitsTable ft(orbitPath, FITS::ReadOnly, "OBJECTS");
	inputObjectTable = ft.extract();
      }

      // Sum the Accumulators ???
      
      int objectTableRow = 0; // Current location in object table

      // Step through orbits
      for (int row=0; row<inputOrbitTable.nrows(); row++) {
	auto orb = new FitResult;
	orb->inputFile = orbitPath;
	orb->inputID = row;
	inputOrbitTable.readCell(orb->inputIndependent, "NUNIQUE", row);
	inputOrbitTable.readCell(orb->fpr, "FPR", row);
	
	// Collect known transients for this orbit.
	while (true) {
	  LONGLONG id;
	  inputObjectTable.readCell(id, "ORBITID", objectTableRow);
	  if (id!=orb->inputID) break;
	  int objectID;
	  inputObjectTable.readCell(objectID, "OBJECTID", objectTableRow);
	  orb->detectionIDs.push_back(objectID);
	  ++objectTableRow;
	}

	// Set some initial flags
	orb->friendGroup = -1;
	orb->skippedExposures = false;
	orb->hasOverlap = false;
	orb->isSecure = false;
	orb->changedDetectionList = false;

	// Iterate fit and search
	int MAX_FIT_ITERATIONS=5;
	for (int iter=0; iter<MAX_FIT_ITERATIONS; iter++) {
	  if (!orb->fitAndFind(ephem, tdb0, transients, et, pool))
	    break;
	}

	orbits.push_back(orb);
      }
    }

    // Purge duplicates and subsets, group overlapping sets by friends
    // of friends.
    list<list<FitResult*>> friendGroups;
    
    cerr << "# Starting the purge" << endl;
    for (auto ptr1 : orbits) {
      list<FitResult*> itsFriends; // Collect all friends of this one.
      itsFriends.push_back(ptr1);

      // Loop through previous orbits, finding friends/supersets/subsets
      for (auto groupPtr = friendGroups.begin();
	   groupPtr != friendGroups.end(); ) {
	bool mergeGroup = false; // Set if overlap with any friend
	bool isSubset = false; // Set if redundant orbit
	// Start loop over members of friend group
	for (auto pptr2 = groupPtr->begin();
	     pptr2 != groupPtr->end(); ) {
	  auto ptr2 = *pptr2; // pptr2 is iterator to pointer...
	  if (ptr1->intersects(*ptr2)) {
	    // Merge friend groups
	    mergeGroup = true;
	    // Check each one as subset
	    if (ptr2->includes(*ptr1)) {
	      // ptr1 is a subset, delete it, no need to continue
	      itsFriends.remove(ptr1);
	      isSubset = true;
	      break;
	    } else if (ptr1->includes(*ptr2)) {
	      // ptr2 is a subset, delete it from its group
	      pptr2 = groupPtr->erase(pptr2);
	      // And delete the whole thing
	      delete ptr2;
	    } else {
	      // Proceed to next 2nd orbit
	      ++pptr2;
	    }
	  } else {
	    // No intersection, proceed to next 2nd orbit
	    ++pptr2;
	  }
	} // End loop over members of friend group
	if (mergeGroup) {
	  // Absorb 2nd group and delete it
	  itsFriends.insert(itsFriends.end(),
			    groupPtr->begin(), groupPtr->end());
	  groupPtr = friendGroups.erase(groupPtr);
	} else {
	  // Advance to next group
	  ++groupPtr;
	}
	if (isSubset) {
	  // ptr1 is a duplicate, no need to proceed through other groups
	  delete ptr1;
	  break;
	}
      } // End of loop over friend groups.
      // Add this group to friend list
      friendGroups.push_back(itsFriends);
    } // End loop over all orbits.

    // Now all orbits are grouped into friend sets.
    // Do any purging of overlapping orbits ??
    
    // Create result arrays/tables ??
    
    // Loop through surviving orbits
    int groupCounter = -1;
    for (auto& group : friendGroups) {
      ++groupCounter;
      for (auto& orb : group) {
	// Process each output orbit:
	orb->friendGroup = groupCounter;
	// save detections, get residuals, chisq per point, and band, mag, magerr
	// ???
	cout << orb->inputID
	     << " " << orb->detectionIDs.size()
	     << " " << fixed << setprecision(2) << setw(5) << orb->chisq
	     << " " << orb->fpr
	     << " " << orb->friendGroup
	     << " " << setprecision(6) << setw(8) << orb->abg[ABG::G];
	for (auto id : orb->detectionIDs)
	  cout << " " << id;
	cout << endl;

	delete orb;
      }
    }
    
  } catch (std::runtime_error& e) {
    quit(e);
  }
  exit(0);
}

    

