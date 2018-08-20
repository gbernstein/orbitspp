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

const int DEBUGLEVEL=1;
const int MIN_UNIQUE=4; // min number of independent detections to retain orbit

// Time that must pass between exposures to be considered independent detections
// (e.g. when asteroids or defects would have moved out of linking range)
const double INDEPENDENT_TIME_INTERVAL = 0.1*DAY;

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
  int inputUnique; // Different nights in original fit
  double fpr;       // False positive rate in original search

  vector<int> detectionIDs;  // IDs of fitted detections, *ascending*
  // The observations in ICRS coords:
  int nUnique;      // Number of distinct detection times
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
  // ?? maybe a little
  fit.setBindingConstraint(1.);
  fit.setLinearOrbit();
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
  nUnique = 0; // reset counter of unique times
  vector<double> timesOfDetections; // used for counting unique

  if (DEBUGLEVEL>1)
    cerr << "# ...finding detections" << endl;

  for (int i=0; i<nHits; i++) {
    Opportunity& opp = opportunities[i];
    opp.orbitPred = astrometry::Gnomonic(xPred[i], yPred[i], frame.orient);
    opp.orbitCov(0,0) = covxxPred[i];
    opp.orbitCov(0,1) = covxyPred[i];
    opp.orbitCov(1,0) = covxyPred[i];
    opp.orbitCov(1,1) = covyyPred[i];

    // Check in detail for error ellipse crossing a CCD
    opp.ccdnums = opp.eptr->whichCCDs(opp.orbitPred, opp.orbitCov*SEARCH_CHISQ);

    // Load transient and corner information into each exposure, using frame
    transients.fillExposure(frame, opp.eptr);

    DVector allChi = opp.eptr->chisq(xPred[i], yPred[i],
				     covxxPred[i], covyyPred[i], covxyPred[i]);
    for (int iTrans=0; iTrans<allChi.size(); iTrans++) {
      if (allChi[iTrans]<SEARCH_CHISQ && opp.eptr->valid[iTrans]) {
	opp.hasDetection = true;
	newDetectionIDs.push_back(opp.eptr->id[iTrans]);
	// Is this a unique detection?
	bool isUnique = true;
	for (auto tdb : timesOfDetections) {
	  if ( abs(tdb-tdbAll[i]) < INDEPENDENT_TIME_INTERVAL) {
	    isUnique = false;
	    break;
	  }
	}
	if (isUnique)
	  nUnique++;
	else
	  timesOfDetections.push_back(tdbAll[i]);
      }
    }
  }

  // Compare new and old detection lists
  std::sort(newDetectionIDs.begin(), newDetectionIDs.end());
  std::sort(detectionIDs.begin(), detectionIDs.end());

  bool changed = (newDetectionIDs!=detectionIDs);

  if (DEBUGLEVEL>1) {
    set<int> s1;
    s1.insert(detectionIDs.begin(), detectionIDs.end());
    set<int> s2;
    s2.insert(newDetectionIDs.begin(), newDetectionIDs.end());
    cerr << "# ...detections lost: ";
    set<int> result;
    set_difference(s1.begin(), s1.end(), s2.begin(), s2.end(),
		   std::inserter(result, result.end()));
    for (auto id: result) cerr << " " << id;
    cerr << endl;
    cerr << "# ...detections gained: ";
    result.clear();
    set_difference(s2.begin(), s2.end(), s1.begin(), s1.end(),
		   std::inserter(result, result.end()));
    for (auto id: result) cerr << " " << id;
    cerr << endl;
  }
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
  string outPath; // location of output file

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

      parameters.addMember("o",&outPath, 0,
			   "Output orbit/object table");
      parameters.addMember("transientFile",&transientPath, def,
			   "transient catalog (null=>environment)", "");
      parameters.addMember("ephemerisFile",&ephemerisPath, def,
			   "SPICE file (null=>environment)", "");
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
			   "TDB of reference time (=orbit epoch), yrs since J2000", 16.);
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

    img::Image<float> fpImage; // Summed false-positive image
    
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
      {
	// Sum the false-positive estimates
	img::FitsImage<float> fi(orbitPath, FITS::ReadOnly, "FP");
	if (!fpImage.getBounds()) {
	  // first image
	  fpImage = fi.extract();
	} else {
	  fpImage += fi.extract();
	}
      }
	
      int objectTableRow = 0; // Current location in object table

      // Step through orbits
      for (int row=0; row<inputOrbitTable.nrows(); row++) {
	auto orb = new FitResult;
	orb->inputFile = orbitPath;
	orb->inputID = row;
	inputOrbitTable.readCell(orb->inputUnique, "NUNIQUE", row);
	orb->nUnique = orb->inputUnique;
	inputOrbitTable.readCell(orb->fpr, "FPR", row);
	
	// Collect known transients for this orbit.
	while (true) {
	  LONGLONG id;
	  inputObjectTable.readCell(id, "ORBITID", objectTableRow);
	  if (id!=orb->inputID) break;
	  int objectID;
	  inputObjectTable.readCell(objectID, "OBJECTID", objectTableRow);
	  if (objectID<0) continue;  // negative ID is non-detection
	  orb->detectionIDs.push_back(objectID);
	  ++objectTableRow;
	}

	// Set some initial flags
	orb->friendGroup = -1;
	orb->changedDetectionList = false;

	bool success=false; // We will throw away fitting failures
	// Iterate fit and search
	int MAX_FIT_ITERATIONS=5;
	try {
	  for (int iter=0; iter<MAX_FIT_ITERATIONS; iter++) {
	    if (!orb->fitAndFind(ephem, tdb0, transients, et, pool)) {
	      success = true;
	      break;
	    }
	    // Also quit if we no longer have enough detections
	    if (orb->nUnique < MIN_UNIQUE) 
	      break;
	    if (iter==MAX_FIT_ITERATIONS-1) cerr << "# WARNING: Not converging at " << orb->inputID << endl;
	  }
	} catch (std::runtime_error& e) {
	  cerr << "# Fitting failure " << orb->inputFile << "/" << orb->inputID
	       << " nUnique " << orb->nUnique << " arc " << orb->arc << endl;
	  cerr << e.what() << endl;
	}

	if (success)
	  orbits.push_back(orb);
	else
	  delete orb;
      }
    }

    // Purge duplicates and subsets, group overlapping sets by friends
    // of friends.
    list<list<FitResult*>> friendGroups;
    
    cerr << "# Starting the purge" << endl;
    for (auto ptr1 : orbits) {
      if (ptr1->nUnique < MIN_UNIQUE) { 
	// Throw it away
	delete ptr1;
	continue;
      }

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
	    // Check relations between pairs
	    if ((*ptr1)==(*ptr2)) {
	      // If both are same orbit, keep lowest FPR
	      if (ptr1->fpr > ptr2->fpr) {
		// Delete ptr1
		itsFriends.remove(ptr1);
		isSubset = true;
		break;
	      } else {
		// Delete ptr2
		pptr2 = groupPtr->erase(pptr2);
		delete ptr2;
	      }		
	    } else if (ptr2->includes(*ptr1)) {
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
    
    // Save for each orbit:
    vector<string> startFile; // File holding source orbit
    vector<LONGLONG> startID; // and its id
    vector<int> nDetect; // Number of detections
    vector<int> nUnique; // Number of independent nights
    vector<double> chisq;
    vector<vector<double>> abg;
    vector<vector<double>> abgInvCov;
    vector<vector<double>> elements;  // Elements
    vector<vector<double>> elementCov;
    vector<double> arc;
    vector<double> fpr;
    vector<bool> overlap;
    vector<bool> changed;
    vector<int> friendGroup;
    vector<vector<double>> outFrame; // Reference frame, stored as 7 numbers
    
    // Save for each detection or possible missed detection CCD
    vector<LONGLONG>       orbitID;     // Row number of orbit
    vector<int>            objectID;    // ID from transient table, -1 if no detection
    vector<double>         tdb;         // TDB of exposure
    vector<int>            expnum;
    vector<short int>      ccdnum;      // exposure/ccd of (possible non-)detection
    vector<vector<double>> prediction;  //RA, Dec of predicted position
    vector<vector<double>> predCov;     //ICRS covariance matrix of prediction
    vector<vector<double>> detection;   //RA, Dec of detection (if any)
    vector<vector<double>> detectCov;   //Cov matrix of detection
    vector<vector<double>> residual;    //Det - pred, in ICRS tangent plane
    vector<double>         detectChisq; //Det - pred chisq (meas errors only)
    vector<string>         band;
    vector<double>         mag;
    vector<double>         sn;          //S/N level of flux detection
    
    // Loop through surviving orbits
    int groupCounter = -1;
    for (auto& group : friendGroups) {
      ++groupCounter;
      for (auto& orb : group) {
	// Process each output orbit:
	orb->friendGroup = groupCounter;
	startFile.push_back(orb->inputFile);
	startID.push_back(orb->inputID);
	nDetect.push_back(orb->detectionIDs.size());
	nUnique.push_back(orb->nUnique);
	chisq.push_back(orb->chisq);
	fpr.push_back(orb->fpr);
	vector<double> v(6);
	for (int i=0; i<6; i++) v[i] = orb->abg[i];
	abg.push_back(v);
	for (int i=0; i<6; i++) v[i] = orb->el[i];
	elements.push_back(v);
	v.resize(36);
	for (int i=0; i<6; i++)
	  for (int j=0; j<6; j++) v[i*6+j] = orb->invcov(i,j);
	abgInvCov.push_back(v);
	for (int i=0; i<6; i++)
	  for (int j=0; j<6; j++) v[i*6+j] = orb->elCov(i,j);
	elementCov.push_back(v);
	arc.push_back(orb->arc);
	overlap.push_back(group.size() > 1); // Does orbit overlap others?
	changed.push_back(orb->changedDetectionList);
	v.resize(7);
	orb->frame.orient.getPole().getLonLat(v[0],v[1]);
	v[2] = frame.orient.getPA();
	v[3] = frame.origin.getVector()[0];
	v[4] = frame.origin.getVector()[1];
	v[5] = frame.origin.getVector()[2];
	v[6] = frame.tdb0;
	outFrame.push_back(v);
	
	// save detections and non-detections, transforming into ICRS
	// First make a multimap holding expnum/detection pairs
	std::multimap<int,int> detectionFinder;
	for (int id : orb->detectionIDs) {
	  int expnum = transients.getValue<int>("EXPNUM",id);
	  detectionFinder.emplace(expnum,id);
	}
	// Now loop through all exposures.  Different if it has
	// detections or not.
	for (auto& opp: orb->opportunities) {
	  const Exposure& expo = *(opp.eptr);

	  // Transform orbit prediction and error into
	  // ICRS (cov in local gnomonic, ICRS aligned)
	  astrometry::SphericalICRS radecPred(opp.orbitPred);
	  astrometry::Orientation localICRS(radecPred);
	  astrometry::Gnomonic gn(radecPred, localICRS);
	  Matrix22 covICRS = gn.convertWithCovariance(opp.orbitPred,
						      opp.orbitCov);
	  if (opp.hasDetection) {
	    // Produce one output line per detection
	    auto dptr = detectionFinder.begin(); // declare outside of while
	    while ( (dptr=detectionFinder.find(expo.expnum))!=detectionFinder.end()) {
	      int id = dptr->second;
	      detectionFinder.erase(dptr);

	      expnum.push_back(expo.expnum);
	      tdb.push_back(expo.tdb);
	      band.push_back(et.band(expo.expnum));
	      vector<double> v(2);
	      radecPred.getLonLat(v[0],v[1]);
	      prediction.push_back(v);
	      v.resize(3);
	      v[0] = covICRS(0,0); v[1]=covICRS(1,1); v[2]=covICRS(0,1);
	      predCov.push_back(v);
	      orbitID.push_back(startID.size()-1);  // How big are orbit arrays so far?
	      objectID.push_back(id);
	      ccdnum.push_back( transients.getValue<short int>("CCDNUM",id));
	      auto obs = transients.getObservation(id, ephem, et);
	      v.resize(2);
	      obs.radec.getLonLat(v[0],v[1]);
	      detection.push_back(v);
	      v.resize(3);
	      v[0] = obs.cov(0,0); v[1]=obs.cov(1,1); v[2]=obs.cov(0,1);
	      detectCov.push_back(v);
	      // Calculate residual from prediction in local gnomonic
	      gn.convertFrom(obs.radec);
	      v.resize(2);
	      gn.getLonLat(v[0],v[1]);
	      residual.push_back(v);
	      // And chisq for this point
	      double detCov = obs.cov(0,0)*obs.cov(1,1)-obs.cov(1,0)*obs.cov(0,1);
	      detectChisq.push_back( (v[0]*v[0]*obs.cov(1,1)
				+ v[1]*v[1]*obs.cov(0,0) +
				- 2.*v[1]*v[0]*obs.cov(1,1))/detCov);
	      mag.push_back(transients.getValue<double>("MAG",id));
	      sn.push_back(transients.getValue<double>("FLUX_AUTO",id)
			   / transients.getValue<double>("FLUXERR_AUTO",id) );
	    }
	  } else {
	    // No detection in this exposure.  Output a row for every CCD it might be on.
	    // Possibly none if orbit doesn't really touch a CCD here.
	    for (short int ccd : opp.ccdnums) {
	      expnum.push_back(expo.expnum);
	      tdb.push_back(expo.tdb);
	      band.push_back(et.band(expo.expnum));
	      vector<double> v(2);
	      radecPred.getLonLat(v[0],v[1]);
	      prediction.push_back(v);
	      v.resize(3);
	      v[0] = covICRS(0,0); v[1]=covICRS(1,1); v[2]=covICRS(0,1);
	      predCov.push_back(v);
	      orbitID.push_back(startID.size()-1);  // How big are orbit arrays so far?
	      objectID.push_back(-1); // Negative object ID signals no detection
	      ccdnum.push_back(ccd);
	      detection.push_back(vector<double>(2,0.));
	      detectCov.push_back(vector<double>(3,0.));
	      residual.push_back(vector<double>(2,0.));
	      detectChisq.push_back(0.);
	      mag.push_back(0.);
	      sn.push_back(0.);
	    }
	  }
	} // End loop over opportunities

	// Should not be anything left in the detectionFinder:
	for (auto pr : detectionFinder) {
	  cerr << "#WEIRD: detection " << pr.second
	       << " in exposure " << pr.first
	       << " was not output for orbit " << orb->inputID
	       << endl;
	}

	// Done with this orbit
	delete orb;
      } // close orbit loop
    } // close friend-group loop
    
    // Now write tables to output file
    {
      // First the orbit table - make new file
      FITS::FitsTable ft(outPath, FITS::ReadWrite | FITS::Create | FITS::OverwriteFile, "ORBITS");
      img::FTable table = ft.use();
      table.addColumn(startFile,"OLDFILE");
      table.addColumn(startID,"OLDID");
      table.addColumn(nDetect,"NDETECT");
      table.addColumn(nUnique,"NUNIQUE");
      table.addColumn(chisq,"CHISQ");
      table.addColumn(outFrame,"FRAME");
      table.addColumn(abg,"ABG");
      table.addColumn(abgInvCov,"ABGINVCOV");
      table.addColumn(elements,"ELEMENTS");
      table.addColumn(elementCov,"ELEMENTCOV");
      table.addColumn(arc,"ARC");
      table.addColumn(fpr,"FPR");
      table.addColumn(changed,"CHANGED");
      table.addColumn(overlap,"OVERLAP");
      table.addColumn(friendGroup,"GROUP");
    }
    {
      // Now add detection table to the same file
      FITS::FitsTable ft(outPath, FITS::ReadWrite | FITS::Create, "OBJECTS");
      img::FTable table = ft.use();
      table.addColumn(orbitID,"ORBITID");
      table.addColumn(objectID,"OBJECTID");
      table.addColumn(tdb,"TDB");
      table.addColumn(expnum,"EXPNUM");
      table.addColumn(ccdnum,"CCDNUM");
      table.addColumn(band,"BAND");
      table.addColumn(prediction,"PREDICT");
      table.addColumn(predCov,"PREDCOV");
      table.addColumn(detection,"DETECT");
      table.addColumn(detectCov,"DETCOV");
      table.addColumn(residual,"RESIDUAL");
      table.addColumn(detectChisq,"CHISQ");
      table.addColumn(mag,"MAG");
      table.addColumn(sn,"SN");
    }
    {
      // And an image of the false positive rate
      img::FitsImage<float> fi(outPath,FITS::ReadWrite | FITS::Create, "FP");
      fi.copy(fpImage);
    }
    
  } catch (std::runtime_error& e) {
    quit(e);
  }
  exit(0);
}


    

