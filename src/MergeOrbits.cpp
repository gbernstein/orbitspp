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

#ifdef _OPENMP
#include "omp.h"
#endif

using namespace std;
using namespace orbits;

const int DEBUGLEVEL= 0;
const int MIN_UNIQUE_FIT=4; // min number of independent detections to fit
const double FIELD_RADIUS = 1.15*DEGREE; // A bit larger to be complete

// Time that must pass between exposures to be considered independent detections
// (e.g. when asteroids or defects would have moved out of linking range)
const double INDEPENDENT_TIME_INTERVAL = 0.1*DAY;

const string usage =
  "MergeOrbits: program to combine multiple orbit files (produced by\n"
  "   GrowOrbits or by this program) into a single file with duplicates\n"
  "   removed.  Every input orbit is refit, with minimal priors, to all\n"
  "   detections listed in the transientFile.  The output FITS file\n"
  "   contains two tables and one image.  The ORBITS table contains\n"
  "   one row per orbit that survives the merge.  The OBJECTS table\n"
  "   contains one row per either a detection that is assigned to \n"
  "   an orbit, or an exposure/CCD combination that is within the\n"
  "   error ellipse for a fitted orbit but has no detection on the\n"
  "   exposure, i.e. a missed opportunity.  The FP extension is\n"
  "   an image that holds a 2d table of estimated total number of false\n"
  "   positive detections, summed over all the searches creating the\n"
  "   input files.\n"
  "\n"
  "   As well as re-finding all detections (or non-detections)\n"
  "   consistent with each input orbit and re-fitting the orbits,\n"
  "   this program discards all refit orbits whose detection list\n"
  "   is a duplicate or subset of those in another orbit.  Orbits\n"
  "   which share one or more detections are marked as \"overlap\"\n"
  "   and put into groups by a friends-of-friends algorithm.\n"
  "   [??? add culling of these friend groups later ???]\n"
  "\n"
  "Usage: MergeOrbits <orbitfile1> [orbitfile2] ... [-parameter value]\n"
  "   <orbitfileN> are any number >=1 of outputs from GrowOrbits\n"
  "   [-parameter value] are additional parameter name/value pairs.\n"
  "   Parameters are listed below.\n"
  "\n"
  "Input files formats: Any output from this program or GrowOrbits\n"
  "        is valid input here.  The three extension names described\n"
  "        above must be present.  Only the following columns are\n"
  "        used from the ORBITS table: FPR, NUNIQUE\n"
  "        & used from the OBJECTS table: ORBITID, OBJECTID\n"
  "Output file format: These are the columns in the ORBITS table.\n"
  "Parentheses give the FITS column type.\n"
  "   OLDFILE:    (A) The input file from which the orbit came\n"
  "   OLDID:      (K) The row number of the orbit in that file\n"
  "   FPR:        (D) Estimated false positive rate for finding\n"
  "               the final detection of the orbit in its original\n"
  "               search (no new FPR is calculated during merge; if\n"
  "               two orbits from different searches are merged we\n"
  "               retain the lower FPR)\n"
  "   NDETECT:    (J) Total number of detections found for this orbit\n"
  "   NUNIQUE:    (J) Total number of \"unique\" detections, i.e. not\n"
  "               within a few hours of each other\n"
  "   ARC:        (D) Time interval between first and last detection (yrs)\n"
  "   ARCCUT:     (D) Shortest arc obtained from dropping any one night\n"
  "   CHANGED:    (L) Flag is set if detection list was different when\n"
  "               re-searching the transient catalog for this orbit.\n"
  "   OVERLAP:    (L) Flag is set if this orbit shares detection(s) with\n"
  "               any other surviving orbit.\n"
  "   MAXOVERLAP: (J) Largest NUNIQUE value for any other orbit which\n"
  "               shares a detection with this one.\n"
  "   GROUP:      (J) Integer code that will be common to all orbits\n"
  "               sharing detections with this one (by friends-of-\n"
  "               friends algorithm)\n"
  "   FRAME:      (7D) The reference frame for this orbit's ABG, stored\n"
  "               as an array of seven numbers giving\n"
  "               [RA0, DEC0, PA0, X0, Y0, Z0, TDB0] where first 3\n"
  "               are the spherical coord orientation (rad), second 3\n"
  "               are ICRS origin (in AU), and last is reference TDB\n"
  "               in years past J2000\n"
  "   ABG:        (6D) The best-fit orbit, in alpha/beta/gamma values\n"
  "   ABGINVCOV:  (36D) Inverse covariance matrix of ABG after fit\n"
  "   ELEMENTS:   (6D) Best-fit orbital elements (radians)\n"
  "               [A E I LAN AOP TOP]\n"
  "   ELEMENTCOV: (36D) Covariance matrix of orbital elements.\n"
  "\n"
  "The columns in the OBJECTS table are:\n"
  "   ORBITID:    (K) Row number in ORBITS table of the orbit to which\n"
  "               this (non-) detection belongs.\n"
  "   OBJECTID:   (J) Row number of detection in transient table.\n"
  "               If <0, there is no detection, and this row is\n"
  "               indicating and exposure/ccd pair which is within\n"
  "               the search ellipse for the orbit on an exposure\n"
  "               with no detection. Then OBJECTID is negative of\n"
  "               closest detection to prediction, and CHISQ gives its\n"
  "               distance in a chisq sense\n"
  "   EXPNUM:     (J) Exposure number of the (non-)detection\n"
  "   CCDNUM:     (I) Device (CCD) number of the (non-) detection\n"
  "   TDB:        (D) Time of midpoint of exposure (years past J2000)\n"
  "   BAND:       (A) Filter band name of the observation\n"
  "   PREDICT:    (2D) RA, Dec (radians) of orbit prediction of position\n"
  "   PREDCOV:    (3D) Covariance matrix elements (xx,yy,xy) of \n"
  "               uncertainty ellipse for prediction of orbit fit.\n"
  "               The x axis is to ICRS east, y to N, values in\n"
  "               rad^2.\n"
  "   DETECT:     (2D) RA,Dec of detected object (if any)\n"
  "   DETCOV:     (3D) Measurement covariance matrix for DETECT, same\n"
  "               format as PREDCOV.\n"
  "   RESIDUAL:   (2D) (detected-predicted) position, in radians,\n"
  "               same x/y system.\n"
  "   CHISQ:      (D) The chisq signficance of residual compared to\n"
  "               DETCOV ellipse\n"
  "   MAG:        (D) Magnitude of detection\n"
  "   SN:         (D) S/N ratio of detection (flux / fluxerr)\n"
  "\n"
  "The contents of the FP image are documented elsewhere.  Parameters for this program:\n"
  " -secureUnique  is number of unique nights needed to define an orbit as \"secure\".\n"
  " -secureFPR     is maximum FPR needed to define an orbit as \"secure\".\n"
  " -secureArccut  is minimum arccut needed to define an orbit as \"secure\".\n";
// ?? more parameter descriptions needed...
  

struct FitResult {
  // Everything we need to know about an orbit
  FitResult(): trajectory(nullptr) {};
  ~FitResult() {if (trajectory) delete trajectory;}
  
  string inputFile; // Orbit file it came from
  LONGLONG inputID;      // Starting orbit ID
  int inputUnique; // Different nights in original fit
  double fpr;       // False positive rate in original search

  vector<int> detectionIDs;  // IDs of fitted detections, *ascending*
  // The observations in ICRS coords:
  int nUnique;      // Number of distinct detection times
  vector<Observation, Eigen::aligned_allocator<Observation>> obsICRS;  
  double arc;       // Time span from first to last detection
  double arccut;    // Shortest time span after we cut out any one night.

  Trajectory* trajectory;  // Fitted trajectory to the object
  
  class Opportunity {
  public:
    Opportunity(): eptr(nullptr), orbitCov(0.), hasDetection(false) {}
    // Class representing an exposure that could have seen this
    Exposure* eptr;
    astrometry::Gnomonic orbitPred;  // Orbit prediction (with frame Orientation)
    Matrix22 orbitCov;   // Model covariance in frame
    vector<int> ccdnums;
    bool hasDetection;   // set if there is a detection on the exposure
    double nearestChisq;      // Chisq to model of best detection on exposure.
    int  nearestID;      // object number of nearest
    EIGEN_NEW
  };
  vector<Opportunity, Eigen::aligned_allocator<Opportunity>> opportunities;

  double chisq;     // Chisq of best fit
  Frame frame;       // Frame for its ABG
  ABG abg;          // orbit
  ABGCovariance abgcov; // and (inverse) covariance
  Elements el;
  ElementCovariance elCov;

  int friendGroup;  // Number of its overlap group (-1=loner)
  list<FitResult*> *itsGroup; // Will point to group of intersecting orbits
  int mostPopulousOverlap; // Highest nUnique for any orbit it shares detection with
  
  bool changedDetectionList; // set if detections were added / dropped
  bool secure;      // Object is sure enough to claim the detections
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

  bool testSecure() {
    // Does this meet the criteria for a "secure" orbit?
    // Appropriately sets the boolean flag and returns its value
    // Note need to set the class variables that define this test
    secure = nUnique >= secureUnique
             && fpr <= secureFPR
             && arccut >= secureArccut;
    return secure;
  }

  static void setSecure(int secureUnique_, double secureFPR_, double secureArccut_) {
    // Set the class-wide criteria for secure detection
    secureUnique = secureUnique_;
    secureFPR = secureFPR_;
    secureArccut = secureArccut_;
  }
  static void setSearchChisq(double chisq) {
    searchChisq = chisq;
  }
private:
  static int secureUnique;
  static double secureFPR;
  static double secureArccut;
  static double searchChisq;  // Maximum chisq to consider a match in 2d
};

int FitResult::secureUnique = 8;
double FitResult::secureFPR = 0.001;
double FitResult::secureArccut = 0.7;
double FitResult::searchChisq = 9.;

bool
FitResult::fitAndFind(const Ephemeris& ephem,
		      double tdb0,
		      TransientTable& transients,
		      ExposureTable& exposureTable,
		      vector<Exposure*>& pool) {
  // First fit an orbit to the detectionID's currently in the vector.
  // Then search all exposures for any additional transients.
  // Return true if any detections changed.

  // Acquire Observation for each detected transient
  obsICRS.clear();
  if (DEBUGLEVEL>1)
#pragma omp critical (io)
    cerr << "#" << inputID << "  fitAndFind: Reading detection info "
	 << inputID << " has " << detectionIDs.size()
	 << " " << nUnique << endl;
  for (auto id : detectionIDs) {
    obsICRS.push_back(transients.getObservation(id, ephem, exposureTable));
  }

  if (DEBUGLEVEL>1)
#pragma omp critical (io)
    cerr << "#" << inputID << "  ...Calculating arc length" << endl;
  // Calculate arc, and arc after dropping one night
  {
    double t1 = obsICRS.front().tdb;
    double t2 = t1;
    for (auto& obs : obsICRS) {
      if (obs.tdb < t1) t1 = obs.tdb;
      if (obs.tdb > t2) t2 = obs.tdb;
    }
    arc = t2 - t1;
    // Now go back and find second-highest/lowest
    double t1a = t2;
    double t2a = t1;
    for (auto& obs : obsICRS) {
      if (obs.tdb-t1 > 0.5*DAY &&  obs.tdb < t1a) t1a = obs.tdb;
      if (t2-obs.tdb > 0.5*DAY &&  obs.tdb > t2a) t2a = obs.tdb;
    }
    arccut = MIN(t2-t1a, t2a-t1);
  }

  if (DEBUGLEVEL>1)
#pragma omp critical (io)
    cerr << "#" << inputID << "  ...reference frame" << endl;

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
#pragma omp critical (io)
    cerr << "#" << inputID << "  ...fitting in frame at "
	 << frame.orient.getPole() << endl;

  // Execute fit
  Fitter fit(ephem, Gravity::GIANTS);
  for (auto& obs : obsICRS)
    fit.addObservation(obs);
  fit.setFrame(frame);
  // ! No priors on the orbit.
  // ?? maybe a little
  fit.setBindingConstraint(4.); // ?? stronger binding  than 1.
  fit.setLinearOrbit();
  fit.setLinearOrbit();
  if (DEBUGLEVEL>1)
#pragma omp critical (io)
    {
      cerr << "#" << inputID << "  linear chisq " << fit.getChisq() << " abg " << fit.getABG(true)
	   << " a " << fit.getElements()[Elements::A] << endl;
    }
    fit.newtonFit();

  if (DEBUGLEVEL>1)
#pragma omp critical (io)
    cerr << "#" << inputID << "  ...Fit complete" << endl;

  // Save fit results
  chisq = fit.getChisq();
  abg = fit.getABG();
  abgcov = fit.getInvCovarABG();
  el = fit.getElements();
  elCov = fit.getElementCovariance();
  trajectory = new Trajectory(fit.getTrajectory());
  
  if (DEBUGLEVEL>1)
#pragma omp critical (io)
    cerr << "#" << inputID << "  ...chisq " << chisq << " abg " << abg
	 << " a,e " << el[Elements::A] << "," << el[Elements::E] << endl;

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
  const double maxChordSq = pow(2. * sin(FIELD_RADIUS/2.),2.);

  DMatrix target = fit.getTrajectory().observe(tdbAll, earthAll);
  // Find those within radius
  BVector hits = (target - axisAll).rowwise().squaredNorm().array() < maxChordSq;
  int nHits = hits.array().count();
  
  if (DEBUGLEVEL>1)
#pragma omp critical (io)
    cerr << "#" << inputID << "  ...exposure search has "
	 << nHits << " hits" << endl;

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
#pragma omp critical (io)
    cerr << "#" << inputID << "  ...getting predictions" << endl;

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
#pragma omp critical (io)
    cerr << "#" << inputID << "  ...finding detections" << endl;

  for (int i=0; i<nHits; i++) {
    Opportunity& opp = opportunities[i];
    opp.orbitPred = astrometry::Gnomonic(xPred[i], yPred[i], frame.orient);
    opp.orbitCov(0,0) = covxxPred[i];
    opp.orbitCov(0,1) = covxyPred[i];
    opp.orbitCov(1,0) = covxyPred[i];
    opp.orbitCov(1,1) = covyyPred[i];

    // Check in detail for error ellipse crossing a CCD
    try {
      opp.ccdnums = opp.eptr->whichCCDs(opp.orbitPred,
					opp.orbitCov*searchChisq);
    } catch (astrometry::AstrometryError &m) {
      // Catch prediction that has moved out of hemisphere
#pragma omp critical (io)
      cerr << "WARNING: Orbit " << inputID << " AstrometryError" << endl;
      opp.nearestID = 0;
      opp.nearestChisq = 0;
      continue;  // Move on to next exposure
    }

    // Load transient information into each exposure, if not done already
    transients.fillExposure(opp.eptr);

    DVector allChi = opp.eptr->chisq(frame,
				     xPred[i], yPred[i],
				     covxxPred[i], covyyPred[i], covxyPred[i]);
    if (DEBUGLEVEL>2) {
      double t = 0.5*(covxxPred[i]+covyyPred[i]);
      double e = hypot(0.5*(covxxPred[i]-covyyPred[i]),covxyPred[i]);
      cerr << "...exposure " << opp.eptr->expnum << " a/b "
	   << setprecision(3) << sqrt(t+e)/ARCSEC << "/" << sqrt(t-e)/ARCSEC
	   << " finds: ";
    }
    // We will record nature of the closest transient on the exposure
    if (allChi.size()>0) {
      opp.nearestID = opp.eptr->id[0];
      opp.nearestChisq = allChi[0];
    } else {
      opp.nearestID = 0;
      opp.nearestChisq = 0;
    }
    for (int iTrans=0; iTrans<allChi.size(); iTrans++) {
      if (allChi[iTrans] < opp.nearestChisq) {
	opp.nearestID = opp.eptr->id[iTrans];
	opp.nearestChisq = allChi[iTrans];
      }
      
      if (allChi[iTrans]<searchChisq && opp.eptr->valid[iTrans]) {
	opp.hasDetection = true;
	newDetectionIDs.push_back(opp.eptr->id[iTrans]);
	if (DEBUGLEVEL>2) cerr << opp.eptr->id[iTrans] << " ";
	// Is this a unique detection?
	bool isUnique = true;
	for (auto tdb : timesOfDetections) {
	  if ( abs(tdb-tdbAll[i]) < INDEPENDENT_TIME_INTERVAL) {
	    isUnique = false;
	    break;
	  }
	}
	if (isUnique) {
	  timesOfDetections.push_back(tdbAll[i]);
	  nUnique++;
	}
      }
    }
    if (DEBUGLEVEL>2)
      cerr << " nearest " << opp.nearestID << "@" << opp.nearestChisq << endl;

  }

  // Compare new and old detection lists
  std::sort(newDetectionIDs.begin(), newDetectionIDs.end());
  std::sort(detectionIDs.begin(), detectionIDs.end());

  bool changed = (newDetectionIDs!=detectionIDs);

  if (DEBUGLEVEL>2) 
#pragma omp critical(io)
    {
    set<int> s1;
    s1.insert(detectionIDs.begin(), detectionIDs.end());
    set<int> s2;
    s2.insert(newDetectionIDs.begin(), newDetectionIDs.end());
    cerr << "#" << inputID << "  ...detections lost: ";
    set<int> result;
    set_difference(s1.begin(), s1.end(), s2.begin(), s2.end(),
		   std::inserter(result, result.end()));
    for (auto id: result) cerr << " " << id;
    cerr << endl;
    cerr << "#" << inputID << "  ...detections gained: ";
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
  if (DEBUGLEVEL>1)
#pragma omp critical(io)
    cerr << "#" << inputID << " FitAndFind complete" << changed << endl;
  return changed;
}

class OrbitFactory {
public:
  OrbitFactory(const vector<string>& orbitFiles_): fileNumber(-1),
						   orbitFiles(orbitFiles_) {}
  FitResult* next() {
    // Return pointer to next orbit to analyze, nullptr if we
    // are done.

    // Open a new file if needed
    while (fileNumber<0 || orbitTableRow >= inputOrbitTable.nrows()) {
      fileNumber++;
      if (fileNumber >= orbitFiles.size()) {
	// No more input files, done.
	return nullptr;
      }
      orbitPath = orbitFiles[fileNumber];
      cerr << "# Reading orbits from " << orbitPath << endl;
      // Open the tables/files
      {
        FITS::FitsTable ft(orbitPath, FITS::ReadOnly,"ORBITS");
	inputOrbitTable = ft.extract();
      }
      {
        FITS::FitsTable ft(orbitPath, FITS::ReadOnly, "OBJECTS");
	inputObjectTable = ft.extract();
      }
      // Reset the row counters
      objectTableRow = 0;
      orbitTableRow = 0;
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
    }

    // Now acquire the next row from the table
    auto orb = new FitResult;
    orb->inputFile = orbitPath;
    orb->inputID = orbitTableRow;
    inputOrbitTable.readCell(orb->inputUnique, "NUNIQUE", orbitTableRow);
    orb->nUnique = orb->inputUnique;
    inputOrbitTable.readCell(orb->fpr, "FPR", orbitTableRow);
    ++orbitTableRow;
    // Collect known transients for this orbit.
    while (objectTableRow < inputObjectTable.nrows()) {
      LONGLONG id;
      inputObjectTable.readCell(id, "ORBITID", objectTableRow);
      if (id!=orb->inputID)
	break;  // This row belongs to another orbit
      int objectID;
      inputObjectTable.readCell(objectID, "OBJECTID", objectTableRow);
      if (objectID>=0) {
	// negative ID is non-detection
	orb->detectionIDs.push_back(objectID);
      }
      ++objectTableRow;
    }
    return orb;
  }

  img::Image<float> getFPR() const {return fpImage;}
						 
private:
  int fileNumber;
  const vector<string>& orbitFiles;
  img::FTable inputOrbitTable;
  img::FTable inputObjectTable;
  string orbitPath;
  int objectTableRow;
  int orbitTableRow;
  img::Image<float> fpImage; // Summed false-positive image
};

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

  int minUnique;  // Min number unique nights to retain orbit for output and FoF
  double maxFPR;  // Max (initial) FPR value to retain
  // Parameters defining a secure orbit:
  int secureUnique;
  double secureFPR;
  double secureArccut;
  double searchChisq;
  // ??? make unique/FPR pairs

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
      parameters.addMember("radius",&radius, def | lowopen | up,
			   "Radius of search area (deg)", 85.,0.,85.);
      parameters.addMember("minUnique",&minUnique, def | low,
			   "Minimum unique nights in output orbit", 5,MIN_UNIQUE_FIT);
      parameters.addMember("maxFPR",&maxFPR, def | lowopen,
			   "Maximum FPR input orbit to use", 0.03,0.);
      parameters.addMember("secureUnique",&secureUnique, def | low,
			   "minimum unique nights for secure detection", 8, 5);
      parameters.addMember("secureFPR",&secureFPR, def,
			   "maximum FPR for secure detection", 0.001);
      parameters.addMember("secureArccut",&secureArccut, def | low,
			   "minimum Arccut for secure detection", 0.7);
      parameters.addMember("searchChisq",&searchChisq, def | low,
			   "max 2d chisq for match to orbit", 9., 4.);
      
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
    ephem.cacheGiants();

    //** Preload cache for giants to
    // avoid a bug in the SharedLUT
    for (auto body : {SUN, JUPITER, SATURN, URANUS, NEPTUNE}) {
      // Ask for positions at starts of 2012, 2020
      ephem.position(body, 12.);
      ephem.position(body, 20.);
    }
    cerr << "# Primed ephemeris cache" << endl;
    /**/
    // Establish a reference frame
    Frame frame;
    {
      astrometry::SphericalICRS pole(ra0*DEGREE, dec0*DEGREE);
      frame.orient = astrometry::Orientation(pole);
      frame.orient.alignToEcliptic();
      frame.tdb0 = tdb0;
      frame.origin = ephem.observatory(807, tdb0);
    }

    // Set up the definition of secure orbit
    FitResult::setSecure(secureUnique, secureFPR, secureArccut);

    FitResult::setSearchChisq(searchChisq);

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
    
    // Create the orbit factory
    OrbitFactory factory(orbitFiles);

    // Start reading input orbits, giving each one its own
    // reference frame
    list<FitResult*> orbits;

    cerr << "# Starting linking" << endl;

    // Start an outer loop in which we'll queue up a pile of
    // orbits to distribute.
#ifdef _OPENMP
    int blocksize = omp_get_max_threads() * 32;
#else
    int blocksize = 32;
#endif

    while (true) {
      vector<FitResult*> blockOfOrbits;
      do {
	auto orb = factory.next();
	if (orb==nullptr) 
	  break;  // No more orbits to acquire

	// Some initial quality checks; first FPR
	if (orb->fpr > maxFPR) {
	  // Don't use this one.
	  //**/cerr << "Deleting for FPR" << endl;
	  delete orb;
	  continue;
	}

	// Also require a minimum number of detections to proceed
	if (orb->detectionIDs.size() < 4) {
	  //**/cerr << "Deleting for size" << endl;
	  delete orb;
	  continue;
	}

	blockOfOrbits.push_back(orb);
      } while (blockOfOrbits.size() < blocksize);

      int imax = blockOfOrbits.size();
      /**/cerr << "Acquired block of " << imax << " input orbits" << endl;
      /**/cerr << "   Kept so far: " << orbits.size() << endl;
      if (imax==0)
	break;  // Done with all input!

      // Now distribute the orbits among processors
#pragma omp parallel for schedule(dynamic,1)
      for (int i=0; i<imax; i++) 
      {
	auto orb = blockOfOrbits[i];
	/***
#pragma omp critical(io)
	cerr << ">Thread " << omp_get_thread_num()
	     << " starting " << orb->inputID
	     << " i " << i
	     << endl;
	/**/
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
	    if (orb->nUnique < MIN_UNIQUE_FIT)  {
	      break;
	    }
	    if (iter==MAX_FIT_ITERATIONS-1)
	      cerr << "# WARNING: Not converging at " << orb->inputID << endl;
	  }
	} catch (Fitter::NonConvergent& e) {
#pragma omp critical(io)
	  {
	    cerr << "# Fitting failure " << orb->inputFile << "/" << orb->inputID
		 << " nUnique " << orb->nUnique << " arc " << orb->arc << endl;
	    cerr << e.what() << endl;
	  }
	}
	if (success)
#pragma omp critical(orbitpush)
	  {
	    orbits.push_back(orb);
	  }
	else {
	    delete orb;
	}
	  
      } // End orbit loop (and parallel section)
      //**/cerr << "Completed orbit block" << endl;
    } // End block loop

    cerr << "# Processed " << orbits.size() << " orbits" << endl;

    // Determine and count secure orbits
    set<int> secureDetections;
    int nSecure = 0;
    for (auto ptr : orbits) {
      if (ptr->testSecure()) {
	secureDetections.insert(ptr->detectionIDs.begin(),
				ptr->detectionIDs.end());
	++nSecure;
      }
    }
      
    cerr << "# Found " << nSecure
	 << " secure orbits with " << secureDetections.size()
	 << " distinct detections" << endl;
    
    // Sweep through all orbits, deleting those that
    // are too small or share detections with secure
    // orbits.

    // This multimap will contain all detection/orbit 
    // pairs.  Sorting by detection id's eliminates the
    // need to compare all pairs of orbits/detections
    // to make the subset and friend identifications.
    // Once this map is constructed, it becomes
    // the repository of all valid FitResult pointers.

    class Obj2Orbit: public multimap<int, FitResult*> {
    public:
      typedef multimap<int,FitResult*>::iterator iter;
      // Add all the detections of a FitResult to the map
      void add(FitResult* fptr) {
	for (auto id : fptr->detectionIDs)
	  emplace(id,fptr);
      }
      // Remove all detections from FitResult in this
      // element, delete the FitResult,
      // and return iterator to next element of the multimap.
      multimap<int,FitResult*>::iterator remove(iter killit) {
	auto out = killit;
	++out;  // Save this as we will soon delete this guy
	FitResult* fptr = killit->second;
	for (auto id : fptr->detectionIDs) {
	  // Search all map elements having this ID as key
	  // for the one we want to kill.
	  auto pr = this->equal_range(id);
	  for (iter ptr = pr.first;
	       ptr != pr.second;
	       ++ptr) {
	    if (ptr->second == fptr) {
	      // This is the one to delete
	      // Be careful if it was going to be the output
	      if (ptr == out)
		out = this->erase(ptr);
	      else
		this->erase(ptr);
	      break;
	    }
	  }
	}
	delete fptr;
	return out;
      }
    } obj2orbit;
    
    cerr << "# Purging insufficient orbits and shards" << endl;
    int nLeft = 0;
    for (auto ptr : orbits) {
      if (ptr->nUnique < minUnique) { 
	// Throw it away
	delete ptr;
	continue;
      }

      // is this a mis-track of a secure orbit?
      if (!ptr->secure) {
	bool kill = false;
	for (auto id : ptr->detectionIDs)
	  if (secureDetections.count(id)>0) {
	    kill = true;
	    break;
	  }
	if (kill) {
	  delete ptr;
	  continue;
	}
      }

      nLeft++;
      obj2orbit.add(ptr); // Valid orbit, we will continue with it.
      ptr->itsGroup = nullptr;  // Prepare for grouping
      ptr->mostPopulousOverlap = 0; // Will track quality of friends
    }
    secureDetections.clear();  // Done with this list.
    
    cerr << "# Remaining orbits: " << nLeft << endl;
    
    // Purge subsets or duplicates.  We only bother
    // to compare orbits that share a detection.
    cerr << "# Starting purge of duplicates/subsets" << endl;
    for (auto ptr1 = obj2orbit.begin();
	 ptr1 != obj2orbit.end(); ) {
      // Compare to all orbits sharing detection,
      // which will be adjacent in multimap
      auto ptr2 = ptr1;
      ++ptr2;
      bool kill1 = false; // Set if we should delete ptr1's orbit
      while (ptr2!=obj2orbit.end() && ptr2->first==ptr1->first) {
	// These orbits share a detection.  Compare.
	if (*(ptr1->second)==*(ptr2->second)) {
	  // If both are same orbit, keep lowest FPR
	  if (ptr1->second->fpr > ptr2->second->fpr) {
	    // Delete ptr1
	    kill1 = true;
	    break;
	  } else {
	    // Delete ptr2
	    ptr2 = obj2orbit.remove(ptr2);
	    nLeft--;
	  }
	} else if (ptr2->second->includes(*(ptr1->second))) {
	  // ptr1 is a subset, delete it, no need to continue
	  kill1 = true;
	  break;
	} else if (ptr1->second->includes(*(ptr2->second))) {
	  // ptr2 is a subset, delete it
	  ptr2 = obj2orbit.remove(ptr2);
	  nLeft--;
	} else {
	  // Neither contains the other, proceed.
	  ++ptr2;
	}
      }
      if (kill1) {
	ptr1 = obj2orbit.remove(ptr1);
	nLeft--;
      } else {
	++ptr1;
      }
    }
    
    cerr << "# Remaining orbits: " << nLeft << endl;
    
    // Now we group orbits into friends-of-friends.
    cerr << "# Beginning FoF grouping" << endl;
    typedef list<FitResult*> Group;

    // Maintain a set of all the Groups.
    set<Group*> allGroups;

    for (auto ptr1 = obj2orbit.begin();
	 ptr1 != obj2orbit.end(); ) {
      // If this orbit is not in a friend group, make one
      // for it:
      if (!ptr1->second->itsGroup) {
	ptr1->second->itsGroup = new Group;
	ptr1->second->itsGroup->push_back(ptr1->second);
	allGroups.insert(ptr1->second->itsGroup);
      }
      auto destGroup = ptr1->second->itsGroup;  // put all friends here
      // Loop over all other orbits sharing this detection
      auto ptr2 = ptr1;
      auto endptr = obj2orbit.upper_bound(ptr1->first);
      for (++ptr2; ptr2!=endptr; ++ptr2) {
	// Any two overlaps raise each other's maxPopulousOverlap
	ptr1->second->mostPopulousOverlap =
	  std::max(ptr1->second->mostPopulousOverlap, ptr2->second->nUnique);
	ptr2->second->mostPopulousOverlap =
	  std::max(ptr2->second->mostPopulousOverlap, ptr1->second->nUnique);

	// Now merge friends groups if different
	auto srcGroup = ptr2->second->itsGroup;
	if (!srcGroup) {
	  // ptr2 is not yet a member of a group, just add it here
	  ptr2->second->itsGroup = destGroup;
	  destGroup->push_back(ptr2->second);
	} else if (srcGroup!=destGroup) {
	  // absorb the srcGroup members
	  for (auto ptr3 : *srcGroup) {
	    ptr3->itsGroup = destGroup; // Redirect each
	  }
	  destGroup->insert(destGroup->end(), srcGroup->begin(), srcGroup->end());
	  // Source group no longer exists
	  allGroups.erase(srcGroup);
	}
      }
      ptr1 = endptr;
    }
    
    cerr << "# Creating outputs" << endl;

    // Now all orbits are grouped into friend sets.
    // Do any further purging within overlap Groups ??
    
    // Save for each orbit:
    vector<string> startFile; // File holding source orbit
    vector<LONGLONG> startID; // and its id
    vector<int> nDetect; // Number of detections
    vector<int> nUnique; // Number of independent nights
    vector<double> chisq;
    vector<vector<double>> abg;
    vector<vector<double>> abgCov;
    vector<vector<double>> elements;  // Elements
    vector<vector<double>> elementCov;
    vector<double> arc;
    vector<double> arccut;
    vector<double> fpr;
    vector<bool> overlap;
    vector<int> maxOverlap;
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
    vector<double>         distance;    //Observer-target distance
    vector<double>         sunDistance; //Sun-target distance
    vector<double>         phase;       //Illumination phase
    vector<double>         phasePA;     // PA of illuminated hemisphere
    
    // Loop through surviving orbits to collect outputs
    int groupCounter = -1;
    for (auto& group : allGroups) {
      ++groupCounter;
      for (auto orb : *group) {
	// Process each output orbit:
	orb->friendGroup = groupCounter;
	friendGroup.push_back(orb->friendGroup);
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
	  for (int j=0; j<6; j++) v[i*6+j] = orb->abgcov(i,j);
	abgCov.push_back(v);
	for (int i=0; i<6; i++)
	  for (int j=0; j<6; j++) v[i*6+j] = orb->elCov(i,j);
	elementCov.push_back(v);
	arc.push_back(orb->arc);
	arccut.push_back(orb->arccut);
	overlap.push_back(group->size() > 1); // Does orbit overlap others?
	maxOverlap.push_back(orb->mostPopulousOverlap);
	changed.push_back(orb->changedDetectionList);
	v.resize(7);
	orb->frame.orient.getPole().getLonLat(v[0],v[1]);
	v[2] = orb->frame.orient.getPA();
	v[3] = orb->frame.origin.getVector()[0];
	v[4] = orb->frame.origin.getVector()[1];
	v[5] = orb->frame.origin.getVector()[2];
	v[6] = orb->frame.tdb0;
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

	  // Get the phase angles and distances for this opportunity
	  Circumstances cc = orb->trajectory->getCircumstances(expo.tdb,
							       expo.earthICRS.getVector());
	  
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
	      // Add observing circumstances
	      distance.push_back(cc.obsDistance);
	      sunDistance.push_back(cc.sunDistance);
	      phase.push_back(cc.phaseAngle);
	      phasePA.push_back(cc.phasePA);
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
	      ccdnum.push_back(ccd);
	      // Make negative object number for nearest detection, give its chisq to orbit
	      if (opp.nearestID>0) {
		objectID.push_back(-opp.nearestID); // Negative object ID signals no detection
	      } else {
		// If the nearest was zero or there were none, put this in:
		objectID.push_back(-1); // I know, this could mean nearest is ID=1...
	      }
	      detectChisq.push_back(opp.nearestChisq);
	      detection.push_back(vector<double>(2,0.));
	      detectCov.push_back(vector<double>(3,0.));
	      residual.push_back(vector<double>(2,0.));
	      mag.push_back(0.);
	      sn.push_back(0.);

	      // Add observing circumstances
	      distance.push_back(cc.obsDistance);
	      sunDistance.push_back(cc.sunDistance);
	      phase.push_back(cc.phaseAngle);
	      phasePA.push_back(cc.phasePA);
	      
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
    cerr << "# Write table of " << startFile.size() << " orbits " << endl;
    {
      // First the orbit table - make new file
      FITS::FitsTable ft(outPath, FITS::ReadWrite | FITS::Create | FITS::OverwriteFile, "ORBITS");
      img::FTable table = ft.use();
      table.addColumn(startFile,"OLDFILE");
      table.addColumn(startID,"OLDID");
      table.addColumn(nDetect,"NDETECT");
      table.addColumn(nUnique,"NUNIQUE");
      table.addColumn(chisq,"CHISQ");
      table.addColumn(outFrame,"FRAME",7);
      table.addColumn(abg,"ABG",6);
      table.addColumn(abgCov,"ABGINVCOV",36);
      table.addColumn(elements,"ELEMENTS",6);
      table.addColumn(elementCov,"ELEMENTCOV",36);
      table.addColumn(arc,"ARC");
      table.addColumn(arccut,"ARCCUT");
      table.addColumn(fpr,"FPR");
      table.addColumn(changed,"CHANGED");
      table.addColumn(overlap,"OVERLAP");
      table.addColumn(maxOverlap,"MAXOVERLAP");
      table.addColumn(friendGroup,"GROUP");
    }
    cerr << "# Write table of " << orbitID.size() << " detections " << endl;
    {
      // Now add detection table to the same file
      FITS::FitsTable ft(outPath, FITS::ReadWrite | FITS::Create, "OBJECTS");
      img::FTable table = ft.use();
      table.addColumn(orbitID,"ORBITID");
      table.addColumn(objectID,"OBJECTID");
      table.addColumn(tdb,"TDB");
      table.addColumn(expnum,"EXPNUM");
      table.addColumn(ccdnum,"CCDNUM");
      table.addColumn(band,"BAND",1,1);
      table.addColumn(prediction,"PREDICT",2);
      table.addColumn(predCov,"PREDCOV",3);
      table.addColumn(detection,"DETECT",2);
      table.addColumn(detectCov,"DETCOV",3);
      table.addColumn(residual,"RESIDUAL",2);
      table.addColumn(detectChisq,"CHISQ");
      table.addColumn(mag,"MAG");
      table.addColumn(sn,"SN");
      table.addColumn(distance,"RANGE");
      table.addColumn(sunDistance,"SOLARD");
      table.addColumn(phase,"PHASE");
      table.addColumn(phasePA,"PHASEPA");
    }
    cerr << "# Write FP image" <<  endl;
    {
      // And an image of the false positive rate
      img::FitsImage<float> fi(outPath,FITS::ReadWrite | FITS::Create, "FP");
      fi.copy(factory.getFPR());
    }

    // Clean up.
    // Drain the exposure pool
    for (auto& ptr : pool) delete ptr;
    
  } catch (std::runtime_error& e) {
    quit(e);
  }
  exit(0);
}


    

