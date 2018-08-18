// Program to do cleanup of orbit sets


#include <iostream>

#include "StringStuff.h"
#include "Astrometry.h"
#include "Fitter.h"
#include "Elements.h"
#include "Exposures.h"
#include "Pset.h"
#include "FitsTable.h"
#include "Eigen/StdVector"

using namespace std;
using namespace orbits;

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

struct Detection {
  Exposure* eptr;
  int objectRow;
  Detection(Exposure* eptr_, int objectRow_): eptr(eptr_), objectRow(objectRow_) {}
  Detection(): eptr(nullptr), objectRow(0) {}
  int id() const {return eptr->id[objectRow];}
  bool operator<(const Detection & rhs) const {return id() < rhs.id();}
};

struct FitResult {
  // Everything we need to know about an orbit
  int orbitID;      // Starting orbit ID
  int nIndependent; // Different nights in fit
  double fpr;       // False positive rate

  vector<int> detectionIDs;  // IDs of fitted detections, *ascending*
  // The observations in ICRS coords:
  vector<Observation, Eigen::aligned_allocator<Observation>> obsICRS;  
  double arc;       // Time span from first to last detection

  class Opportunity {
    // Class representing an exposure that could have seen this
    Exposure* eptr;
    SphericalICRS radec;  // Orbit prediction
    double xPred, yPred;  // prediction in frame
    double covXX, covYY, covXY;  // Model covariance in frame
    int ccdnum;
  };
  vector<Opportunity, Eigen::aligned_allocator<Opportunity>> opportunities;

  double chisq;     // Chisq of best fit
  ReferenceFrame frame;  // Frame for its ABG
  ABG abg;          // orbit
  ABGCovariance invcov; // and covariance
  Elements el;
  ElementCovariance elCov;

  int friendGroup;  // Number of its overlap group (-1=loner)
  bool skippedExposures; // True if some exposures were skipped for large errors
  bool hasOverlap;  // True if this shares detections with another FitResult
  bool isSecure;    // No doubt that this is real, monopolize detections.
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

void
FitResult::fitAndFind(const Ephemeris& ephem,
		      double tdb0,
		      TransientTable& transients,
		      ExposureTable& exposures,
		      vector<Exposure*>& pool) {
  // First fit an orbit to the detectionID's currently in the vector.
  // Then search all exposures for any additional transients
  
  // Acquire Observation for each detected transient
  obsICRS.clear();
  for (auto id : detectionIDs) {
    obsICRS.push_back(transients.getObservation(id, ephem, exposures));
  }

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
  // Choose a reference frame for this orbit.  Use provided TDB,
  // and mean RA/Dec.  Get this mean by summing 3d coords
  {
    Vector3 sum(0.);
    for (auto& obs : obsICRS)
      sum += obs.radec.getUnitVector();
    astrometry::SphericalICRS pole(CartesianICRS(sum));
    frame.orient = astrometry::Orientation(pole);
    frame.orient.alignToEcliptic();
    frame.tdb = tdb0;
    frame.origin = ephem.observatory(807, frame.tdb);
  }

  // Execute fit
  Fitter fit(ephem, Gravity::GIANTS);
  for (auto& obs : obsICRS)
    fit.addObservation(obs);
  fit.setFrame(frame);
  // ! No priors on the orbit.
  fit.setLinearOrbit();
  fit.newtonFit();

  // Save fit results
  chisq = fit.getChisq();
  abg = fit.getABG();
  invcov = fit.getInvCovarABG();
  el = fit.getElements();
  elCov = fit.getElementCovariance();

  // Now go fishing for exposures this orbit crosses
  int nExposures = exposures.size();
  DVector tdbAll(nExposures);
  DMatrix earthAll(nExposures,3);  // ICRS observatory position
  DMatrix axisAll(nExposures,3);   // ICRS direction cosines of pointings
  for (int i=0; i<nExposures; i++) {
    tdbAll[i] = exposures[i]->tdb;
    earthAll.row(i) = exposures[i]->earthICRS.getVector().transpose();
    axisAll.row(i) = exposures[i]->localICRS.getPole().getUnitVector().transpose();
  }

  // convert field radius to a maximum chord length between
  // unit vectors
  const double fieldRadius = 1.15*DEGREE; // A bit larger to be complete
  const double maxChordSq = pow(2. * sin(fieldRadius/2.),2.);

  DMatrix target = fit.getTrajectory().observe(tdbAll, earthAll);
  // Find those within radius
  BVector hits = (target - axis).rowwise().squaredNorm().array() < maxChordSq;
  int nHits = hits.array().count();
  
  // Collect information on all candidate exposures (??? repetitive reading...)
  opportunities.clear();
  tdbAll.resize(nHits);
  earthAll.resize(nHits,3);
  for (int i=0; i<hits.size(); i++) {
    if (hits[i]) {
      Opportunity opp;
      opp.eptr = exposures[i];
      tdbAll[opportunities.size()] = exposures[i]->tdb;
      earthAll.row(opportunities.size()) = exposures[i]->earthICRS.getVector().transpose();
      // Load transient and corner information into each exposure
      transients.fillExposure(frame, opp.eptr);
    }
  }
  // ??? Did we get all the exposures from original orbit?
  tdbAll.array() -= frame.tdb0;
  
  // Predict orbit position at each exposure
  DVector xPred(nHits);
  DVector yPred(nHits);
  DVector covxxPred(nHits);
  DVector covyyPred(nHits);
  DVector covxyPred(nHits);
  fit.predict(tdbAll, earthAll, &xPred, &yPred,
	      &covxxPred, &covyyPred, &covxyPred);

  // Save info on every exposure, find new detections
  for (int i=0; i<nHits; i++) {
    Opportunity& opp = opportunities[i];
    opp.xPred = xPred[i];
    opp.yPred = yPred[i];
    opp.covXX = covxxPred[i];
    opp.covYY = covyyPred[i];
    opp.covXY = covxyPred[i];

    // ?? Check in detail for error ellipse crossing a CCD

    opp.radec = astrometry::SphericalICRS(astrometry::Gnomonic(opp.xPred,
							       opp.yPred,
							       frame.orient));

    opp.ccdnum = opp.eptr->whichCCD(opp.radec);
    
  // Find detections - did any change?? - do CCDNUM match prev?
  // Save residuals for each detection
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

  double ra0;    // Center of field (ignore exposures >60 degrees away)
  double dec0;
  double searchRadius;
  double tdb0;   // Reference time (=epoch of orbits or state vectors)
  // Else orbital elements

  const int obscode=807;
  const double fieldRadius = 1.1*DEGREE;
  
  Pset parameters;
   
  try {
    {
      const int def=PsetMember::hasDefault;
      const int low=PsetMember::hasLowerBound;
      const int up=PsetMember::hasUpperBound;
      const int lowopen = low | PsetMember::openLowerBound;
      const int upopen = up | PsetMember::openUpperBound;

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
      parameters.addMember("radius",&searchRadius, def || lowopen || up,
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
    vector<Exposure*> exposures = et.getPool(frame, ephem, searchRadius*DEGREE);
    int nExposures = exposures.size();
    DVector tdb(nExposures);
    DMatrix earth(nExposures,3);  // ICRS observatory position
    DMatrix axis(nExposures,3);   // ICRS direction cosines of pointings
    for (int i=0; i<nExposures; i++) {
      tdb[i] = exposures[i]->tdb;
      earth.row(i) = exposures[i]->earthICRS.getVector().transpose();
      axis.row(i) = exposures[i]->axisICRS.getUnitVector().transpose();
    }

    cerr << "# Read " << exposures.size() << " exposures" << endl;
    
    // Load CCD corners
    {
      CornerTable tt(cornerPath);
      for (auto& eptr : exposures)
	tt.fillExposure(eptr);
    }

    // convert field radius to a maximum chord length between
    // unit vectors
    const double maxChordSq = pow(2. * sin(fieldRadius/2.),2.);

    // Start reading input orbits, giving each one its own
    // reference frame
    list<FitResult*> orbits;

    for (auto orbitPath : orbitFiles) {
      // Open the tables/files

      // Sum the Accumulators
      
      for (int row=0; row<orbitTable.nrows(); row++) {
	// Collect known transients

	// Make a ReferenceFrame

	// Refit orbit, no priors

	// Find all exposure crossings, 

	// (Re-)Find all transients along orbit

	// Re-fit if any new detections
	// Update exposure crossings, CCD, cov
	// Save new chisq
      }
    }

    // Purge duplicates and subsets

    // Establish overlapping groups, assign group numbers
	
    // Loop through surviving orbits {

      // Save orbit info, including elements

      // save detections, get residuals, chisq per point, and band, mag, magerr
    
      // save non-detections
    // }

    // Write orbit and object tables
    
    // Save the summed accumulator

    
    vector<LONGLONG> idOut;
    vector<int> expnumOut;
    vector<int> ccdnumOut;
    vector<double> tdbOut;
    vector<double> raOut;
    vector<double> decOut;
    

    // Begin reading input lines
    string buffer;
    while (stringstuff::getlineNoComment(cin, buffer)) {
      State xv;
      LONGLONG orbitID;
      xv.tdb = tdb0;
      istringstream iss(buffer);
      if (readState) {
	// Read state and orbitID
	iss >> orbitID
	    >> xv.x[0] >> xv.x[1] >> xv.x[2]
	    >> xv.v[0] >> xv.v[1] >> xv.v[2];
      } else {
	// Read elements and orbitID
	Elements el;
	iss >> orbitID >> el;
	// convert to state
	xv = getState(el, tdb0);
      }

      // Predict for all relevant exposures - get ICRS direction cosines
      Trajectory orbit(ephem, xv, Gravity::GIANTS);
      DMatrix target = orbit.observe(tdb, earth);
      // Find those within radius
      BVector hits = (target - axis).rowwise().squaredNorm().array() < maxChordSq;
      
      for (int i=0; i<nExposures; i++) {
	if (!hits[i]) continue;

	astrometry::SphericalICRS radec;
	radec.setUnitVector(target.row(i).transpose());
	double ra,dec;
	radec.getLonLat(ra,dec);

	if (useStdout) {
	  cout << setw(8) << orbitID
	       << " " << setw(6) << exposures[i]->expnum
	       << " " << fixed << setprecision(8) << setw(11) << exposures[i]->tdb
	       << " " << setprecision(7) << setw(11) << ra/DEGREE
	       << " " << showpos << setprecision(7) << setw(11) << dec/DEGREE
	       << " " << noshowpos << setw(2) << exposures[i]->whichCCD(radec)
	       << endl;
	} else {
	  idOut.push_back(orbitID);
	  expnumOut.push_back(exposures[i]->expnum);
	  ccdnumOut.push_back(exposures[i]->whichCCD(radec));
	  tdbOut.push_back(exposures[i]->tdb);
	  raOut.push_back(ra/DEGREE);
	  decOut.push_back(dec/DEGREE);
	}
      }
    } // End the input reading loop.

    if (!useStdout) {
      // Write the results to a FITS file
      FITS::FitsTable ft(positionPath, FITS::Create + FITS::OverwriteFile);
      img::FTable out = ft.use();
      out.addColumn(idOut,"ORBITID");
      out.addColumn(expnumOut,"EXPNUM");
      out.addColumn(ccdnumOut,"CCDNUM");
      out.addColumn(tdbOut,"TDB");
      out.addColumn(raOut,"RA");
      out.addColumn(decOut,"DEC");
    }
  } catch (std::runtime_error& e) {
    quit(e);
  }
  exit(0);
}

    

