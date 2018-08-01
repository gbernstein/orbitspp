// Manipulate information on DECam exposures
#ifndef EXPOSURES_H
#define EXPOSURES_H

#include <list>
#include "Astrometry.h"
#include "OrbitTypes.h"
#include "PlaneGeometry.h"
#include "Ephemeris.h"
#include "Fitter.h"
#include "FitsTable.h"

namespace orbits {

  class Exposure {
    // Information on an exposure.  All quantities in radians internally.
  public:
    int expnum;
    double mjd;
    double tdb;
    double tobs;  // Time relative to Frame tdb0
    Vector3 earth;  // Observatory position in Frame
    Point axis;   // Optical axis coords, in Frame gnomonic
    Matrix22 atmosphereCov;  // Atmospheric astrometric covariance for exposure.

    double detectionDensity; // Density of transients, per sr
    DMatrix xy;   // Transient coordinates, in Frame
    DVector covXX;
    DVector covXY;
    DVector covYY;
    vector<int> ccdnum; // CCD of origin of transients
    vector<int> id;  // ID or row number of transient in original file
    BVector valid;   // Is this a possible new TNO?

    // Return chisq of transients determined
    // by the given error ellipse.  Transients' measurement errors are added in.
    DVector chisq(double x0, double y0, double covxx, double covyy, double covxy) const;
    
    // ??? Add per-CCD info on boundaries
    EIGEN_NEW;
  };

  class ExposureTable {
    // Class that reads information on DES exposures from a standard FITS file
    // of binary tables.  Information available quickly via expnum is MJD and the
    // ICRS position of the telescope at this MJD.
  public:
    ExposureTable(string exposureFile="");
    // Open exposure table.  Path to FITS file will be read from environment variable
    // DES_EXPOSURE_TABLE if none is given here.

    bool isAstrometric(int expnum) const; // Is astrometric solution / covariance available?
    double mjd(int expnum) const;  // Return MJD_MID
    astrometry::CartesianICRS observatory(int expnum) const; //ICRS observatory position
    Matrix22 atmosphereCov(int expnum) const; // Return atmospheric position covariance (radians), Diag(-1) if unknown
    // Get both at once, return false if expnum is not known:
    bool observingInfo(int expnum, double& mjd, astrometry::CartesianICRS& observatory) const;
  private:
    img::FTable astrometricTable;
    img::FTable nonAstrometricTable;
    map<int,int> astrometricIndex;   // These two map from expnum into table row number.
    map<int,int> nonAstrometricIndex;
  };
  
  // This function finds all exposures from the FITS table that could
  // possibly contain a bound TNO that is in the given
  // range of alpha, beta, gamma at the frame reference epoch.
  // Then it loads the transient lists for each
  extern
  std::vector<Exposure*>
  selectExposures(const Frame& frame,   // Starting coordinates, time
		  const Ephemeris& ephem,  // Ephemeris to use
		  double gamma0,        // Center and width of range 
		  double dGamma,        // of gamma to cover
		  double searchRadius,  // Range of starting coords to cover (radians)
		  string transientFile="data/zone029.transients.fits",  // File of transients
		  string exposureFile="data/y4a1.exposure.positions.fits",  // File of exposure data
		  double fieldRadius = 1.1*DEGREE); // Circumscribed field radius (radians)
  // ?? Allow exclusion of filters?? seasons??
  // ?? Add CCD corners, detections ??

  class Node {
    // Node of a k-d tree of Exposures over time and position
  public:
    typedef vector<const Exposure*>::iterator dataIter;
    Node( dataIter _begin, dataIter _end,
	  const Ephemeris& ephem,
	  const Frame& frame);
    ~Node() {
      if (left) delete left;
      if (right) delete right;
    }
    
    enum Coordinate {TIME, X, Y};  // The three axes of the tree space
    
    list<const Exposure*> find(const Fitter& path) const;
    // Return a list of exposures that are close to the path.

    int countNodes() const {
      // Give total count of nodes rooted to this one (inclusive)
      if (left) return left->countNodes() + right->countNodes() + 1;
      else return 1;
    }

    static void setSpeed(double speed_) {speed=speed_;} // Rough max TNO speed, rad/yr
    static void setFieldRadius(double radius_) {fieldRadius = radius_;} // radians
    static void setObservatory(int obscode_) {obscode=obscode_;}
    static Node* buildTree(dataIter begin_, dataIter end_,
			   const Ephemeris& ephem,
			   const Frame& frame);
    // Builds the whole tree out of time-ordered data.
    // Be sure to set speed and radius first.  Returns root node pointer.
    // Needs 
    EIGEN_NEW
    
  private:
    dataIter begin;  // Range of data associated with this node
    dataIter end;
    Node* left;	// Children (nullptrs if leaf node)
    Node* right;
    double tStart;  // Time (relative to reference frame tdb0) of
    double tStop;   // first and last exposure in set
    DMatrix corners;  // Node bounding box, 4x2 matrix with (LL, LU, UU, UL) in rows.

    Coordinate splitOn;  // Which coordinate does the split happen on?
    double splitValue;   // Value for split

    static double speed;  // Typical object speed, used to choose time vs distance split
    static double fieldRadius; // Size of exposure FOV.  Also sets the scale for radius.
    static int obscode;  // Observatory code for earth positions.
    
    DVector tobs; // Time (relative to tdb0) of begin, midpoint, and end of Node's exposures.
    DMatrix earth; // Position of observatory at begin, mid, and end times.
    
    void split(const Ephemeris& ephem,
	       const Frame& frame);
    // Split all the way down to leaves

  };
  
  class DESTree {
    /* Serves as the root of a kD-tree of exposures that we can search using an orbit.
     * Specialized to DES in that input exposure list is divided into survey years
     * first, then the binary-tree Nodes begin below there.
     */
  public:
    DESTree(std::vector<const Exposure*>& exposurePointers,
	    const Ephemeris& ephem,
	    const Frame& frame,
	    double gamma0);
    ~DESTree();
    
    int countNodes() const;
    list<const Exposure*> find(const Fitter& path) const;
    // Return a list of exposures that are close to the path.
  private:
    std::list<Node*> years;
    static std::list<double> tdb_splits;  // tdb's at which we split DES seasons
  };

} // end namespace orbits
#endif 
