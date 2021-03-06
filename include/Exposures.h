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

#ifdef _OPENMP
#include "omp.h"
#endif

namespace orbits {

  class Exposure {
    // Information on an exposure.  All quantities in radians internally.
  public:
    Exposure(): isFilled(false) {
#ifdef _OPENMP
      omp_init_lock(&fillLock);
#endif
    }
    ~Exposure() {
#ifdef _OPENMP
      omp_destroy_lock(&fillLock);
#endif
    }

    int expnum;
    double mjd;
    double tdb;
    double tobs;  // Time relative to Frame tdb0
    bool astrometric;  // True if exposure is part of astrometric solution

    // Coordinate system centered on optic axis, aligned to ICRS.:
    astrometry::Orientation localICRS; 
    astrometry::CartesianICRS earthICRS; // Position of observatory
    Vector3 earth;  // Observatory position in Frame
    Point axis;   // Optical axis coords, in Frame gnomonic
    Matrix22 atmosphereCov;  // Atmospheric astrometric covariance for exposure.

    // Info on transients:
    double detectionDensity; // Density of transients, per sr
    DMatrix icrs;  // 3xN matrix of ICRS coordinates on unit sphere
    DVector covXX; // cov matrix of detections in an exposure-centered Frame
    DVector covXY;
    DVector covYY;
    vector<int> ccdnum; // CCD of origin of transients
    vector<int> id;  // ID or row number of transient in original file
    BVector valid;   // Is this a possible new TNO?

    
    // map icrs positions into some frame.  Return the two vectors filled with the
    // coordinates of all detections.  
    void mapToFrame(const Frame& frame, DVector& x, DVector& y) const;

    // Use this to set internal xy to values for a frame
    void setFrameXY(const Frame& frame);
    DMatrix xy;   // Transient coordinates, in some user-defined Frame

    // Return chisq of transients determined
    // by the given error ellipse.  Transients' measurement errors are added in.
    DVector chisq(double x0, double y0, double covxx, double covyy, double covxy) const;
    // This one specifies the Frame in which x0,y0 are given, if it may not
    // be whichever one was used to set up xy array.
    DVector chisq(const Frame& frame,
		  double x0, double y0, double covxx, double covyy, double covxy) const;
    
    // Info on CCDs:
    vector<int> devices; // List of ccdnums of devices
    // Device outlines, given in Gnomonic coords about axis aligned to RA/Dec
    vector<ConvexPolygon> deviceBoundsLocalICRS; 
    
    // Return ccdnum of device containing the coordinates. 0 if none.
    int whichCCD(const astrometry::SphericalICRS& radec) const;
    // Return vector of ccdnums that the specified error ellipse
    // touches.  Input coordinates come in some other system and
    // covariance matrix is assumed to be in this system too.
    vector<int> whichCCDs(const astrometry::SphericalCoords& posn,
			   const Matrix22& cov) const;
    EIGEN_NEW;

  private:
    friend class TransientTable;
    bool isFilled;  //
#ifdef _OPENMP
    omp_lock_t fillLock;  // Lock for filling detection structures
#endif
  };

  class ExposureTable {
    // Class that reads information on DES exposures from a standard FITS file
    // of binary tables.  Information available quickly via expnum is MJD and the
    // ICRS position of the telescope at this MJD.
  public:
    ExposureTable(string exposureFile=""); // ?? astrometricOnly=false);
    // Open exposure table.  Path to FITS file will be read from environment variable
    // DES_EXPOSURE_TABLE if none is given here.

    int size() const {return astrometricTable.nrows() + nonAstrometricTable.nrows();}
    int expnum(int index) const;     // Return exposure number at index
    bool isAstrometric(int expnum) const; // Is astrometric solution / covariance available?
    string band(int expnum) const;
    double mjd(int expnum) const;  // Return MJD_MID
    astrometry::CartesianICRS observatory(int expnum) const; //ICRS observatory position
    astrometry::SphericalICRS axis(int expnum) const; // Exposure pointing
    // Get all 3 of above at once, return false if expnum is not known:
    bool observingInfo(int expnum, double& mjd, astrometry::CartesianICRS& observatory,
		       astrometry::SphericalICRS& axis) const;
    Matrix22 atmosphereCov(int expnum) const; // Return atmospheric position covariance (radians), Diag(-1) if unknown

    // This method returns all exposures from the FITS table that could
    // possibly contain a bound TNO that is in the given
    // range of gamma and within searchRadius of frame origin at the frame reference epoch.
    std::vector<Exposure*> getPool(const Frame& frame,
				   const Ephemeris& ephem,
				   double gamma0,        // Center and width of range 
				   double dGamma,        // of gamma to cover
				   double searchRadius,  // Range of starting coords to cover (radians)
				   bool astrometricOnly = true,  // Only look for astrometric exposures
				   double fieldRadius = 1.1*DEGREE) const; // Circumscribed FOV radius (radians)
    // ?? Allow exclusion of filters?? seasons??

    // This version simply returns all exposures centered within a given distance
    // of the frame origin.
    std::vector<Exposure*> getPool(const Frame& frame,
				   const Ephemeris& ephem,
				   double searchRadius,  // Range of starting coords to cover (radians)
				   bool astrometricOnly = true) const;  // Only look for astrometric exposures

  private:
    img::FTable astrometricTable;
    img::FTable nonAstrometricTable;
    map<int,int> astrometricIndex;   // These two map from expnum into table row number.
    map<int,int> nonAstrometricIndex;
    vector<int> astrometricExpnum;   // Maps from row number (=index) to expnum
    vector<int> nonAstrometricExpnum;
    // Build Exposure object for given exposure number, putting exposure-level
    // info into specified reference frame
    Exposure* buildExposure(int exposureNumber,
			    const Frame& frame,
			    const Ephemeris& ephem)  const;
  };
  
  class TransientTable {
    // Class holding a table of transient detections.
  public:
    // Null string input will look to an environment variable.
    TransientTable(const string transientFile="");
    bool hasExpnum(int expnum) const;  // Are there transients for this exposure?
    // Fill the exposure structure with info for all its transients, if not done already.
    // Return false if there are none.  Thread-safe.
    bool fillExposure(Exposure* eptr) const;
    // Return an Observation object for the transient in selected row
    Observation getObservation(int objectID,
			       const Ephemeris& ephem,
			       ExposureTable& exposures) const;
    // Return value of some named quantity
    template<class T>
    T getValue(const string& column, int row) const {
      T out;
      transientTable.readCell(out, column, row);
      return out;
    }
  private:
    img::FTable transientTable;
    img::FTable transientIndex;
    std::map<int,int> findExpnum;  // Lookup table from expnum to transientIndex row number
  };
  
  class CornerTable {
    // Class holding a table of CCD corner locations.
  public:
    // Null string input will look to an environment variable.
    CornerTable(const string cornerFile="");
    bool hasExpnum(int expnum) const;  // Are there corners for this exposure?
    // Fill the exposure structure with info for its CCDs.
    // Return false if there are none.
    bool fillExposure(Exposure* eptr) const;
  private:
    img::FTable cornerTable;
    std::multimap<int,int> findExpnum;  // Lookup table from expnum to cornerTable row numbers
    static std::map<string,int> detpos2ccdnum; // DECam table of DETPOS vs CCDNUM.
  };

  class Node {
    // Node of a k-d tree of Exposures over time and position
  public:
    typedef vector<Exposure*>::iterator dataIter;
    Node( dataIter _begin, dataIter _end,
	  const Ephemeris& ephem,
	  const Frame& frame);
    ~Node() {
      if (left) delete left;
      if (right) delete right;
    }
    
    enum Coordinate {TIME, X, Y};  // The three axes of the tree space
    
    list<Exposure*> find(const Fitter& path) const;
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

    // Don't need EIGEN_NEW since no static-sized Eigen members.
    
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
    DESTree(std::vector<Exposure*>& exposurePointers,
	    const Ephemeris& ephem,
	    const Frame& frame,
	    double gamma0);
    ~DESTree();
    
    int countNodes() const;
    list<Exposure*> find(const Fitter& path) const;
    // Return a list of exposures that are close to the path.
  private:
    std::list<Node*> years;
    static std::list<double> tdb_splits;  // tdb's at which we split DES seasons
  };

} // end namespace orbits
#endif 
