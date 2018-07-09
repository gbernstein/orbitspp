// Manipulate information on DECam exposures
#ifndef EXPOSURES_H
#define EXPOSURES_H

#include "Astrometry.h"
#include "OrbitTypes.h"
#include "PlaneGeometry.h"
#include "Ephemeris.h"
#include "Fitter.h"

namespace orbits {

  class Exposure {
  public:
    int expnum;
    double mjd;
    double tdb;
    astrometry::Vector3 earth;  // Observatory position in Frame
    Point axis;   // Optical axis coords, in Frame gnomonic, degrees

    double detectionDensity; // Density of transients, per sq degree

    // ??? Add per-CCD info on boundaries, detection lists.
  };

  class Node {
    // Node of a k-d tree of Exposures over time and position
  public:
    typedef vector<const Exposure*>::iterator dataIter;
    Node( dataIter _begin, dataIter _end);
    ~Node() {
      if (left) delete left;
      if (right) delete right;
    }
    
    enum Coordinate {TIME, X, Y};  // The three axes of the tree space
    
    Point lower;  // Corners of the node
    Point upper;
    double tStart;  // Time (relative to reference frame tdb0) of
    double tStop;   // first and last exposure in set

    void split();
    // Split all the way down to leaves, starting with 

    list<const Exposure*> find(const Fitter& path) const;
    // Return a list of exposures that are close to the path.

    static void setSpeed(double speed_) {speed=speed_;}
    static void setFieldRadius(double radius_) {fieldRadius = radius_;}
    
  private:
    dataIter begin;  // Range of data associated with this node
    dataIter end;
    Node* left;	// Children (nullptrs if leaf node)
    Node* right;
    Coordinate splitOn;  // Which coordinate does the split happen on?
    double splitValue;   // Value for split

    static double speed;  // Typical object speed, used to choose time vs distance split
    static double fieldRadius; // Size of exposure FOV.  Also sets the scale for radius.
  };
  
  extern
  std::vector<Exposure>
  selectExposures(const Frame& frame,   // Starting coordinates, time
		  const Ephemeris& ephem,  // Ephemeris to use
		  double gamma0,        // Center and width of range 
		  double dGamma,        // of gamma to cover
		  double searchRadius,  // Range of starting coords to cover
		  string exposureFile="data/y4a1.exposure.positions.fits",  // File of exposure data
		  double fieldRadius = 1.1); // Circumscribed field radius (degrees)
  // ?? Allow exclusion of filters??
  // ?? Add CCD corners, detections ??
  

} // end namespace orbits
#endif 
