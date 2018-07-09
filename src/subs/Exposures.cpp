#include <algorithm>
#include "Exposures.h"
#include "FitsTable.h"
#include "Trajectory.h"

using namespace orbits;

std::vector<Exposure>
orbits::selectExposures(const Frame& frame,   // Starting coordinates, time
			const Ephemeris& ephem,  
			double gamma0,        // Center and width of range 
			double dGamma,        // of gamma to cover
			double searchRadius,  // Range of starting coords to cover
			string exposureFile,
			double fieldRadius) { 


  std::vector<Exposure> out;  // This will be the returned array

  // Read the exposure information table
  img::FTable exposureTable;
  try {
    FITS::FitsTable ft(exposureFile);
    exposureTable = ft.extract();
  } catch (FITS::FITSError& m) {
    cerr << "Error reading exposure table" << endl;
    quit(m);
  }

  // Create trajectories for initially stationary targets to
  // calculate parallax + gravity displacements
  double gammaMin = gamma0 - 0.5*dGamma;
  double gammaMax = gamma0 + 0.5*dGamma;

  // Make orbits that start from standstill
  ABG abg;
  abg[ABG::G]=gammaMin;

  astrometry::Vector3 x0,v0;
  abg.getState(0.,x0, v0);
  State s0;
  s0.tdb = frame.tdb0;
  s0.x = frame.toICRS(x0);
  s0.v = frame.toICRS(v0,true);  // Should be zero but do this anyway
  Trajectory trajectoryMin(ephem, s0, Gravity::BARY);

  abg[ABG::G]=gammaMax;
  abg.getState(0.,x0, v0);
  s0.x = frame.toICRS(x0);
  s0.v = frame.toICRS(v0,true);  // Should be zero but do this anyway
  Trajectory trajectoryMax(ephem, s0, Gravity::BARY);

  // Step through all exposures
  for (int iexp = 0; iexp<exposureTable.nrows(); iexp++) {
    // Skip exposures not pointed anywhere near right direction
    double ra, dec;
    exposureTable.readCell(ra, "ra", iexp);
    exposureTable.readCell(dec, "dec", iexp);
    astrometry::SphericalICRS pointing(ra*DEGREE, dec*DEGREE);

    if (pointing.distance(frame.orient.getPole()) > 60*DEGREE)  // 60 degree cut!!
      continue;
    
    // Fill Exposure structure, converting to Frame
    Exposure expo;
    {
      double mjd;
      exposureTable.readCell(mjd, "mjd_mid", iexp);
      expo.mjd = mjd;
      expo.tdb = ephem.mjd2tdb(mjd);
    }

    {
      int expnum;
      exposureTable.readCell(expnum, "expnum", iexp);
      expo.expnum = expnum;
    }
    
    {
      astrometry::Gnomonic xy(pointing, frame.orient);
      double x,y;
      xy.getLonLat(x,y);
      expo.axis[0] = x / DEGREE;
      expo.axis[1] = y / DEGREE;
    }

    astrometry::CartesianICRS earthICRS;
    {
      std::vector<double> observatory;
      exposureTable.readCell(observatory,"observatory",iexp);
      for (int i=0; i<3; i++) earthICRS[i]=observatory[i];
      expo.earth = frame.fromICRS(earthICRS.getVector());
    }
    
    // Get the center positions of the orbits at this time
    double xMin, xMax, yMin, yMax;
    astrometry::SphericalICRS posn =
      trajectoryMin.observe(expo.tdb,earthICRS);
    astrometry::Gnomonic(posn, frame.orient).getLonLat(xMin,yMin);
    posn = trajectoryMax.observe(expo.tdb, earthICRS);
    astrometry::Gnomonic(posn, frame.orient).getLonLat(xMax,yMax);

    Point mean;
    mean[0] = 0.5*(xMin+xMax) / DEGREE;
    mean[1] = 0.5*(yMin+yMax) / DEGREE;
    double gammaRadius = 0.5*hypot(xMax-xMin, yMax-yMin) / DEGREE;

    // Build matching radius, starting with
    // the radii of the exposure and the source region
    double radius = fieldRadius + searchRadius;
    // Add the allowance for varying gamma:
    radius += gammaRadius;

    // Calculate unbinding velocity
    const double BINDING_FACTOR = 1.1;
    double bind = BINDING_FACTOR * sqrt(8.) * PI * pow(gammaMax,1.5) * (expo.tdb-frame.tdb0);
    radius += asin(bind)/DEGREE;

    // Now test distance
    if ( mean.distanceSq(expo.axis) < radius*radius)
      out.push_back(expo);
  }

  return out;
}

////////////////////////////////////////////////////////////////////////
/* kD-tree nodes */
////////////////////////////////////////////////////////////////////////

// Tree-wide static variables:
double
Node::speed;

double
Node::fieldRadius;

Node::Node(dataIter _begin, dataIter _end): begin(_begin), end(_end),
					    left(nullptr), right(nullptr),
					    splitOn(TIME)
{
  // Determine bounds of member exposures
  lower = (*begin)->axis;
  upper = (*begin)->axis;
  tStart = (*begin)->tdb; // ??? reference to tdb0??
  auto last = end;
  --last;
  tStop = (*end)->tdb; // ??? assumes buffer is time-ordered
  for (auto i=begin; i!=end; ++i) {
    lower[0] = min(lower[0], (*i)->axis[0]);
    upper[0] = max(upper[0], (*i)->axis[0]);
    lower[1] = min(lower[1], (*i)->axis[1]);
    upper[1] = max(upper[1], (*i)->axis[1]);
  }
}

// Unary predicates for partitioning the exposures
class TimeSplit {
public:
  TimeSplit(double value_): value(value_) {}
  bool operator()(const Exposure* exptr) const {
    return exptr->tdb < value;
  }
private:
  double value;
};
class XSplit {
public:
  XSplit(double value_): value(value_) {}
  bool operator()(const Exposure* exptr) const {
    return exptr->axis[0] < value;
  }
private:
  double value;
};
class YSplit {
public:
  YSplit(double value_): value(value_) {}
  bool operator()(const Exposure* exptr) const {
    return exptr->axis[1] < value;
  }
private:
  double value;
};

void
Node::split() {
  // Find largest dimension
  double dx = upper[0] - lower[0];
  double dy = upper[1] - lower[1];
  double dt = (tStop-tStart)*speed;
  double splittingThreshold = 0.5 * fieldRadius;
  if (dx < splittingThreshold
      && dy < splittingThreshold
      && dt < splittingThreshold)
    return; // No need for further splits; leaf node.
  
  // Now split, including stable partition of parent array such
  // that leaf nodes remain in time order.
  if (dt > dx && dt > dy) {
    splitOn = TIME;
    splitValue = 0.5*(tStart + tStop);
    TimeSplit splitter(splitValue);
    auto mid = stable_partition(begin, end, splitter);
    left = new Node(begin,mid);
    right = new Node(mid, end);
  } else if (dx > dy) {
    splitOn = X;
    splitValue = 0.5*(lower[0]+upper[0]);
    XSplit splitter(splitValue);
    auto mid = stable_partition(begin, end, splitter);
    left = new Node(begin,mid);
    right = new Node(mid, end);
  } else {
    splitOn = Y;
    splitValue = 0.5*(lower[1]+upper[1]);
    YSplit splitter(splitValue);
    auto mid = stable_partition(begin, end, splitter);
    left = new Node(begin,mid);
    right = new Node(mid, end);
  }
  // Recursive split of children
  left->split();
  right->split();
}
