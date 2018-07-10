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
      expo.tobs = ephem.mjd2tdb(mjd) - frame.tdb0;
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
Node::speed = 360. / 40.;
// Default apparent motion to reflex at 40 AU, in degrees/yr.

double
Node::fieldRadius = 1.1;  // Default to DECam radius, in degrees

int
Node::obscode = 807;  // Observatory code for earth positions.

Node::Node(dataIter _begin, dataIter _end,
	   const Ephemeris& ephem,
	   const Frame& frame): begin(_begin), end(_end),
				left(nullptr), right(nullptr),
				splitOn(TIME),
				tobs(3,0.),
				earth(3,3,0.)
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
  // Set up the start/middle/end observing times and locations
  tobs[0] = tStart;
  tobs[1] = 0.5*(tStart+tStop);
  tobs[2] = tStop;
  earth.row(0) = (*begin)->earth.transpose();
  earth.row(2) = (*last)->earth.transpose();
  // Use ephemeris to get observer position at midpoint, since there
  // is not going to be an exposure at this time.
  // Get observatory position
  astrometry::CartesianICRS observatory = ephem.observatory(obscode,
							    tobs[1]+frame.tdb0);
  earth.row(1) = frame.fromICRS(observatory.getVector()).transpose();
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
Node::split(const Ephemeris& ephem,
	    const Frame& frame) {
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
  dataIter mid;
  if (dt > dx && dt > dy) {
    splitOn = TIME;
    splitValue = 0.5*(tStart + tStop);
    TimeSplit splitter(splitValue);
    mid = stable_partition(begin, end, splitter);
  } else if (dx > dy) {
    splitOn = X;
    splitValue = 0.5*(lower[0]+upper[0]);
    XSplit splitter(splitValue);
    mid = stable_partition(begin, end, splitter);
  } else {
    splitOn = Y;
    splitValue = 0.5*(lower[1]+upper[1]);
    YSplit splitter(splitValue);
    mid = stable_partition(begin, end, splitter);
  }
  // Recursive split of children
  left = new Node(begin,mid,ephem,frame);
  right = new Node(mid, end,ephem,frame);
  left->split(ephem,frame);
  right->split(ephem,frame);
}

list<const Exposure*>
Node::find(const Fitter& path) const {
  // Get locations and error ellipses for orbit and temporal endpoints, midpoint.
  DVector x(3),y(3),covxx(3),covyy(3),covxy(3);
  path.predict(tobs, earth, &x, &y, &covxx, &covyy, &covxy);

  // Turn these into a center and radius for search.

  // Add field radius to search

  list<const Exposure*> out;
  // If no spatial intersection:
  return out;

  // If this is a leaf node, return everything as a possibility

  // If Node is entirely contained in search, return everything

  // If partial coverage, descend tree
  out = left->find(path);
  out.splice(out.end(),right->find(path));

  return out;
}

  
