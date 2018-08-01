#include <algorithm>
#include "Exposures.h"
#include "FitsTable.h"
#include "Trajectory.h"

using namespace orbits;

// Name of environment variable giving path to observatories file
const string EXPOSURE_ENVIRON = "DES_EXPOSURE_TABLE";

ExposureTable::ExposureTable(string exposureFile) {
  string path = exposureFile;
  if (path.empty())  {
    // Find FITS file name from environment variable
    char *kpath = getenv(EXPOSURE_ENVIRON.c_str());
    if (kpath == NULL) 
      throw std::runtime_error("No path given for DES exposure table FITS file");
    path = kpath;
  }
  // ??? Could make this faster by saving columns as arrays, save FTable overhead.
  
  try {
    FITS::FitsTable ft1(path,FITS::ReadOnly,1);
    astrometricTable = ft1.extract();
    FITS::FitsTable ft2(path,FITS::ReadOnly,2);
    nonAstrometricTable = ft2.extract();
  } catch (FITS::FITSError& m) {
    cerr << "Error reading DES exposure table from " << path << endl;
    quit(m);
  }

  // Build indices into each table from maps
  vector<int> expnum;
  astrometricTable.readCells(expnum,"EXPNUM");
  for (int i=0; i<expnum.size(); i++)
    astrometricIndex[expnum[i]] = i;
  nonAstrometricTable.readCells(expnum,"EXPNUM");
  for (int i=0; i<expnum.size(); i++)
    nonAstrometricIndex[expnum[i]] = i;
}

bool
ExposureTable::isAstrometric(int expnum) const {
  return astrometricIndex.count(expnum);
}

double
ExposureTable::mjd(int expnum) const {
  auto ptr = astrometricIndex.find(expnum);
  int index;
  double mjd;
  if (ptr==astrometricIndex.end()) {
    // Try non-ast table; throw exception if not there.
    index = nonAstrometricIndex.at(expnum);
    nonAstrometricTable.readCell(mjd,"MJD_MID",index);
  } else {
    index = ptr->second;
    astrometricTable.readCell(mjd,"MJD_MID",index);
  }
  return mjd;
}

astrometry::CartesianICRS
ExposureTable::observatory(int expnum) const {
  auto ptr = astrometricIndex.find(expnum);
  int index;
  vector<double> v3;
  if (ptr==astrometricIndex.end()) {
    // Try non-ast table; throw exception if not there.
    index = nonAstrometricIndex.at(expnum);
    nonAstrometricTable.readCell(v3,"OBSERVATORY",index);
  } else {
    index = ptr->second;
    astrometricTable.readCell(v3,"OBSERVATORY",index);
  }
  return astrometry::CartesianICRS(v3[0],v3[1],v3[2]);
}

Matrix22
ExposureTable::atmosphereCov(int expnum) const {
  auto ptr = astrometricIndex.find(expnum);
  int index;
  Matrix22 out(0.);
  if (ptr==astrometricIndex.end()) {
    // Negative diagonals for unknown cov
    out(0,0) = out(1,1) = -1.;
  } else {
    vector<double> cov;
    index = ptr->second;
    astrometricTable.readCell(cov,"COV",index);
    const double masSq = pow(0.001*ARCSEC,2.); // Table is in milliarcsec
    out(0,0) = cov[0]*masSq;
    out(1,1) = cov[1]*masSq;
    out(0,1) = cov[2]*masSq;
    out(1,0) = cov[2]*masSq;
  }
  return out;
}

// Get both at once, return false if expnum is not known:
bool
ExposureTable::observingInfo(int expnum,
			     double& mjd,
			     astrometry::CartesianICRS& observatory) const {
  auto ptr = astrometricIndex.find(expnum);
  int index;
  vector<double> v3;
  if (ptr==astrometricIndex.end()) {
    // Try non-ast table
    ptr = nonAstrometricIndex.find(expnum);
    if (ptr==nonAstrometricIndex.end())
      // No data
      return false;
    index = ptr->second;
    nonAstrometricTable.readCell(mjd,"MJD_MID",index);
    nonAstrometricTable.readCell(v3,"OBSERVATORY",index);
  } else {
    index = ptr->second;
    astrometricTable.readCell(mjd,"MJD_MID",index);
    astrometricTable.readCell(v3,"OBSERVATORY",index);
  }
  for (int i=0; i<3; i++) observatory[i] = v3[i];
  return true;
}

std::vector<const Exposure*>
orbits::selectExposures(const Frame& frame,   // Starting coordinates, time
			const Ephemeris& ephem,  
			double gamma0,        // Center and width of range 
			double dGamma,        // of gamma to cover
			double searchRadius,  // Range of starting coords to cover
			string transientFile,
			string exposureFile,
			double fieldRadius) { 


  std::vector<const Exposure*> out;  // This will be the returned array
  std::list<Exposure*> firstCut; // Make a list first, then toss those with no transient data.
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
    // FITS table has RA, Dec in degrees.
    double ra, dec;
    exposureTable.readCell(ra, "ra", iexp);
    exposureTable.readCell(dec, "dec", iexp);
    astrometry::SphericalICRS pointing(ra*DEGREE, dec*DEGREE);

    if (pointing.distance(frame.orient.getPole()) > 60*DEGREE)  // 60 degree cut!!
      continue;
    
    // Fill Exposure structure, converting to Frame
    auto expo = new Exposure;
    {
      double mjd;
      exposureTable.readCell(mjd, "mjd_mid", iexp);
      expo->mjd = mjd;
      expo->tdb = ephem.mjd2tdb(mjd);
      expo->tobs = ephem.mjd2tdb(mjd) - frame.tdb0;
    }

    {
      int expnum;
      exposureTable.readCell(expnum, "expnum", iexp);
      expo->expnum = expnum;
    }
    
    {
      astrometry::Gnomonic xy(pointing, frame.orient);
      double x,y;
      xy.getLonLat(x,y);
      expo->axis[0] = x;
      expo->axis[1] = y;
      std::vector<double> cov;
      // Read atmospheric covariance and rotate into Frame
      exposureTable.readCell(cov, "cov", iexp);
      Matrix22 covRADec;
      covRADec(0,0) = cov[0];
      covRADec(1,1) = cov[1];
      covRADec(1,0) = cov[2];
      covRADec(0,1) = cov[2];
      covRADec *= pow(0.001*ARCSEC,2.);  // milliarcsec->rad
      expo->atmosphereCov = frame.fromICRS(covRADec);
    }

    astrometry::CartesianICRS earthICRS;
    {
      std::vector<double> observatory;
      exposureTable.readCell(observatory,"observatory",iexp);
      for (int i=0; i<3; i++) earthICRS[i]=observatory[i];
      expo->earth = frame.fromICRS(earthICRS.getVector());
    }
    
    // Get the center positions of the orbits at this time
    double xMin, xMax, yMin, yMax;
    astrometry::SphericalICRS posn =
      trajectoryMin.observe(expo->tdb,earthICRS);
    astrometry::Gnomonic(posn, frame.orient).getLonLat(xMin,yMin);
    posn = trajectoryMax.observe(expo->tdb, earthICRS);
    astrometry::Gnomonic(posn, frame.orient).getLonLat(xMax,yMax);

    Point mean;
    mean[0] = 0.5*(xMin+xMax);
    mean[1] = 0.5*(yMin+yMax);
    double gammaRadius = 0.5*hypot(xMax-xMin, yMax-yMin);

    // Build matching radius, starting with
    // the radii of the exposure and the source region
    double radius = fieldRadius + searchRadius;
    // Add the allowance for varying gamma:
    radius += gammaRadius;

    // Calculate unbinding velocity
    const double BINDING_FACTOR = 1.1;
    double bind = BINDING_FACTOR * sqrt(8.) * PI * pow(gammaMax,1.5) * (expo->tdb-frame.tdb0);
    radius += asin(bind);

    // Now test distance
    if ( mean.distanceSq(expo->axis) < radius*radius)
      firstCut.push_back(expo);
    else
      delete expo;
  }

  // Now read in the transients for the exposures we're keeping.
  img::FTable transientTable;
  img::FTable transientIndex;
  try {
    FITS::FitsTable ft(transientFile, FITS::ReadOnly, "DATA");
    std::vector<string> columnsOfInterest = {"RA","DEC","CCDNUM","ERRAWIN_WORLD"};
    transientTable = ft.extract(0,-1,columnsOfInterest);
    FITS::FitsTable ft2(transientFile, FITS::ReadOnly, "INDEX");
    transientIndex = ft2.extract();
  } catch (FITS::FITSError& m) {
    cerr << "Error reading transient table" << endl;
    quit(m);
  }

  // Build an index into the index:
  std::map<int,int> findExpnum;
  int expnum;
  for (int row=0; row<transientIndex.nrows(); row++) {
    transientIndex.readCell(expnum, "EXPNUM", row);
    findExpnum[expnum] = row;
  }

  for (auto expoptr : firstCut) {
    if (findExpnum.count(expoptr->expnum)==0) {
      // No transient file entry for this exposure.  Flush it.
      delete expoptr;
      continue;
    }

    int begin,end;
    transientIndex.readCell(begin, "FIRST", findExpnum[expoptr->expnum]);
    transientIndex.readCell(end,  "LAST", findExpnum[expoptr->expnum]);
    float density;
    transientIndex.readCell(density, "DENSITY", findExpnum[expoptr->expnum]);
    expoptr->detectionDensity = density / (DEGREE*DEGREE); // Change per sq deg to per sr
    int nTransients = end-begin;
    expoptr->xy.resize(nTransients, 2);
    expoptr->covXX.resize(nTransients);
    expoptr->covXY.resize(nTransients);
    expoptr->covYY.resize(nTransients);
    expoptr->ccdnum.resize(nTransients);
    expoptr->id.resize(nTransients);
    expoptr->valid = BVector(nTransients,true); //Everyone is valid to start with.

    std::vector<double> ra;
    transientTable.readCells(ra, "RA", begin, end);
    std::vector<double> dec;
    transientTable.readCells(dec, "DEC", begin, end);
    std::vector<double> sig;
    transientTable.readCells(sig, "ERRAWIN_WORLD", begin, end);
    std::vector<short int> ccd;
    transientTable.readCells(ccd, "CCDNUM", begin, end);
    
    // Fill in individual detections' properties
    // Their coordinates and sigma are in degrees.
    DMatrix xyRADec(nTransients,2,0.);
    for (int i=0; i<nTransients; i++) {
      expoptr->id[i] = i+begin;
      expoptr->ccdnum[i] = ccd[i];
      expoptr->covXX[i] = sig[i]*sig[i]*DEGREE*DEGREE + expoptr->atmosphereCov(0,0);
      expoptr->covXY[i] = expoptr->atmosphereCov(1,0);
      expoptr->covYY[i] = sig[i]*sig[i]*DEGREE*DEGREE + expoptr->atmosphereCov(1,1);
      astrometry::SphericalICRS radec(ra[i]*DEGREE, dec[i]*DEGREE);
      astrometry::Gnomonic xy(radec, frame.orient);
      double x,y;
      xy.getLonLat(x,y);
      expoptr->xy(i,0) = x;
      expoptr->xy(i,1) = y;
    }

    // Add to output
    out.push_back(expoptr);
  } // End exposure loop

  return out;
}

////////////////////////////////////////////////////////////////////////
/* kD-tree nodes */
////////////////////////////////////////////////////////////////////////

// Tree-wide static variables:
double
Node::speed = TPI / 40.;
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
				earth(3,3,0.),
				corners(4,2,0.)
{
  // Determine bounds of member exposures.
  // Rows 0 and 2 of corners are lower and upper bounds
  corners.row(0) = (*begin)->axis.transpose();
  corners.row(2) = (*begin)->axis.transpose();
  tStart = (*begin)->tobs;
  auto last = end;
  --last;
  tStop = (*last)->tobs;  // Time ordering assumed.
  for (auto i=begin; i!=end; ++i) {
    corners(0,0) = min(corners(0,0), (*i)->axis[0]);
    corners(2,0) = max(corners(2,0), (*i)->axis[0]);
    corners(0,1) = min(corners(0,1), (*i)->axis[1]);
    corners(2,1) = max(corners(2,1), (*i)->axis[1]);
  }
  // Make the other two corners:
  corners(1,0) = corners(0,0);
  corners(1,1) = corners(2,1);
  corners(3,0) = corners(2,0);
  corners(3,1) = corners(0,1);
  
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
    return exptr->tobs < value;
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
  double dx = corners(2,0)-corners(0,0);
  double dy = corners(2,1)-corners(0,1);
  double dt = (tStop-tStart)*speed;
  double splittingThreshold = 0.5 * fieldRadius;
  if (dx < splittingThreshold
      && dy < splittingThreshold
      && dt < splittingThreshold) {
    /**cerr << "Leaf node: " << end-begin << " elements" 
	     << " T " << tStart << "-" << tStop
	     << " X " << corners(0,0) << "-" << corners(2,0)
	     << " Y " << corners(0,1) << "-" << corners(2,1)
	     << endl; **/
    return; // No need for further splits; leaf node.
  }
  // Now split, including stable partition of parent array such
  // that leaf nodes remain in time order.
  dataIter mid;
  if (dt > dx && dt > dy) {
    splitOn = TIME;
    splitValue = 0.5*(tStart + tStop);
    TimeSplit splitter(splitValue);
    mid = stable_partition(begin, end, splitter);
    //**/cerr << ".Splitting at time " << splitValue << endl;
  } else if (dx > dy) {
    splitOn = X;
    splitValue = 0.5*(corners(0,0)+corners(2,0));
    XSplit splitter(splitValue);
    mid = stable_partition(begin, end, splitter);
    //**/cerr << ".Splitting at X " << splitValue << endl;
  } else {
    splitOn = Y;
    splitValue = 0.5*(corners(0,1)+corners(2,1));
    YSplit splitter(splitValue);
    mid = stable_partition(begin, end, splitter);
    //**/cerr << ".Splitting at Y " << splitValue << endl;
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
  Point ctr;
  Vector2 motion;
  // Model motion as a quadratic, [xy] = a + bt + ct^2,
  // where t goes from -1 to +1 over time interval
  double a = x[1];
  double b = 0.5 * (x[2]-x[0]);
  double c = 0.5 * (x[0] + x[2] - 2*a);
  if (c ==0. || abs(b / (2*c))>1.) {
    // Endpoints are the maximum
    ctr[0] = 0.5*(x[2]+x[0]);
    motion[0] = 0.5*abs(x[2]-x[0]);
  } else if (c<0) {
    // There is a max during interval
    double xMax = a - b*b/(4*c);
    double xMin = b>0 ? x[0] : x[2];
    ctr[0] = 0.5*(xMax+xMin);
    motion[0] = 0.5*(xMax-xMin);
  } else {
    // There is a minimum during interval
    double xMin = a - b*b/(4*c);
    double xMax = b>0 ? x[2] : x[0];
    ctr[0] = 0.5*(xMax+xMin);
    motion[0] = 0.5*(xMax-xMin);
  }
  // Y motion:
  a = y[1];
  b = 0.5 * (y[2]-y[0]);
  c = 0.5 * (y[0] + y[2] - 2*a);
  if (c ==0. || abs(b / (2*c))>1.) {
    // Endpoints are the maximum
    ctr[1] = 0.5*(y[2]+y[0]);
    motion[1] = 0.5*abs(y[2]-y[0]);
  } else if (c<0) {
    // There is a max during interval
    double xMax = a - b*b/(4*c);
    double xMin = b>0 ? y[0] : y[2];
    ctr[1] = 0.5*(xMax+xMin);
    motion[1] = 0.5*(xMax-xMin);
  } else {
    // There is a minimum during interval
    double xMin = a - b*b/(4*c);
    double xMax = b>0 ? y[2] : y[0];
    ctr[1] = 0.5*(xMax+xMin);
    motion[1] = 0.5*(xMax-xMin);
  }
  // Find maximum of major axes
  DVector tr2 = 0.5*(covxx + covyy);
  DVector det = covxx.array()*covyy.array() - covxy.array()*covxy.array();
  double maxasq = (tr2.array() + sqrt(tr2.array()*tr2.array() - det.array())).maxCoeff();

  // Add motion, error, and field radii together to get match radius
  double matchRadius = hypot(motion[0],motion[1]) + sqrt(maxasq) + fieldRadius;

  list<const Exposure*> out;
  // If no spatial intersection: (ignoring curved boundaries at corners)
  if (ctr[0] >= corners(2,0) + matchRadius ||
      ctr[0] <= corners(0,0) - matchRadius ||
      ctr[1] >= corners(2,1) + matchRadius ||
      ctr[1] <= corners(0,1) - matchRadius) {
    // No intersection, return empty list
    return out;
  }

  if (!left) {
    // If this is a leaf node, return everything as a possibility
    out.insert(out.end(), begin, end);
    //**/cerr << "Leaf node in with " << end-begin << " exposures" << endl;
    return out;
  }

  // If Node is entirely contained in search - check all 4 corners
  if ( (ctr.distanceSq(corners).array() < (matchRadius*matchRadius)).all()) {
    // Return all exposures
    out.insert(out.end(), begin, end);
    //**/cerr << "Total node in with " << end-begin << " exposures" << endl;
    return out;
  }
    
  // If Node is only partially in search region
  out = left->find(path);
  out.splice(out.end(),right->find(path));

  return out;
}


Node*
Node::buildTree(dataIter begin_, dataIter end_,
		const Ephemeris& ephem,
		const Frame& frame) {
  auto root = new Node(begin_, end_, ephem, frame);
  root->split(ephem,frame);
  return root;
}

std::list<double>
DESTree::tdb_splits = {13.5,14.5,15.5,16.5,17.5,18.5};  // July 1 of each year splits DES seasons

DESTree::DESTree(std::vector<const Exposure*>& exposurePointers,
		 const Ephemeris& ephem,
		 const Frame& frame,
		 double gamma0) {
  // Set up the Node class
  Node::setSpeed(TPI * gamma0); // Speed is max reflex at gamma0, rad/yr
  Node::setFieldRadius(1.1*DEGREE);    // DECam radius
  Node::setObservatory(807);    // CTIO
  
  auto begin = exposurePointers.begin();
  auto end = exposurePointers.end();

  // Make a distinct tree for each DES season
  for (auto tdb : tdb_splits) {
    TimeSplit splitter(tdb-frame.tdb0);
    auto mid = stable_partition(begin, end, splitter);
    if (mid-begin > 0) {
      // We have a tree to build for this year.
      years.push_back( Node::buildTree(begin, mid, ephem, frame));
      begin = mid;
    }
  }
  // Add any left after last split time
  if (end-begin > 0) 
    years.push_back( Node::buildTree(begin, end, ephem, frame));
  /**/cerr << "DESTree has " << years.size() << " seasons" << endl;
}

DESTree::~DESTree() {
  // Delete all years' trees
  for (auto n : years) delete n;
}

int
DESTree::countNodes() const {
  int count=0;
  for (auto n : years)
    count += n->countNodes();
  return count;
}

list<const Exposure*>
DESTree::find(const Fitter& path) const {
  list<const Exposure*> out;
  for (auto n : years) 
    out.splice(out.end(),n->find(path));
  return out;
}

DVector
Exposure::chisq(double x0, double y0, double covxx, double covyy, double covxy) const {
  DVector dx = xy.col(0).array() - x0;
  DVector dy = xy.col(1).array() - y0;
  DVector cxx = covXX.array() + covxx;
  DVector cxy = covXY.array() + covxy;
  DVector cyy = covYY.array() + covyy;

  DVector det = cxx.array() * cyy.array() - cxy.array()*cxy.array();
  DVector chisq = dx.array() * (dx.array()*cyy.array() - dy.array()*cxy.array())
    + dy.array() * (dy.array()*cxx.array() - dx.array()*cxy.array());
  chisq.array() /= det.array();
  return chisq;
}
  
