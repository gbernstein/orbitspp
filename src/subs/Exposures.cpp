#include <algorithm>
#include "Exposures.h"
#include "FitsTable.h"
#include "Trajectory.h"
#include "StringStuff.h"

using namespace orbits;

// Name of environment variable giving path to observatories file
const string EXPOSURE_ENVIRON = "DES_EXPOSURE_TABLE";
const string TRANSIENT_ENVIRON = "DES_TRANSIENT_TABLE";
const string CORNERS_ENVIRON = "DES_CORNER_TABLE";

// Ignore any transient whose ERRAWIN_WORLD is above this:
const double MAXIMUM_TRANSIENT_POSITION_ERROR = 1.*ARCSEC;

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
  astrometricTable.readCells(astrometricExpnum,"EXPNUM");
  for (int i=0; i<astrometricExpnum.size(); i++)
    astrometricIndex[astrometricExpnum[i]] = i;
  nonAstrometricTable.readCells(nonAstrometricExpnum,"EXPNUM");
  for (int i=0; i<nonAstrometricExpnum.size(); i++)
    nonAstrometricIndex[nonAstrometricExpnum[i]] = i;
}

bool
ExposureTable::isAstrometric(int expnum) const {
  return astrometricIndex.count(expnum);
}

int
ExposureTable::expnum(int index) const {
  if (index<0) {
    throw std::runtime_error("ExposureTable index out of range");
  } else if (index < astrometricExpnum.size()) {
    return astrometricExpnum[index];
  }
  index -= astrometricExpnum.size();
  if (index < nonAstrometricExpnum.size()) {
    return nonAstrometricExpnum[index];
  } else {
    throw std::runtime_error("ExposureTable index out of range");
  }
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

string
ExposureTable::band(int expnum) const {
  auto ptr = astrometricIndex.find(expnum);
  int index;
  string band;
  if (ptr==astrometricIndex.end()) {
    // Try non-ast table; throw exception if not there.
    index = nonAstrometricIndex.at(expnum);
    nonAstrometricTable.readCell(band,"FILTER",index);
  } else {
    index = ptr->second;
    astrometricTable.readCell(band,"FILTER",index);
  }
  return band;
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

astrometry::SphericalICRS
ExposureTable::axis(int expnum) const {
  auto ptr = astrometricIndex.find(expnum);
  int index;
  double ra, dec;
  if (ptr==astrometricIndex.end()) {
    // Try non-ast table; throw exception if not there.
    index = nonAstrometricIndex.at(expnum);
    nonAstrometricTable.readCell(ra,"RA",index);
    nonAstrometricTable.readCell(dec,"DEC",index);
  } else {
    index = ptr->second;
    astrometricTable.readCell(ra,"RA",index);
    astrometricTable.readCell(dec,"DEC",index);
  }
  return astrometry::SphericalICRS(ra*DEGREE, dec*DEGREE);
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
			     astrometry::CartesianICRS& observatory,
			     astrometry::SphericalICRS& axis) const {
  auto ptr = astrometricIndex.find(expnum);
  int index;
  vector<double> v3;
  double ra, dec;
  if (ptr==astrometricIndex.end()) {
    // Try non-ast table
    ptr = nonAstrometricIndex.find(expnum);
    if (ptr==nonAstrometricIndex.end())
      // No data
      return false;
    index = ptr->second;
    nonAstrometricTable.readCell(mjd,"MJD_MID",index);
    nonAstrometricTable.readCell(v3,"OBSERVATORY",index);
    nonAstrometricTable.readCell(ra,"RA",index);
    nonAstrometricTable.readCell(dec,"DEC",index);
  } else {
    index = ptr->second;
    astrometricTable.readCell(mjd,"MJD_MID",index);
    astrometricTable.readCell(v3,"OBSERVATORY",index);
    astrometricTable.readCell(ra,"RA",index);
    astrometricTable.readCell(dec,"DEC",index);
  }
  for (int i=0; i<3; i++) observatory[i] = v3[i];
  axis = astrometry::SphericalICRS(ra*DEGREE, dec*DEGREE);
  return true;
}

std::vector<Exposure*>
ExposureTable::getPool(const Frame& frame,   // Starting coordinates, time
		       const Ephemeris& ephem,  
		       double gamma0,        // Center and width of range 
		       double dGamma,        // of gamma to cover
		       double searchRadius,  // Range of starting coords to cover
		       bool astrometricOnly,
		       double fieldRadius) const {

  std::vector<Exposure*> out;  // This will be the returned array

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
  for (int iexp = 0; iexp<size(); iexp++) {
    int num = expnum(iexp);
    // Skip if this is not astrometric and we don't want it.
    if (astrometricOnly && !isAstrometric(num))
      continue;
    // Skip exposures not pointed anywhere near right direction
    if (axis(num).distance(frame.orient.getPole()) > 60*DEGREE)  // 60 degree cut!!
      continue;
    
    // Read all the exposure's info:
    auto expo = buildExposure(num, frame, ephem);
    
    // Get the center positions of the orbits at this time
    double xMin, xMax, yMin, yMax;
    astrometry::SphericalICRS posn =
      trajectoryMin.observe(expo->tdb,expo->earthICRS);
    astrometry::Gnomonic(posn, frame.orient).getLonLat(xMin,yMin);
    posn = trajectoryMax.observe(expo->tdb, expo->earthICRS);
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
      out.push_back(expo);
    else
      delete expo;
  }

  return out;
}

std::vector<Exposure*>
ExposureTable::getPool(const Frame& frame,   // Starting coordinates, time
		       const Ephemeris& ephem,  
		       double searchRadius,  // Range of coords to cover
		       bool astrometricOnly) const {

  std::vector<Exposure*> out;  // This will be the returned array

  // Step through all exposures
  for (int iexp = 0; iexp<size(); iexp++) {
    int num = expnum(iexp);
    // Skip if this is not astrometric and we don't want it.
    if (astrometricOnly && !isAstrometric(num))
      continue;
    // Skip exposures not pointed anywhere near right direction
    if (axis(num).distance(frame.orient.getPole()) > searchRadius)  
      continue;
    
    // Read all the exposure's info:
    out.push_back(buildExposure(num, frame, ephem));
  }
  return out;
}

Exposure*
ExposureTable::buildExposure(int exposureNumber,
			     const Frame& frame,
			     const Ephemeris& ephem) const {
    astrometry::SphericalICRS axisICRS;
    astrometry::CartesianICRS earthICRS;
    double mjd;
    observingInfo(exposureNumber, mjd, earthICRS, axisICRS);
    
    // Fill Exposure structure, converting to Frame
    auto expo = new Exposure;
    expo->expnum = exposureNumber;
    expo->mjd = mjd;
    expo->tdb = ephem.mjd2tdb(mjd);
    expo->tobs = ephem.mjd2tdb(mjd) - frame.tdb0;
    expo->earthICRS = earthICRS;
    expo->earth = frame.fromICRS(earthICRS.getVector());
    expo->astrometric = isAstrometric(exposureNumber);
    expo->localICRS.setPole(axisICRS);
    expo->localICRS.alignToICRS();
    
    {
      // Put exposure axis into Frame projection
      astrometry::Gnomonic xy(axisICRS, frame.orient);
      double x,y;
      xy.getLonLat(x,y);
      expo->axis[0] = x;
      expo->axis[1] = y;
    }
    {
      // Read atmospheric covariance and rotate into Frame
      Matrix22 covRADec = atmosphereCov(exposureNumber);
      expo->atmosphereCov = frame.fromICRS(covRADec);
    }
    return expo;
}
  
int
Exposure::whichCCD(const astrometry::SphericalICRS& radec) const {
  // Convert radec to projection about axis
  astrometry::Gnomonic gn(radec, localICRS);
  double xx,yy;
  gn.getLonLat(xx,yy);
  Point p(xx,yy);
  for (int i=0; i<deviceBoundsLocalICRS.size(); i++)
    if (deviceBoundsLocalICRS[i].inside(p))
      return devices[i];
  return 0;
}

vector<int>
Exposure::whichCCDs(const astrometry::SphericalCoords& radec,
		    const Matrix22& cov) const {
  // Convert input coord to projection about axis,
  // obtaining information about local derivatives
  astrometry::Gnomonic gn(localICRS);
  Matrix22 dLocaldInput;
  gn.convertFrom(radec, dLocaldInput);
  // New covariance matrix:
  Matrix22 localCov = dLocaldInput * cov * dLocaldInput.transpose();
  double xx,yy;
  gn.getLonLat(xx,yy);
  Point p(xx,yy);
  Ellipse e(p,localCov);
  vector<int> out;
  for (int i=0; i<deviceBoundsLocalICRS.size(); i++) 
    if (e.intersects(deviceBoundsLocalICRS[i])) 
      out.push_back(devices[i]);
  return out;
}

TransientTable::TransientTable(const string transientFile) {
  string path = transientFile;
  if (path.empty())  {
    // Find FITS file name from environment variable
    char *kpath = getenv(TRANSIENT_ENVIRON.c_str());
    if (kpath == NULL) 
      throw std::runtime_error("No path given for DES transient table FITS file");
    path = kpath;
  }
  // Read the tables
  FITS::FitsTable ft(path, FITS::ReadOnly, "DATA");
  transientTable = ft.extract();
  FITS::FitsTable ft2(path, FITS::ReadOnly, "INDEX");
  transientIndex = ft2.extract();

  // Build an index into the index:
  int expnum;
  for (int row=0; row<transientIndex.nrows(); row++) {
    transientIndex.readCell(expnum, "EXPNUM", row);
    findExpnum[expnum] = row;
  }
}

bool
TransientTable::hasExpnum(int expnum) const {
  return findExpnum.count(expnum);
}

Observation
TransientTable::getObservation(int objectID,
			       const Ephemeris& ephem,
			       ExposureTable& exposures) const {
  if (objectID<0 || objectID>transientTable.nrows())
    throw std::runtime_error("TransientTable::getObservation objectID out of range");

  Observation out;
  double ra, dec, sig;
  transientTable.readCell(ra,"RA", objectID);
  transientTable.readCell(dec,"DEC", objectID);
  out.radec = astrometry::SphericalICRS(ra*DEGREE, dec*DEGREE);

  int expnum;
  transientTable.readCell(expnum,"EXPNUM",objectID);

  transientTable.readCell(sig, "ERRAWIN_WORLD", objectID);
  double noiseCov = sig*sig*DEGREE*DEGREE;
  out.cov = exposures.atmosphereCov(expnum);
  out.cov(0,0) += noiseCov;
  out.cov(1,1) += noiseCov;

  out.tdb = ephem.mjd2tdb(exposures.mjd(expnum));
  out.observer = exposures.observatory(expnum);

  return out;
}  


bool
TransientTable::fillExposure(const Frame& frame,
			     Exposure* eptr) const {
  auto iter = findExpnum.find(eptr->expnum);
  if (iter == findExpnum.end())
    return false; // No entries for this.
  int index = iter->second;
  int begin,end;
  transientIndex.readCell(begin, "FIRST", index);
  transientIndex.readCell(end,  "LAST", index);
  float density;
  transientIndex.readCell(density, "DENSITY", index);
  eptr->detectionDensity = density / (DEGREE*DEGREE); // Change per sq deg to per sr

  std::vector<double> ra;
  transientTable.readCells(ra, "RA", begin, end);
  std::vector<double> dec;
  transientTable.readCells(dec, "DEC", begin, end);
  std::vector<double> sig;
  transientTable.readCells(sig, "ERRAWIN_WORLD", begin, end);
  std::vector<short int> ccd; 
  transientTable.readCells(ccd, "CCDNUM", begin, end);
    
  // Count number of usable detections
  BVector use(end-begin);
  for (int i=0; i<end-begin; i++)
    use[i] = sig[i]*DEGREE < MAXIMUM_TRANSIENT_POSITION_ERROR;
  
  // Fill in individual detections' properties
  // Their coordinates and sigma are in degrees.
  int nTransients = use.array().count();
  eptr->xy.resize(nTransients, 2);
  eptr->covXX.resize(nTransients);
  eptr->covXY.resize(nTransients);
  eptr->covYY.resize(nTransients);
  eptr->ccdnum.resize(nTransients);
  eptr->id.resize(nTransients);
  eptr->valid = BVector(nTransients,true); //Everyone is valid to start with.

  for (int i=0, j=0; j<end-begin; j++) {
    if (!use[j]) continue;
    eptr->id[i] = j+begin;
    eptr->ccdnum[i] = ccd[j];
    eptr->covXX[i] = sig[j]*sig[j]*DEGREE*DEGREE + eptr->atmosphereCov(0,0);
    eptr->covXY[i] = eptr->atmosphereCov(1,0);
    eptr->covYY[i] = sig[j]*sig[j]*DEGREE*DEGREE + eptr->atmosphereCov(1,1);
    astrometry::SphericalICRS radec(ra[j]*DEGREE, dec[j]*DEGREE);
    try {
      astrometry::Gnomonic xy(radec, frame.orient);
      double x,y;
      xy.getLonLat(x,y);
      eptr->xy(i,0) = x;
      eptr->xy(i,1) = y;
      eptr->valid[i] = true;
    } catch (astrometry::AstrometryError& a) {
      // Get here for projection out of gnomonic's hemisphere
      eptr->xy(i,0) = 0.;
      eptr->xy(i,1) = 0.;
      eptr->valid[i] = false;
    }
    i++;
  }
  return true;
}

////////////////////////////////////////////////////////////////////////
// Table for corner locations
////////////////////////////////////////////////////////////////////////
std::map<string,int>
CornerTable::detpos2ccdnum;

CornerTable::CornerTable(const string cornerFile) {
  string path = cornerFile;
  if (path.empty())  {
    // Find FITS file name from environment variable
    char *kpath = getenv(CORNERS_ENVIRON.c_str());
    if (kpath == NULL) 
      throw std::runtime_error("No path given for DES corner table FITS file");
    path = kpath;
  }
  // Read the table
  FITS::FitsTable ft(path, FITS::ReadOnly, 1);
  cornerTable = ft.extract();

  // Build an index from expnum to rows of the table
  vector<int> expnum;
  cornerTable.readCells(expnum,"EXPNUM");
  
  for (int row=0; row<expnum.size(); row++) {
    findExpnum.emplace(expnum[row],row); // ??
  }

  // DETPOS table:
  if (detpos2ccdnum.empty()) {
    detpos2ccdnum["N1"]= 32;
    detpos2ccdnum["N2"]= 33;
    detpos2ccdnum["N3"]= 34;
    detpos2ccdnum["N4"]= 35;
    detpos2ccdnum["N5"]= 36;
    detpos2ccdnum["N6"]= 37;
    detpos2ccdnum["N7"]= 38;
    detpos2ccdnum["N8"]= 39;
    detpos2ccdnum["N9"]= 40;
    detpos2ccdnum["N10"]= 41;
    detpos2ccdnum["N11"]= 42;
    detpos2ccdnum["N12"]= 43;
    detpos2ccdnum["N13"]= 44;
    detpos2ccdnum["N14"]= 45;
    detpos2ccdnum["N15"]= 46;
    detpos2ccdnum["N16"]= 47;
    detpos2ccdnum["N17"]= 48;
    detpos2ccdnum["N18"]= 49;
    detpos2ccdnum["N19"]= 50;
    detpos2ccdnum["N20"]= 51;
    detpos2ccdnum["N21"]= 52;
    detpos2ccdnum["N22"]= 53;
    detpos2ccdnum["N23"]= 54;
    detpos2ccdnum["N24"]= 55;
    detpos2ccdnum["N25"]= 56;
    detpos2ccdnum["N26"]= 57;
    detpos2ccdnum["N27"]= 58;
    detpos2ccdnum["N28"]= 59;
    detpos2ccdnum["N29"]= 60;
    detpos2ccdnum["N31"]= 62;
    detpos2ccdnum["S1"]= 25;
    detpos2ccdnum["S2"]= 26;
    detpos2ccdnum["S3"]= 27;
    detpos2ccdnum["S4"]= 28;
    detpos2ccdnum["S5"]= 29;
    detpos2ccdnum["S6"]= 30;
    detpos2ccdnum["S7"]= 31;
    detpos2ccdnum["S8"]= 19;
    detpos2ccdnum["S9"]= 20;
    detpos2ccdnum["S10"]= 21;
    detpos2ccdnum["S11"]= 22;
    detpos2ccdnum["S12"]= 23;
    detpos2ccdnum["S13"]= 24;
    detpos2ccdnum["S14"]= 13;
    detpos2ccdnum["S15"]= 14;
    detpos2ccdnum["S16"]= 15;
    detpos2ccdnum["S17"]= 16;
    detpos2ccdnum["S18"]= 17;
    detpos2ccdnum["S19"]= 18;
    detpos2ccdnum["S20"]=  8;
    detpos2ccdnum["S21"]=  9;
    detpos2ccdnum["S22"]= 10;
    detpos2ccdnum["S23"]= 11;
    detpos2ccdnum["S24"]= 12;
    detpos2ccdnum["S25"]=  4;
    detpos2ccdnum["S26"]=  5;
    detpos2ccdnum["S27"]=  6;
    detpos2ccdnum["S28"]=  7;
    detpos2ccdnum["S29"]=  1;
    detpos2ccdnum["S30"]=  2;
    detpos2ccdnum["S31"]=  3;
    detpos2ccdnum["N30"]= 61;
  }
}

bool
CornerTable::hasExpnum(int expnum) const {
  return findExpnum.count(expnum)>0;
}

bool
CornerTable::fillExposure(Exposure* eptr) const {
  auto range = findExpnum.equal_range(eptr->expnum);
  if (range.first==range.second)
    return false; // No entries for this.
  vector<double> ra;
  vector<double> dec;
  string detpos;
  for (auto ptr = range.first; ptr!=range.second; ++ptr) {
    int row = ptr->second;
    cornerTable.readCell(detpos, "DETPOS", row);
    stringstuff::stripWhite(detpos);
    cornerTable.readCell(ra,  "RA", row);
    cornerTable.readCell(dec,  "DEC", row);
  
    eptr->devices.push_back(detpos2ccdnum.at(detpos));

    vector<Point> vertices;
    double xx,yy;
    for (int i=0; i<4; i++) {
      // Convert vertex to local gnomonic system:
      astrometry::Gnomonic gn(astrometry::SphericalICRS(ra[i]*DEGREE, dec[i]*DEGREE),
			      eptr->localICRS);
      gn.getLonLat(xx,yy);
      vertices.push_back(Point(xx,yy));
    }
    eptr->deviceBoundsLocalICRS.push_back(ConvexPolygon(vertices));
  }
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
  bool operator()(Exposure* exptr) const {
    return exptr->tobs < value;
  }
private:
  double value;
};
class XSplit {
public:
  XSplit(double value_): value(value_) {}
  bool operator()(Exposure* exptr) const {
    return exptr->axis[0] < value;
  }
private:
  double value;
};
class YSplit {
public:
  YSplit(double value_): value(value_) {}
  bool operator()(Exposure* exptr) const {
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

list<Exposure*>
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

  list<Exposure*> out;
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

DESTree::DESTree(std::vector<Exposure*>& exposurePointers,
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

list<Exposure*>
DESTree::find(const Fitter& path) const {
  list<Exposure*> out;
  for (auto n : years) 
    out.splice(out.end(),n->find(path));
  return out;
}

DVector
Exposure::chisq(double x0, double y0, double covxx, double covyy, double covxy) const {
  if (xy.rows()==0) return DVector();
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
  
