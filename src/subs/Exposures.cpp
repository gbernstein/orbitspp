#include "Exposures.h"
#include "FitsTable.h"
#include "Trajectory.h"

using namespace orbits;

std::vector<Exposure>
selectExposures(const Frame& frame,   // Starting coordinates, time
		const Ephemeris& ephem,  
		double gamma0,        // Center and width of range 
		double dGamma,        // of gamma to cover
		double searchRadius,  // Range of starting coords to cover
		string exposureFile,
		double fieldRadius) { // Circumscribed field radius (degrees)


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
