// Pieces of fitting routines

#include "Fitter.h"
#include "StringStuff.h"

using namespace std;
using namespace orbits;

Fitter::Fitter(const Ephemeris& eph_): eph(eph_),
				       inertialTrajectory(nullptr),
				       fullTrajectory(nullptr),
				       useGiants(false),
				       energyConstraintFactor(0.),
				       gammaPriorSigma(0.) {}
void
Fitter::readObservations(istream& is) {
  string line;
  while (stringstuff::getlineNoComment(is, line)) {
    observations.push_back( MPCObservation(line));
  };
}

// Routine used to put observations into time order
bool
observationTimeOrder(const MPCObservation& m1,
		     const MPCObservation& m2)
{
  return m1.mjd < m2.mjd;
}

void
Fitter::setFrame(const Frame& f_) {
  f = f_;
  // Recalculate all coordinates in this frame;
  int n = observations.size();
  // Put observations into time order.
  std::sort(observations.begin(), observations.end(),
	    observationTimeOrder);
  tdb.resize(n);
  thetaX.resize(n);
  thetaY.resize(n);
  invcovXX.resize(n);
  invcovXY.resize(n);
  invcovYY.resize(n);
  xE.resize(3,n);
  xGrav.resize(3,n);

  astrometry::Gnomonic projection(f.orient, true);  // share Orientation
  astrometry::CartesianCustom projection3d(f);
  for (int i=0; i<n; i++) {
    const MPCObservation& obs = observations[i];
    tdb[i] = eph.mjd2tdb(obs.mjd);
    astrometry::Matrix22 partials(0.);
    // Convert angles to our frame and get local partials for covariance
    projection.convertFrom(obs.radec, partials);
    double x,y;
    projection.getLonLat(x,y);
    thetaX[i] = x;
    thetaY[i] = y;
    // Get inverse covariance
    astrometry::Matrix22 cov(0.);
    cov(1,1) = cov(0,0) = obs.sigma*obs.sigma;
    cov = partials.transpose() * cov * partials;
    double det = cov(0,0)*cov(1,1) - cov(0,1)*cov(1,0);
    invcovXX[i] = cov(1,1)/det;
    invcovXY[i] = -cov(1,0)/det;
    invcovYY[i] = cov(0,0)/det;

    // Get observatory position in our frame.
    astrometry::CartesianICRS obspos = eph.observatory(obs.obscode, tdb[i]);
    projection3d.convertFrom(obspos);
    xE.row(i) = projection3d.getVector().transpose();
  }
}

    
