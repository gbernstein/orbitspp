// Pieces of fitting routines

#include "Fitter.h"
#include "StringStuff.h"
#include "AstronomicalConstants.h"

#include <iomanip>

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
  dt.resize(n);
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
    dt[i] = tdb[i] - f.tdb0;
    astrometry::Matrix22 partials(0.);
    // Convert angles to our frame and get local partials for covariance
    projection.convertFrom(obs.radec, partials);
    double x,y;
    projection.getLonLat(x,y);
    thetaX[i] = x;
    thetaY[i] = y;
    // Get inverse covariance
    astrometry::Matrix22 cov(0.);
    cov(1,1) = cov(0,0) = pow(obs.sigma*ARCSEC,2);
    double det = cov(0,0)*cov(1,1) - cov(0,1)*cov(1,0);
    invcovXX[i] = cov(1,1)/det;
    invcovXY[i] = -cov(1,0)/det;
    invcovYY[i] = cov(0,0)/det;

    // Get observatory position in our frame.
    astrometry::CartesianICRS obspos = eph.observatory(obs.obscode, tdb[i]);
    projection3d.convertFrom(obspos);
    xE.col(i) = projection3d.getVector();

    /**/cerr << i << " " << std::fixed << setprecision(4) << dt[i]
	     << " " << setprecision(4) << xE(0,i) 
	     << " " << setprecision(4) << xE(1,i) 
	     << " " << setprecision(4) << xE(2,i)
	     << endl;
  }
}

void
Fitter::chooseFrame(int obsNumber) {
  if (obsNumber >= static_cast<int> (observations.size())) 
    throw runtime_error("Request for nonexistent obsNumber in Fitter::chooseFrame()");
  if (observations.empty())
    throw runtime_error("Must have observations to call Fitter::chooseFrame()");
  if (obsNumber < 0) {
    double mjd0 = observations.front().mjd;
    double tsum=0;
    for (auto& obs : observations)
      tsum += obs.mjd - mjd0;
    double target = mjd0 + tsum / observations.size();
    double minDT = 1e9;
    for (int i=0; i<observations.size(); i++) {
      if (abs(observations[i].mjd-target)<minDT) {
	obsNumber = i;
	minDT = abs(observations[i].mjd-target);
      }
    }
  }
  /**/cerr << "Chose observation number " << obsNumber << endl;

  // Now construct ecliptic-aligned frame and set up all observations
  const MPCObservation& obs = observations[obsNumber];
  astrometry::Orientation orient(obs.radec);
  orient.alignToEcliptic();
  double tdb0 = eph.mjd2tdb(obs.mjd);
  astrometry::CartesianICRS origin = eph.observatory(obs.obscode,tdb0);
  setFrame(Frame(origin, orient, tdb0));
}


    
