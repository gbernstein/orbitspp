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

void
Fitter::setLinearOrbit() {

  // Create derivatives of signal wrt abg
  int n = thetaX.size();
  Vector ones(n,1.);
  Matrix dxdp(n,5,0.);
  Matrix dydp(n,5,0.);
  /**/cerr << "here 1" << endl;
  dxdp.col(ABG::A) = ones;
  dxdp.col(ABG::ADOT) = dt;
  dxdp.col(ABG::G) = -xE.row(0).transpose();
  dydp.col(ABG::B) = ones;
  dydp.col(ABG::BDOT) = dt;
  dydp.col(ABG::G) = -xE.row(1).transpose();

  /**/cerr << "here 2" << endl;
  // Create A matrix, b vector
  Matrix tmpx = invcovXX.asDiagonal() * dxdp + invcovXY.asDiagonal() * dydp;
  Matrix tmpy = invcovXY.asDiagonal() * dxdp + invcovYY.asDiagonal() * dydp;
  Vector blin = tmpx.transpose() * thetaX + tmpy.transpose() * thetaY;
  Matrix Alin = dxdp.transpose() * tmpx + dydp.transpose() * tmpy;
  
  /**/cerr << "here 3" << endl;
  // Add gamma constraint, if any
  if (gammaPriorSigma > 0.) {
    double w = pow(gammaPriorSigma, -2);
    blin[ABG::G] += w * gamma0;
    Alin(ABG::G,ABG::G) += w;
  }
  
  /**/cerr << "here 4" << endl;
  // Solve (check degeneracies??)
  auto llt = Alin.llt();
  Vector answer = llt.solve(blin);
  /**/cerr << "here 5" << endl;
  abg[ABG::A] = answer[ABG::A];
  abg[ABG::B] = answer[ABG::B];
  abg[ABG::G] = answer[ABG::G];
  abg[ABG::ADOT] = answer[ABG::ADOT];
  abg[ABG::BDOT] = answer[ABG::BDOT];
  abg[ABG::GDOT] = 0.;

  /**/
  Matrix cov = Alin.inverse();
  Vector sigma = cov.diagonal().cwiseSqrt();
  Matrix corr = cov.cwiseQuotient( sigma * sigma.transpose());
  /**/cerr << scientific << "sigmas:\n" << sigma << endl;
  /**/cerr << "solution covariance:\n" << corr << endl;
}
