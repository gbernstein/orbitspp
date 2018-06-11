// Pieces of fitting routines

#include "Fitter.h"
#include "StringStuff.h"
#include "AstronomicalConstants.h"

#include <iomanip>

using namespace std;
using namespace orbits;

/**/
template<class T>
void
write6(const T& v, std::ostream& os, int precision=6) {
  stringstuff::StreamSaver ss(os);  // Cache stream state
  os << std::fixed << std::showpos << std::setprecision(precision);
  for (int i=0; i<6; i++)
    os << v[i] << " ";
  return;  // Stream state restored on destruction of ss
};
  
Fitter::Fitter(const Ephemeris& eph_,
	       Gravity grav_): eph(eph_),
			       fullTrajectory(nullptr),
			       grav(grav_),
			       energyConstraintFactor(0.),
			       gammaPriorSigma(0.),
			       b(6),
			       A(6,6) {}
void
Fitter::readMPCObservations(istream& is) {
  string line;
  while (stringstuff::getlineNoComment(is, line)) {
    auto obs = mpc2Observation(MPCObservation(line),eph);
    observations.push_back(obs);
  };
}

// Routine used to put observations into time order
bool
observationTimeOrder(const Observation& m1,
		     const Observation& m2)
{
  return m1.tdb < m2.tdb;
}

void
Fitter::setFrame(const Frame& f_) {
  f = f_;
  // Recalculate all coordinates in this frame;
  int n = observations.size();
  // Put observations into time order.
  std::sort(observations.begin(), observations.end(),
    observationTimeOrder);
  tObs.resize(n);
  tEmit.resize(n);
  tdbEmit.resize(n);
  thetaX.resize(n);
  thetaY.resize(n);
  thetaXModel.resize(n);
  thetaYModel.resize(n);
  invcovXX.resize(n);
  invcovXY.resize(n);
  invcovYY.resize(n);
  xE.resize(n,3);
  xGrav.resize(n,3);
  dThetaXdABG.resize(n,6);
  dThetaYdABG.resize(n,6);

  astrometry::Gnomonic projection(f.orient);  
  astrometry::CartesianCustom projection3d(f);
  for (int i=0; i<n; i++) {
    const Observation& obs = observations[i];
    // Set up time variables.  Set light travel time to zero
    tObs[i] = obs.tdb - f.tdb0;
    tEmit[i] = tObs[i];
    tdbEmit[i] = obs.tdb;
    
    // Convert angles to our frame and get local partials for covariance
    astrometry::Matrix22 partials(0.);
    projection.convertFrom(obs.radec, partials);
    double x,y;
    projection.getLonLat(x,y);
    thetaX[i] = x;
    thetaY[i] = y;
    // Alter the partials for the cos(dec) factor in derivatives
    obs.radec.getLonLat(x,y);
    partials(0,0) /= cos(y);
    partials(1,0) /= cos(y);
    // Get inverse covariance
    astrometry::Matrix22 cov = partials * obs.cov * partials.transpose();
    double det = cov(0,0)*cov(1,1) - cov(0,1)*cov(1,0);
    invcovXX[i] = cov(1,1)/det;
    invcovXY[i] = -cov(1,0)/det;
    invcovYY[i] = cov(0,0)/det;

    // Get observatory position in our frame.
    astrometry::CartesianICRS obspos = obs.observer;
    projection3d.convertFrom(obspos);
    xE.row(i) = projection3d.getVector();

    // Initialize gravitational effect to zero
    xGrav.setZero();
    
    /**cerr << i << " " << std::fixed << setprecision(4) << dt[i]
	     << " " << setprecision(4) << xE(0,i) 
	     << " " << setprecision(4) << xE(1,i) 
	     << " " << setprecision(4) << xE(2,i)
	     << " partials: " << partials
	     << endl;
    /**cerr << " Cov: " << std::scientific << cov << endl;
    /**/
  }
}

void
Fitter::chooseFrame(int obsNumber) {
  if (obsNumber >= static_cast<int> (observations.size())) 
    throw runtime_error("Request for nonexistent obsNumber in Fitter::chooseFrame()");
  if (observations.empty())
    throw runtime_error("Must have observations to call Fitter::chooseFrame()");
  if (obsNumber < 0) {
    double tdb0 = observations.front().tdb;
    double tsum=0;
    for (auto& obs : observations)
      tsum += obs.tdb - tdb0;
    double target = tdb0 + tsum / observations.size();
    double minDT = 1e9;
    for (int i=0; i<observations.size(); i++) {
      if (abs(observations[i].tdb-target)<minDT) {
	obsNumber = i;
	minDT = abs(observations[i].tdb-target);
      }
    }
  }
  
  // Now construct ecliptic-aligned frame and set up all observations
  const Observation& obs = observations[obsNumber];
  astrometry::Orientation orient(obs.radec);
  orient.alignToEcliptic();
  double tdb0 = obs.tdb;
  astrometry::CartesianICRS origin = obs.observer;
  setFrame(Frame(origin, orient, tdb0));
}

void
Fitter::calculateOrbitDerivatives() {
  int nobs = observations.size();

  // Calculate positions
  Vector denom = Vector(nobs,1.) + abg[ABG::GDOT]*tEmit
    + abg[ABG::G]*(xGrav.col(2) - xE.col(2));
  // ??? add TMV version
  denom = denom.cwiseInverse();  // Denom is now 1/(z*gamma)
  thetaXModel = abg[ABG::ADOT]*tEmit
    + abg[ABG::G]*(xGrav.col(0) - xE.col(0));
  thetaXModel.array() += abg[ABG::A];
  thetaYModel = abg[ABG::BDOT]*tEmit
    + abg[ABG::G]*(xGrav.col(1) - xE.col(1));
  thetaYModel.array() += abg[ABG::B];
  thetaXModel.array() *= denom.array();
  thetaYModel.array() *= denom.array();

  // Calculate derivatives
  dThetaXdABG.setZero();
  dThetaYdABG.setZero();
  dThetaXdABG.col(ABG::A).setOnes();
  dThetaYdABG.col(ABG::B).setOnes();
  dThetaXdABG.col(ABG::ADOT) = tEmit;
  dThetaYdABG.col(ABG::BDOT) = tEmit;
  dThetaXdABG.col(ABG::G) = xGrav.col(0) - xE.col(0)
    + thetaXModel.cwiseProduct(xE.col(2)-xGrav.col(2));
  dThetaYdABG.col(ABG::G) = xGrav.col(1) - xE.col(1)
    + thetaYModel.cwiseProduct(xE.col(2)-xGrav.col(2));
  dThetaXdABG.col(ABG::GDOT) -= thetaXModel.cwiseProduct(tEmit);
  dThetaYdABG.col(ABG::GDOT) -= thetaYModel.cwiseProduct(tEmit);

  dThetaXdABG.applyOnTheLeft(denom.asDiagonal());
  dThetaYdABG.applyOnTheLeft(denom.asDiagonal());

  return;
}

void
Fitter::iterateTimeDelay() {
  // Use current orbit to update light-travel time
  int nobs = observations.size();
  // Light-travel time is distance/(speed of light)

  Matrix xyz = abg.getXYZ(tEmit) + xGrav - xE;
  // Get sqrt(xyz.xyz):
  Vector dist = xyz.cwiseProduct(xyz).rowwise().sum().cwiseSqrt();
  tEmit = tObs - dist/SpeedOfLightAU;
  tdbEmit.array() = tEmit.array() + f.tdb0;
  return;
}

void
Fitter::calculateGravity() {
  // Calculate the non-inertial terms in trajectory of object
  // using current ABG and gravity settings
  if (grav==Gravity::INERTIAL) {
    xGrav.setZero();
    return;
  }
  // Kill any old trajectory
  if (fullTrajectory) {
    delete fullTrajectory;
    fullTrajectory = nullptr;
  }
  // Get initial condition for integrator from abg
  astrometry::Vector3 x0;
  astrometry::Vector3 v0;
  abg.getState(0., x0, v0);

  // Rotate to ICRS for integrator
  State s0;
  s0.x = astrometry::CartesianICRS(f.orient.m().transpose() * x0
				   + f.origin.getVector());
  s0.v = astrometry::CartesianICRS(f.orient.m().transpose() * v0);
  s0.tdb = f.tdb0;
  
  // Integrate - Trajectory returns 3xN matrix
  fullTrajectory = new Trajectory(eph, s0, grav);
  astrometry::DMatrix xyz = fullTrajectory->position(tdbEmit);
  // Convert back to Fitter reference frame and subtract inertial motion

  xyz.colwise() -= f.origin.getVector();  // ?? Problem mixing static size with dynamic ??
  xGrav = (f.orient.m() * xyz).transpose() - abg.getXYZ(tEmit);
  /**{
    stringstuff::StreamSaver ss(cerr);
    cerr << std::fixed << std::showpos << std::setprecision(4);
    for (int i=0; i<xGrav.rows(); i++)
      cerr << i << " " << tEmit[i]
	   << " " << xE(i,0) << " " << xE(i,1) << " " << xE(i,2)
	   << " " << xGrav(i,0) << " " << xGrav(i,1) << " " << xGrav(i,2)
	   << endl;
    cerr << "Done" << endl;
    } **/
  return;
}

void
Fitter::calculateChisqDerivatives() {
  // Assumes that model values are up to date
  Vector dx = thetaX - thetaXModel;
  Vector dy = thetaY - thetaYModel;
  chisq = dx.transpose() * (invcovXX.asDiagonal() * dx);
  chisq += 2.* dx.transpose() * (invcovXY.asDiagonal() * dy);
  chisq += dy.transpose() * (invcovYY.asDiagonal() * dy);
  Matrix tmp = invcovXX.asDiagonal() * dThetaXdABG;
  b = tmp.transpose() * dx;
  A = dThetaXdABG.transpose() * tmp;
  tmp = invcovXY.asDiagonal() * dThetaXdABG;
  b += tmp.transpose() * dy;
  A += dThetaYdABG.transpose() * tmp;
  tmp = invcovXY.asDiagonal() * dThetaYdABG;
  b += tmp.transpose() * dx;
  A += dThetaXdABG.transpose() * tmp;
  tmp = invcovYY.asDiagonal() * dThetaYdABG;
  b += tmp.transpose() * dy;
  A += dThetaYdABG.transpose() * tmp;

  // Add gamma constraint, if any
  if (gammaPriorSigma > 0.) {
    double w = pow(gammaPriorSigma, -2);
    double dg = (gamma0-abg[ABG::G]);
    chisq += dg * w * dg;
    b[ABG::G] += w * dg;
    A(ABG::G,ABG::G) += w;
  }

  // Add derivatives from binding prior - crude one on GDTO??
  {
    const double BIND_PRIOR = 1.; // Adjust penalty for high velocity
    double w = BIND_PRIOR / (2. * GM * pow(abg[ABG::G],3.));
    chisq += abg[ABG::GDOT] * w * abg[ABG::GDOT];
    b[ABG::GDOT] += -w * abg[ABG::GDOT];
    A(ABG::GDOT,ABG::GDOT) += w;
  }
}
  
void
Fitter::setLinearOrbit(double nominalGamma) {

  // Assume that current gravity is good (probably zero)
  calculateOrbitDerivatives();
  calculateChisqDerivatives();

  // Extract parameters other than GDOT, which we assume is last
  Vector blin = b.subVector(0, ABG::GDOT);
  Matrix Alin = A.subMatrix(0, ABG::GDOT, 0, ABG::GDOT);
  
  
  // Solve (check degeneracies??)
  auto llt = Alin.llt();
  Vector answer = llt.solve(blin);

  // Transfer answer to instance members
  abg.subVector(0,ABG::GDOT) += answer;
  abg[ABG::GDOT] = 0.;
}

void
Fitter::newtonFit(double chisqTolerance) {
  // Assuming that we are coming into this with a good starting point
  //**iterateTimeDelay();
  calculateGravity();
  calculateOrbitDerivatives();
  calculateChisqDerivatives();
  const int MAX_ITERATIONS = 20;

  // Do a series of iterations with time delay and gravity fixed
  int iterations = 0;
  double oldChisq;
  do {
    oldChisq = chisq;
    // Solve (check degeneracies??)
    auto llt = A.llt();
    abg += llt.solve(b);
    /**/{
      /**/cerr << "Iteration " << iterations <<endl;
      auto dd = llt.solve(b);
      cerr << " change: ";
      write6(dd,cerr);
      cerr << endl << " abg: ";
      write6(abg,cerr);
    }
    calculateOrbitDerivatives();
    calculateChisqDerivatives();
    /**/cerr << endl << " New chisq: " << chisq << endl;
    iterations++;
  } while (iterations < MAX_ITERATIONS && abs(chisq-oldChisq) > chisqTolerance);
  if (iterations >= MAX_ITERATIONS)
    throw std::runtime_error("Fitter::newtonFit exceeded max iterations");
   
  // Then converge with time delay and gravity recalculated
  //**iterateTimeDelay();
  calculateGravity();
  calculateOrbitDerivatives();
  calculateChisqDerivatives();
  iterations = 0;
  do {
    oldChisq = chisq;
    // Solve (check degeneracies??)
    auto llt = A.llt();
    abg += llt.solve(b);
    /**/{
      /**/cerr << "Gravity iteration " << iterations <<endl;
      auto dd = llt.solve(b);
      cerr << " change: ";
      write6(dd,cerr);
      cerr << endl << " abg: ";
      write6(abg,cerr);
      cerr << endl;
    }
    //**iterateTimeDelay();
    calculateGravity();
    calculateOrbitDerivatives();
    calculateChisqDerivatives();
    /**/cerr << endl << " New chisq: " << chisq << endl;
    iterations++;
  } while (iterations < MAX_ITERATIONS && abs(chisq-oldChisq) > chisqTolerance);
  if (iterations >= MAX_ITERATIONS)
    throw std::runtime_error("Fitter::newtonFit exceeded max fine iterations");
  return;
}

void
Fitter::printResiduals(std::ostream& os) const {
  stringstuff::StreamSaver ss(os);
  os << "# Residuals: " << endl << std::fixed << std::setprecision(1);
  double chitot = 0;
  double chi;
  Vector dx = thetaX - thetaXModel;
  Vector dy = thetaY - thetaYModel;
  os << "# N   T    dx    dy   chisq" << endl;
  os << "#   (days) (arcsecond)   " << endl;
  for (int i=0; i<dx.size(); i++) {
    chi = dx[i]*dx[i]*invcovXX[i] + 2.*dx[i]*dy[i]*invcovXY[i] + dy[i]*dy[i]*invcovYY[i];
    os << std::setw(3) << i << "  " << std::showpos << std::setw(6) << tObs[i]/DAY
       << " " << dx[i]/ARCSEC << " " << dy[i]/ARCSEC
       << std::noshowpos << " " << chi << endl;
    chitot += chi;
  }
  os << " chisq w/o priors " << chitot << " w/priors: " << chisq << endl;
}

void
Fitter::printCovariance(std::ostream& os) const {
  stringstuff::StreamSaver ss(os);
  os << "# ABG Covariance " << endl << std::scientific << std::setprecision(4) << std::showpos;
  Matrix c = A.inverse();
  for (int i=0; i<6; i++) {
    for (int j=0; j<6; j++)
      os << std::setw(10) << c(i,j) << " ";
    os << endl;
  }
}
