// Pieces of fitting routines

#include "Fitter.h"
#include "StringStuff.h"
#include "AstronomicalConstants.h"
#include "Elements.h"

#include <iomanip>

using namespace std;
using namespace orbits;

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
			       bindingConstraintFactor(0.),
			       gammaPriorSigma(0.),
			       frameIsSet(false),
			       abgIsFit(false),
			       positionsAreValid(false),
			       positionDerivsAreValid(false),
			       chisqIsValid(false),
			       chisqDerivsAreValid(false) {} 

void
Fitter::setABG(const ABG& abg_, const ABGCovariance& cov_) {
  if (!frameIsSet)
    throw std::runtime_error("ERROR: Fitter::setABG called before setFrame");
  abg = abg_;
  A = cov_.inverse();
  newABG(); // update status flags
  // Set up Trajectory if using gravity
  if (grav != Gravity::INERTIAL) createTrajectory();
}

void
Fitter::addObservation(const Observation& obs) {
  if (frameIsSet)
    throw std::runtime_error("ERROR: Fitter::addObservations called after setFrame");
      observations.push_back(obs);
}

void
Fitter::readMPCObservations(istream& is) {
  if (frameIsSet)
    throw std::runtime_error("ERROR: Fitter::readMPCObservations called after setFrame");
  string line;
  while (stringstuff::getlineNoComment(is, line)) {
    auto obs = mpc2Observation(MPCObservation(line),eph);
    observations.push_back(obs);
  };
}

void
Fitter::resizeArrays(int n) {
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
}

void
Fitter::setFrame(const Frame& f_) {
  if (frameIsSet)
    throw std::runtime_error("ERROR: Fitter::setFrame called after already set");
  f = f_;
  frameIsSet = true;
  newData(); // Reset results flags

  // Recalculate all coordinates in this frame;
  int n = observations.size();
  resizeArrays(n);

  astrometry::Gnomonic projection(f.orient);  
  for (int i=0; i<n; i++) {
    const Observation& obs = observations[i];
    // Set up time variables.  Set light travel time to zero
    tObs[i] = obs.tdb - f.tdb0;
    tEmit[i] = tObs[i];
    tdbEmit[i] = obs.tdb;
    
    // Convert angles to our frame and get local partials for covariance
    Matrix22 partials(0.);
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
    Matrix22 cov = partials * obs.cov * partials.transpose();
    double det = cov(0,0)*cov(1,1) - cov(0,1)*cov(1,0);
    invcovXX[i] = cov(1,1)/det;
    invcovXY[i] = -cov(1,0)/det;
    invcovYY[i] = cov(0,0)/det;

    // Get observatory position in our frame.
    xE.row(i) = f.fromICRS(obs.observer.getVector());

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

// Replace observations with data that is already in desired Frame
void
Fitter::setObservationsInFrame(const DVector& tObs_,    // TDB since reference time
			       const DVector& thetaX_,  // Observed positions
			       const DVector& thetaY_,
			       const DVector& covXX_,   // Covariance of observed posns
			       const DVector& covYY_,
			       const DVector& covXY_,
			       const DMatrix& xE_) {    // Observatory posn at observations
  if (!frameIsSet)
    throw std::runtime_error("ERROR: Fitter::setObservationsInFrame() was "
			     "called before setFrame()");
  int n = tObs_.size();
  if (thetaX_.size()!=n ||
      thetaY_.size()!=n ||
      covXX_.size()!=n ||
      covYY_.size()!=n ||
      covXY_.size()!=n ||
      xE_.rows() !=n ||
      xE_.cols() !=3)
    throw std::runtime_error("ERROR: Mismatched array sizes in "
			     "Fitter::setObservationsInFrame()");

  resizeArrays(n);
  newData(); // Reset results flags

  tObs = tObs_;
  tEmit = tObs_; // No light-travel correction to start
  tdbEmit = tEmit.array() + f.tdb0;
  thetaX = thetaX_;
  thetaY = thetaY_;
  DVector det = covXX_.array()*covYY_.array() - covXY_.array()*covXY_.array();
  invcovXX = covYY_.array() / det.array();
  invcovYY = covXX_.array() / det.array();
  invcovXY = -covXY_.array() / det.array();
  xE = xE_;
  xGrav.setZero();
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
Fitter::calculateOrbit(bool doDerivatives) {
  int nobs = tEmit.size();

  // Is is already done?
  if (doDerivatives && positionsAreValid && positionDerivsAreValid) return;
  if (!doDerivatives && positionsAreValid) return;
  
  // Calculate positions
  DVector denom = DVector(nobs,1.) + abg[ABG::GDOT]*tEmit
    + abg[ABG::G]*(xGrav.col(2) - xE.col(2));

  denom = denom.cwiseInverse();  // Denom is now 1/(z*gamma)
  thetaXModel = abg[ABG::ADOT]*tEmit
    + abg[ABG::G]*(xGrav.col(0) - xE.col(0));
  thetaXModel.array() += abg[ABG::A];
  thetaYModel = abg[ABG::BDOT]*tEmit
    + abg[ABG::G]*(xGrav.col(1) - xE.col(1));
  thetaYModel.array() += abg[ABG::B];
  thetaXModel.array() *= denom.array();
  thetaYModel.array() *= denom.array();
  positionsAreValid = true;

  if (doDerivatives) {
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
    positionDerivsAreValid = true;
  }
  
  // Downstream quantities are invalid now
  chisqIsValid = chisqDerivsAreValid = false; 
  return;
}

void
Fitter::iterateTimeDelay() {
  // Use current orbit to update light-travel time
  int nobs = tObs.size();
  // Light-travel time is distance/(speed of light)

  DMatrix xyz = abg.getXYZ(tEmit) + xGrav - xE;
  // Get sqrt(xyz.xyz):
  DVector dist = xyz.cwiseProduct(xyz).rowwise().sum().cwiseSqrt();
  if ( (dist.array()>1e4).any())
    throw std::runtime_error("ERROR: distance >10,000 AU in Fitter::iterateTimeDelay()");
  tEmit = tObs - dist/SpeedOfLightAU;
  tdbEmit.array() = tEmit.array() + f.tdb0;

  // This invalidates the solution
  // Solutions no longer match gravity
  newABG();
  abgIsFit = false;
  
  return;
}

void
Fitter::createTrajectory() {
  // Make a new Trajectory integrator
  // based on current ABG.

  // Kill any old trajectory
  if (fullTrajectory) {
    delete fullTrajectory;
    fullTrajectory = nullptr;
  }
  // Get initial condition for integrator from abg
  Vector3 x0;
  Vector3 v0;
  abg.getState(0., x0, v0);

  // Rotate to ICRS for integrator
  State s0;
  s0.x = astrometry::CartesianICRS(f.toICRS(x0));
  s0.v = astrometry::CartesianICRS(f.toICRS(v0,true));
  s0.tdb = f.tdb0;
  
  fullTrajectory = new Trajectory(eph, s0, grav);
}

void
Fitter::calculateGravity() {
  // Calculate the non-inertial terms in trajectory of object
  // using current ABG and gravity settings
  if (grav==Gravity::INERTIAL) {
    xGrav.setZero();
    return;
  }
  createTrajectory();
  
  DMatrix xyz = fullTrajectory->position(tdbEmit);

  // Convert back to Fitter reference frame and subtract inertial motion
  // Note Frame and Trajectory use 3 x n, we are using n x 3 here.
  xGrav = f.fromICRS(xyz) - abg.getXYZ(tEmit);
  if (false) /**/{
    DMatrix x1 = f.fromICRS(xyz);
    DMatrix x2 = abg.getXYZ(tEmit);
    stringstuff::StreamSaver ss(cerr);
    cerr << std::fixed << std::showpos << std::setprecision(4);
    for (int i=0; i<xGrav.rows(); i++) {
      cerr << i << ": " << tEmit[i] << endl;
      cerr << " " << x1(i,0) << " " << x1(i,1) << " " << x1(i,2) << endl;
      cerr << " " << x2(i,0) << " " << x2(i,1) << " " << x2(i,2) << endl;
      cerr << " " << xGrav(i,0) << " " << xGrav(i,1) << " " << xGrav(i,2) << endl;
    } /**/
  }

  // Solutions no longer match gravity
  newABG();
  abgIsFit = false;

  return;
}

void
Fitter::calculateChisq(bool doDerivatives) {
  if (tObs.size() < 2)
    throw std::runtime_error("ERROR: Fitter::calculateChisq called with <2 observations");

  // Is it already done?
  if (doDerivatives && chisqIsValid && chisqDerivsAreValid) return;
  if (!doDerivatives && chisqIsValid) return;
  
  // Make sure that model positions are current (and derivatives, if needed)
  calculateOrbit(doDerivatives);
    
  DVector dx = thetaX - thetaXModel;
  DVector dy = thetaY - thetaYModel;
  chisq = dx.transpose() * (invcovXX.asDiagonal() * dx);
  chisq += 2.* dx.transpose() * (invcovXY.asDiagonal() * dy);
  chisq += dy.transpose() * (invcovYY.asDiagonal() * dy);

  if (doDerivatives) {
    DMatrix tmp = invcovXX.asDiagonal() * dThetaXdABG;
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
  }
  
  // Add gamma constraint, if any
  if (gammaPriorSigma > 0.) {
    double w = pow(gammaPriorSigma, -2);
    double dg = (gamma0-abg[ABG::G]);
    chisq += dg * w * dg;
    if (doDerivatives) {
      b[ABG::G] += w * dg;
      A(ABG::G,ABG::G) += w;
    }
  }

  // Add derivatives from binding prior - crude one
  // which pushes to v=0, plunging perihelion.
  if (bindingConstraintFactor > 0.) {
    // Do not allow negative chisq:
    double w = bindingConstraintFactor / (2. * GM * pow(abs(abg[ABG::G]),3.));
    chisq += w * (abg[ABG::ADOT] * abg[ABG::ADOT]
		  + abg[ABG::BDOT] * abg[ABG::BDOT]
		  + abg[ABG::GDOT] * abg[ABG::GDOT]);

    if (doDerivatives) {
      double dampFactor = 0.5; // This adjustment prevents the solution from oscillating
      // due to the energy constraint.
      b[ABG::ADOT] += -w * abg[ABG::ADOT] * dampFactor;
      b[ABG::BDOT] += -w * abg[ABG::BDOT] * dampFactor;
      b[ABG::GDOT] += -w * abg[ABG::GDOT] * dampFactor;
      A(ABG::ADOT,ABG::ADOT) += w;
      A(ABG::BDOT,ABG::BDOT) += w;
      A(ABG::GDOT,ABG::GDOT) += w;
    }
  }

  chisqIsValid = true;
  chisqDerivsAreValid = doDerivatives;
}

double
Fitter::getChisq() {
  calculateChisq(false);
  return chisq;
}

void
Fitter::setLinearOrbit() {

  // Get positions/derivatives for current ABG
  calculateChisq(true);

  // Extract parameters other than GDOT, which we assume is last
  DVector blin = b.subVector(0, ABG::GDOT);
  DMatrix Alin = A.subMatrix(0, ABG::GDOT, 0, ABG::GDOT);
  
  // Solve (check degeneracies??)
  auto llt = Alin.llt();
  DVector answer = llt.solve(blin);

  // Transfer answer to instance members
  abg.subVector(0,ABG::GDOT) += answer;
  abg[ABG::GDOT] = 0.;

  newABG(); // Reset results flags
  // But don't consider linear solution to be valid
  abgIsFit = false;

  // Update positions, chisq - no derivatives
  calculateChisq(false);
}

void
Fitter::abgSanityCheck() const {
  const double MAX_GAMMA = 1.; // Any closer than this is failure
  const double MAX_KE = 10.; // Failure when crude |KE/PE| exceeds this

  if (abs(abg[ABG::G]) > MAX_GAMMA) {
    FormatAndThrow<std::runtime_error>()
      << "ERROR: Fitter gamma grows too large: " << abg[ABG::G];
  }
  double ke = (abg[ABG::ADOT] * abg[ABG::ADOT]
	       + abg[ABG::BDOT] * abg[ABG::BDOT]
	       + abg[ABG::GDOT] * abg[ABG::GDOT]) / (2. * GM * pow(abs(abg[ABG::G]),3.));
  if (ke > MAX_KE) {
    FormatAndThrow<std::runtime_error>()
      << "ERROR: Fitter |KE/PE| grows too large: " << ke;
  }
}

void
Fitter::newtonFit(double chisqTolerance, bool dump) {
  // Assuming that we are coming into this with a good starting point
  iterateTimeDelay();
  calculateGravity();
  calculateChisq(true);
  const int MAX_ITERATIONS = 20;

  // cancel validity of results until re-converged
  abgIsFit = false;
  
  // Do a series of iterations with time delay and gravity fixed
  int iterations = 0;
  double oldChisq;
  do {
    oldChisq = chisq;
    // Solve (check degeneracies??)
    auto llt = A.llt();
    abg += llt.solve(b);

    if (dump) {
      cerr << "Iteration " << iterations <<endl;
      auto dd = llt.solve(b);
      cerr << " change: ";
      write6(dd,cerr);
      cerr << endl << " abg:    ";
      write6(abg,cerr);
    }

    // Quit if iterations have gone completely awry
    abgSanityCheck();

    // Recalculate with new ABG
    calculateChisq(true);
    if (dump) cerr << endl << " New chisq: " << chisq << endl;
    iterations++;
  } while (iterations < MAX_ITERATIONS && abs(chisq-oldChisq) > chisqTolerance);
  if (iterations >= MAX_ITERATIONS)
    throw std::runtime_error("Fitter::newtonFit exceeded max iterations");
   
  // Then converge with time delay and gravity recalculated
  iterateTimeDelay();
  calculateGravity();
  calculateChisq(true);
  iterations = 0;
  do {
    oldChisq = chisq;
    // Solve (check degeneracies??)
    auto llt = A.llt();
    abg += llt.solve(b);
    if (dump) {
      cerr << "Gravity iteration " << iterations <<endl;
      auto dd = llt.solve(b);
      cerr << " change: ";
      write6(dd,cerr);
      cerr << endl << " abg:    ";
      write6(abg,cerr);
      cerr << endl;
    }

    // Quit if iterations have gone completely awry
    abgSanityCheck();

    // Update with new ABG
    iterateTimeDelay();
    calculateGravity();
    calculateChisq(true);
    if (dump) cerr << endl << " New chisq: " << chisq << endl;
    iterations++;
  } while (iterations < MAX_ITERATIONS && abs(chisq-oldChisq) > chisqTolerance);
  if (iterations >= MAX_ITERATIONS)
    throw std::runtime_error("Fitter::newtonFit exceeded max fine iterations");
  abgIsFit = true;
  return;
}

void
Fitter::printResiduals(std::ostream& os) {
  // Update residuals and chisq to current ABG if needed
  calculateChisq(false);

  stringstuff::StreamSaver ss(os);
  os << "# Residuals: " << endl << std::fixed << std::setprecision(2);
  double chitot = 0;
  double chi;
  DVector dx = thetaX - thetaXModel;
  DVector dy = thetaY - thetaYModel;
  os << "# N    T      dx      dy    chisq" << endl;
  os << "#    (days)  (arcsecond)   " << endl;
  for (int i=0; i<dx.size(); i++) {
    chi = dx[i]*dx[i]*invcovXX[i] + 2.*dx[i]*dy[i]*invcovXY[i] + dy[i]*dy[i]*invcovYY[i];
    os << std::setw(3) << i << "  "
       << std::showpos << tObs[i]/DAY << " "
       << std::setprecision(3)
       << std::setw(7) << dx[i]/ARCSEC << " "
       << std::setw(7) << dy[i]/ARCSEC << " "
       << std::noshowpos << std::setprecision(2) << chi
       << endl;
    chitot += chi;
  }
  os << " chisq w/o priors: " << chitot
     << " w/priors: " << chisq 
     << " DOF: " << getDOF()
     << endl;
}

Elements
Fitter::getElements() const {
  Vector3 x0;
  Vector3 v0;
  // Get state in our Frame:
  abg.getState(0.,x0,v0);
  // Convert to ICRS
  State s;
  s.x = astrometry::CartesianICRS(f.toICRS(x0));
  s.v = astrometry::CartesianICRS(f.toICRS(v0,true));
  s.tdb = f.tdb0;
  return orbits::getElements(s);
}

ElementCovariance
Fitter::getElementCovariance() const {
  // Derivative of state vector in reference frame:
  Matrix66 dSdABG_frame = abg.getStateDerivatives();
  //**/cerr << "dSdABG_frame: " << endl << dSdABG_frame << endl;

  // Rotate x and v derivatives into ICRS
  Matrix66 dSdABG;
  // Note that Frame is working with Nx3 arrays so we need
  // to put xyz on the 2nd index.
  DMatrix tmp = dSdABG_frame.subMatrix(0,3,0,6).transpose();
  dSdABG.subMatrix(0,3,0,6) = f.toICRS(tmp , true).transpose();
  tmp = dSdABG_frame.subMatrix(3,6,0,6).transpose();
  dSdABG.subMatrix(3,6,0,6) = f.toICRS(tmp, true).transpose();

  //**/cerr << "dSdABG: " << endl << dSdABG << endl;
  
  // Get state in our Frame:
  Vector3 x0;
  Vector3 v0;
  abg.getState(0.,x0,v0);
  // Convert to ICRS
  State s;
  s.x = astrometry::CartesianICRS(f.toICRS(x0));
  s.v = astrometry::CartesianICRS(f.toICRS(v0,true));
  s.tdb = f.tdb0;
  // Get element derivatives
  Matrix66 dEdABG = getElementDerivatives(s) * dSdABG;
  //**/cerr << "dEdABG: " << endl << dEdABG << endl;
  return Matrix66(dEdABG * A.inverse() * dEdABG.transpose());
}
  
// Forecast position using current fit.  Cov matrix elements given if filled:
void
Fitter::predict(const DVector& t_obs,    // Time of observations, relative to tdb0
		const DMatrix& earth,  // Observation coordinates, in our frame, Nx3
		DVector* xOut,         // Angular coordinates, in our frame
		DVector* yOut,
		DVector* covXX,        // Covar matrix of coordinates
		DVector* covYY,
		DVector* covXY) const {

  // Resize arrays if necessary
  int nobs = t_obs.rows();
  if (earth.cols()!=3 || earth.rows()!=nobs)
    throw std::runtime_error("Wrong dimensions for earth positions in Fitter::predict");
  if (!xOut || !yOut)
    throw std::runtime_error("Must provide output xOut and yOut arrays for Fitter::predict");

  bool doDerivs = covXX || covXY || covYY;
  if (doDerivs && (!covXX || !covXY || !covYY))
    throw std::runtime_error("Fitter::predict did not get arrays for all three covariance elements");

  // Iterate 3d position determination and time delay.
  // Use inertial orbit if there is no trajectory.
  DMatrix target;
  DMatrix velocity(nobs,3);
  if (fullTrajectory) {
    // Derive positions from the trajectory, place into our frame
    // Note the Trajectory takes full TDB, not referred to our tdb0:
    DVector tEphem = t_obs.array() + f.tdb0;
    DMatrix v_ICRS;
    target = f.fromICRS(fullTrajectory->position(tEphem,&v_ICRS));
    // Transform body velocity too
    velocity = f.fromICRS(v_ICRS,true);
  } else {
    target = abg.getXYZ(t_obs);
    Vector3 x,v;
    // For an inertial orbit, the velocity is constant:
    abg.getState(0., x, v);
    velocity.rowwise() = v.transpose();
  }
  // Light-time correction.  This approximation holds for the
  // case when acceleration during the light travel can be ignored:
  // light-travel time = d / (c + v_los)
  // where d is distance at time of observation and v_los is
  // line-of-sight velocity.  

  DMatrix dx = target - earth;
  DVector distance = dx.cwiseProduct(dx).rowwise().sum().cwiseSqrt();
  DVector vlos = (velocity.array()*dx.array()).rowwise().sum() / distance.array();
  DVector t_emit = t_obs.array() - distance.array() / (vlos.array() + SpeedOfLightAU);

  // Final calculation of 3d positions.  Save gravity contribution
  DMatrix gravity(nobs,3,0.);
  if (fullTrajectory) {
    DVector tEphem = t_emit.array() + f.tdb0;
    target = f.fromICRS(fullTrajectory->position(tEphem));
    // Subtract inertial motion to get gravity
    gravity = target - abg.getXYZ(t_emit);
  } else {
    target = abg.getXYZ(t_emit);
    // Gravity contribution remains zero in this case.
  }

    // Now calculate angular positions
  DVector denom = DVector(nobs,1.) + abg[ABG::GDOT]*t_emit
    + abg[ABG::G]*(gravity.col(2) - earth.col(2));
  denom = denom.cwiseInverse();  // Denom is now 1/(z*gamma)
  *xOut = abg[ABG::ADOT]*t_emit + abg[ABG::G]*(gravity.col(0) - earth.col(0));
  xOut->array() += abg[ABG::A];
  *yOut = abg[ABG::BDOT]*t_emit + abg[ABG::G]*(gravity.col(1) - earth.col(1));
  yOut->array() += abg[ABG::B];
  xOut->array() *= denom.array();
  yOut->array() *= denom.array();

  if (doDerivs) {
    // Calculate derivatives wrt ABG, using just inertial part
    DMatrix dX(nobs, 6, 0.);
    DMatrix dY(nobs, 6, 0.);
    dX.col(ABG::A).setOnes();
    dY.col(ABG::B).setOnes();
    dX.col(ABG::ADOT) = t_emit;
    dY.col(ABG::BDOT) = t_emit;
    dX.col(ABG::G) = gravity.col(0) - earth.col(0)
      + xOut->cwiseProduct(earth.col(2)-gravity.col(2));
    dY.col(ABG::G) = gravity.col(1) - earth.col(1)
      + yOut->cwiseProduct(earth.col(2)-gravity.col(2));
    dX.col(ABG::GDOT) -= xOut->cwiseProduct(t_emit);
    dY.col(ABG::GDOT) -= yOut->cwiseProduct(t_emit);

    // Now contract with the ABG covariance to get position covariance
    Matrix66 cov = A.inverse();
    DMatrix xTmp = dX * cov;
    DMatrix yTmp = dY * cov;
    *covXX = (dX.array()*xTmp.array()).rowwise().sum();
    *covXY = (dX.array()*yTmp.array()).rowwise().sum();
    *covYY = (dY.array()*yTmp.array()).rowwise().sum();

    // We left off a factor of denom in dX, dY;
    // It's faster to put these in at the end:
    denom.array() *= denom.array();  // Square this
    covXX->array() *= denom.array();
    covXY->array() *= denom.array();
    covYY->array() *= denom.array();
  }  
  return;
}

Fitter*
Fitter::augmentObservation(double tObs_,    // TDB since reference time
			   double thetaX_,  // Observed positions
			   double thetaY_,
			   double covXX_,   // Covariance of observed posns
			   double covYY_,
			   double covXY_,
			   const Vector3& xE_,     // Observatory posn at observations
			   bool newGravity) const {

  // Create new Fitter with room for one more data point
  auto out = new Fitter(eph, grav);
  out->setFrame(f);
  int oldN = tObs.size();
  out->resizeArrays(oldN + 1);
  out->tObs.subVector(0,oldN) = tObs;
  out->tObs[oldN] = tObs_;
  out->thetaX.subVector(0,oldN) = thetaX;
  out->thetaX[oldN] = thetaX_;
  out->thetaY.subVector(0,oldN) = thetaY;
  out->thetaY[oldN] = thetaY_;
  out->invcovXX.subVector(0,oldN) = invcovXX;
  out->invcovYY.subVector(0,oldN) = invcovYY;
  out->invcovXY.subVector(0,oldN) = invcovXY;
  double det = covXX_*covYY_ - covXY_*covXY_;
  out->invcovXX[oldN] = covYY_/det;
  out->invcovYY[oldN] = covXX_/det;
  out->invcovXY[oldN] = -covXY_/det;
  out->xE.subMatrix(0,oldN,0,3) = xE;
  out->xE.row(oldN) = xE_.transpose();
  
  out->xGrav.subMatrix(0,oldN,0,3) = xGrav;
  out->tEmit.subVector(0,oldN) = tEmit;
  out->tdbEmit.subVector(0,oldN) = tdbEmit;
  // Get gravity position and time delay for new point from current fit
  {
    if (!fullTrajectory)
      throw std::runtime_error("ERROR: Fitter::augmentObservation called without "
			       "a valid Trajectory");
    Vector3 xFull = f.fromICRS(fullTrajectory->position(tObs_+f.tdb0).getVector());
    Vector3 xInertial,v;
    abg.getState(tObs_,xInertial,v);
    Vector3 g = xFull-xInertial;
    out->xGrav.row(oldN) = g.transpose();

    Vector3 dx = (xFull-xE_);
    out->tEmit[oldN] = tObs_ - sqrt(dx.dot(dx))/SpeedOfLightAU; 
    out->tdbEmit[oldN] = out->tEmit[oldN] + f.tdb0;
  }

  out->newData(); // Invalidate state of new Fitter
  
  // Do Newton iterations, starting from old solution
  out->abg = abg;
  // And including same priors
  out->gamma0 = gamma0;
  out->gammaPriorSigma = gammaPriorSigma;
  out->bindingConstraintFactor = bindingConstraintFactor;
  
  // Copy old trajectory if not making new one
  if (!newGravity) out->fullTrajectory = new Trajectory(*fullTrajectory);
  out->iterateTimeDelay();
  if (newGravity) out->calculateGravity();
  out->calculateChisq(true);

  // Now run the new fitter through Newton iterations, skipping
  // gravity recalculation if newGravity==false.
  
  int iterations = 0;
  double oldChisq;
  const int MAX_ITERATIONS=5;
  const int CHISQ_TOLERANCE=0.1;
  do {
    oldChisq = out->chisq;
    auto llt = out->A.llt();
    out->abg += llt.solve(out->b);
    // Update with new ABG
    out->iterateTimeDelay();
    if (newGravity) out->calculateGravity();
    out->calculateChisq(true);
    iterations++;
  } while (iterations < MAX_ITERATIONS && abs(chisq-oldChisq) > CHISQ_TOLERANCE);
  // Do not fuss over failure to converge.
  out->abgIsFit = true;

  return out;
}  

  
