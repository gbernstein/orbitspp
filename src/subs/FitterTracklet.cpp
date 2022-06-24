// Pieces of fitting routines

//#include "Fitter.h"
#include "FitterTracklet.h"

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
  
FitterTracklet::FitterTracklet(const Ephemeris& eph_,
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
FitterTracklet::setABG(const ABG& abg_, const ABGCovariance& cov_) {
  if (!frameIsSet)
    throw std::runtime_error("ERROR: FitterTracklet::setABG called before setFrame");
  abg = abg_;
  A = cov_.inverse();
  newABG(); // update status flags
  // Set up Trajectory if using gravity
  if (grav != Gravity::INERTIAL) createTrajectory();
}

void
FitterTracklet::addTracklet(const Tracklet& track) {
  if (frameIsSet)
    throw std::runtime_error("ERROR: FitterTracklet::addObservations called after setFrame");
  tracks.push_back(track);
}



void
FitterTracklet::resizeArrays(int n) {
  tObs1.resize(n);
  tEmit1.resize(n);
  tdbEmit1.resize(n);
  thetaX1.resize(n);
  thetaY1.resize(n);
  invCovArr.resize(10,n);
  xE1.resize(n,3);
  xGrav1.resize(n,3);
  dThetaX1dABG.resize(n,6);
  dThetaY1dABG.resize(n,6);
  tObs2.resize(n);
  tEmit2.resize(n);
  tdbEmit2.resize(n);
  thetaX2.resize(n);
  thetaY2.resize(n);
  xE2.resize(n,3);
  xGrav2.resize(n,3);
  dThetaX2dABG.resize(n,6);
  dThetaY2dABG.resize(n,6);

}

void
FitterTracklet::setFrame(const Frame& f_) {
  if (frameIsSet)
    throw std::runtime_error("ERROR: FitterTracklet::setFrame called after already set");
  f = f_;
  frameIsSet = true;
  newData(); // Reset results flags

  // Recalculate all coordinates in this frame;
  int n = tracks.size();
  resizeArrays(n);

  astrometry::Gnomonic projection(f.orient);  
  for (int i=0; i<n; i++) {

    const Tracklet& tr = tracks[i];

    // Set up time variables.  Set light travel time to zero
    tObs1[i] = tr.tdb1 - f.tdb0;
    tObs2[i] = tr.tdb2 - f.tdb0;
    tEmit1[i] = tObs1[i];
    tEmit2[i] = tObs2[i];

    tdbEmit1[i] = tr.tdb1;
    tdbEmit2[i] = tr.tdb2;


    // Convert angles to our frame and get local partials for covariance
    Matrix22 partials1(0.);
    projection.convertFrom(tr.radec1, partials1);
    double x1,y1;
    projection.getLonLat(x1,y1);
    thetaX1[i] = x1;
    thetaY1[i] = y1;
    // Alter the partials for the cos(dec) factor in derivatives
    tr.radec1.getLonLat(x1,y1);
    partials1(0,0) /= cos(y1);
    partials1(1,0) /= cos(y1);

    Matrix22 partials2(0.);
    projection.convertFrom(tr.radec2, partials2);
    double x2,y2;
    projection.getLonLat(x2,y2);
    thetaX2[i] = x2;
    thetaY2[i] = y2;
    // Alter the partials for the cos(dec) factor in derivatives
    tr.radec2.getLonLat(x2,y2);
    partials2(0,0) /= cos(y2);
    partials2(1,0) /= cos(y2);

    Matrix44 fullpartials(0.);
    fullpartials(0,0) = partials1(0,0);
    fullpartials(1,1) = partials1(1,1);
    fullpartials(0,1) = fullpartials(1,0) = partials1(0,1);
    fullpartials(2,2) = partials2(0,0);
    fullpartials(3,3) = partials2(1,1);
    fullpartials(2,3) = fullpartials(3,2) = partials2(0,1);
    fullpartials(0,2) = fullpartials(2,0) = 0;
    fullpartials(0,3) = fullpartials(3,0) = 0;
    fullpartials(1,2) = fullpartials(2,1) = 0;
    fullpartials(1,3) = fullpartials(3,1) = 0;



    // Get inverse covariance
    Matrix44 cov = fullpartials * tr.cov * fullpartials.transpose();

    Matrix44 invcov = cov.inverse();
    invCovArr(0,i) = invcov(0,0);
    invCovArr(1,i) = invcov(1,1);
    invCovArr(2,i) = invcov(2,2);

    invCovArr(3,i) = invcov(3,3);
    invCovArr(4,i) = invcov(0,1);
    invCovArr(5,i) = invcov(0,2);
    invCovArr(6,i) = invcov(0,3);
    invCovArr(7,i) = invcov(1,2);

    invCovArr(8,i) = invcov(1,3);
    invCovArr(9,i) = invcov(2,3);


    // Get observatory position in our frame.
    xE1.row(i) = f.fromICRS(tr.observer1.getVector());
    xE2.row(i) = f.fromICRS(tr.observer2.getVector());


    // Initialize gravitational effect to zero
    xGrav1.setZero();
    xGrav2.setZero();


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
FitterTracklet::setObservationsInFrame(const DVector& tObs1_,    // TDB since reference time
             const DVector& tObs2_,
			       const DVector& thetaX1_,  // Observed positions
			       const DVector& thetaY1_,
             const DVector& thetaX2_,  // Observed positions
             const DVector& thetaY2_,
			       const DMatrix& covFull_,   // Covariance of observed posns
			       const DMatrix& xE1_,
             const DMatrix& xE2_) {    // Observatory posn at observations
  if (!frameIsSet)
    throw std::runtime_error("ERROR: FitterTracklet::setObservationsInFrame() was "
			     "called before setFrame()");
  int n = tObs1_.size();
  if (thetaX1_.size()!=n ||
      thetaY1_.size()!=n ||
      xE1_.rows() !=n ||
      xE1_.cols() !=3)
    throw std::runtime_error("ERROR: Mismatched array sizes in "
			     "FitterTracklet::setObservationsInFrame()");

  resizeArrays(n);
  newData(); // Reset results flags

  tObs1 = tObs1_;
  tEmit1 = tObs1_; // No light-travel correction to start
  tdbEmit1 = tEmit1.array() + f.tdb0;
  tObs2 = tObs2_;
  tEmit2 = tObs2_; // No light-travel correction to start
  tdbEmit2 = tEmit2.array() + f.tdb0;

  thetaX1 = thetaX1_;
  thetaY1 = thetaY1_;

  thetaX2 = thetaX2_;
  thetaY2 = thetaY2_;

  Matrix44 covar; 
  Matrix44 invcov;
  for (int i; i < n; i++){
    covar(0,0) = covFull_(0,i);
    covar(1,1) = covFull_(1,i);
    covar(2,2) = covFull_(2,i);
    covar(3,3) = covFull_(3,i);
    covar(0,1) = covar(1,0) = covFull_(4,i);
    covar(0,2) = covar(2,0) = covFull_(5,i);
    covar(0,3) = covar(3,0) = covFull_(6,i);
    covar(1,2) = covar(2,1) = covFull_(7,i);
    covar(1,3) = covar(3,1) = covFull_(8,i);
    covar(2,3) = covar(3,2) = covFull_(9,i);
    invcov = covar.inverse();

    invCovArr(0,i) = invcov(0,0);
    invCovArr(1,i) = invcov(1,1);
    invCovArr(2,i) = invcov(2,2);
    invCovArr(3,i) = invcov(3,3);
    invCovArr(4,i) = invcov(0,1);
    invCovArr(5,i) = invcov(0,2);
    invCovArr(6,i) = invcov(0,3);
    invCovArr(7,i) = invcov(1,2);
    invCovArr(8,i) = invcov(1,3);
    invCovArr(9,i) = invcov(2,3);

  }


  xE1 = xE1_;
  xGrav1.setZero();
  xE2 = xE2_;
  xGrav2.setZero();

}

void
FitterTracklet::chooseFrame(int obsNumber) {
  if (obsNumber >= static_cast<int> (tracks.size())) 
    throw runtime_error("Request for nonexistent obsNumber in FitterTracklet::chooseFrame()");
  if (tracks.empty())
    throw runtime_error("Must have tracks to call FitterTracklet::chooseFrame()");
  if (obsNumber < 0) {
    double tdb0 = tracks.front().tdb1;
    double tsum=0;
    for (auto& tr : tracks)
      tsum += tr.tdb1 - tdb0;
    double target = tdb0 + tsum / tracks.size();
    double minDT = 1e9;
    for (int i=0; i<tracks.size(); i++) {
      if (abs(tracks[i].tdb1-target)<minDT) {
	obsNumber = i;
	minDT = abs(tracks[i].tdb1-target);
      }
    }
  }
  
  // Now construct ecliptic-aligned frame and set up all observations
  const Tracklet& tr = tracks[obsNumber];
  astrometry::Orientation orient(tr.radec1);
  orient.alignToEcliptic();
  double tdb0 = tr.tdb1;
  astrometry::CartesianICRS origin = tr.observer1;
  setFrame(Frame(origin, orient, tdb0));
}

void
FitterTracklet::calculateOrbit(bool doDerivatives) {
  int nobs = tEmit1.size();



  // Is is already done?
  if (doDerivatives && positionsAreValid && positionDerivsAreValid) return;
  if (!doDerivatives && positionsAreValid) return;
  
  // Calculate positions
  DVector denom1 = DVector(nobs,1.) + abg[ABG::GDOT]*tEmit1
    + abg[ABG::G]*(xGrav1.col(2) - xE1.col(2));

  DVector denom2 = DVector(nobs,1.) + abg[ABG::GDOT]*tEmit2
    + abg[ABG::G]*(xGrav2.col(2) - xE2.col(2));

  denom1 = denom1.cwiseInverse();  // Denom is now 1/(z*gamma)
  thetaXModel1 = abg[ABG::ADOT]*tEmit1
    + abg[ABG::G]*(xGrav1.col(0) - xE1.col(0));
  thetaXModel1.array() += abg[ABG::A];
  thetaYModel1 = abg[ABG::BDOT]*tEmit1
    + abg[ABG::G]*(xGrav1.col(1) - xE1.col(1));
  thetaYModel1.array() += abg[ABG::B];
  thetaXModel1.array() *= denom1.array();
  thetaYModel1.array() *= denom1.array();


  denom2 = denom2.cwiseInverse();  // Denom is now 1/(z*gamma)
  thetaXModel2 = abg[ABG::ADOT]*tEmit2
    + abg[ABG::G]*(xGrav2.col(0) - xE2.col(0));
  thetaXModel2.array() += abg[ABG::A];
  thetaYModel2 = abg[ABG::BDOT]*tEmit2
    + abg[ABG::G]*(xGrav2.col(1) - xE2.col(1));
  thetaYModel2.array() += abg[ABG::B];
  thetaXModel2.array() *= denom2.array();
  thetaYModel2.array() *= denom2.array();

  positionsAreValid = true;

  if (doDerivatives) {
    // Calculate derivatives
    dThetaX1dABG.setZero();
    dThetaY1dABG.setZero();
    dThetaX1dABG.col(ABG::A).setOnes();
    dThetaY1dABG.col(ABG::B).setOnes();
    dThetaX1dABG.col(ABG::ADOT) = tEmit1;
    dThetaY1dABG.col(ABG::BDOT) = tEmit1;
    dThetaX1dABG.col(ABG::G) = xGrav1.col(0) - xE1.col(0)
      + thetaXModel1.cwiseProduct(xE1.col(2)-xGrav1.col(2));
    dThetaY1dABG.col(ABG::G) = xGrav1.col(1) - xE1.col(1)
      + thetaYModel1.cwiseProduct(xE1.col(2)-xGrav1.col(2));
    dThetaX1dABG.col(ABG::GDOT) -= thetaXModel1.cwiseProduct(tEmit1);
    dThetaY1dABG.col(ABG::GDOT) -= thetaYModel1.cwiseProduct(tEmit1);

    dThetaX1dABG.applyOnTheLeft(denom1.asDiagonal());
    dThetaY1dABG.applyOnTheLeft(denom1.asDiagonal());

    dThetaX2dABG.setZero();
    dThetaY2dABG.setZero();
    dThetaX2dABG.col(ABG::A).setOnes();
    dThetaY2dABG.col(ABG::B).setOnes();
    dThetaX2dABG.col(ABG::ADOT) = tEmit2;
    dThetaY2dABG.col(ABG::BDOT) = tEmit2;
    dThetaX2dABG.col(ABG::G) = xGrav2.col(0) - xE2.col(0)
      + thetaXModel2.cwiseProduct(xE2.col(2)-xGrav2.col(2));
    dThetaY2dABG.col(ABG::G) = xGrav2.col(1) - xE2.col(1)
      + thetaYModel2.cwiseProduct(xE2.col(2)-xGrav2.col(2));
    dThetaX2dABG.col(ABG::GDOT) -= thetaXModel2.cwiseProduct(tEmit2);
    dThetaY2dABG.col(ABG::GDOT) -= thetaYModel2.cwiseProduct(tEmit2);

    dThetaX2dABG.applyOnTheLeft(denom2.asDiagonal());
    dThetaY2dABG.applyOnTheLeft(denom2.asDiagonal());

    positionDerivsAreValid = true;
  }
  
  // Downstream quantities are invalid now
  chisqIsValid = chisqDerivsAreValid = false; 
  return;
}

void
FitterTracklet::iterateTimeDelay() {
  // Use current orbit to update light-travel time
  int nobs = tObs1.size();
  // Light-travel time is distance/(speed of light)

  DMatrix xyz1 = abg.getXYZ(tEmit1) + xGrav1 - xE1;
  // Get sqrt(xyz.xyz):
  DVector dist1 = xyz1.cwiseProduct(xyz1).rowwise().sum().cwiseSqrt();
  if ( (dist1.array()>1e4).any())
    throw NonConvergent("distance >10,000 AU in FitterTracklet::iterateTimeDelay()");
  tEmit1 = tObs1 - dist1/SpeedOfLightAU;
  tdbEmit1.array() = tEmit1.array() + f.tdb0;

  DMatrix xyz2 = abg.getXYZ(tEmit2) + xGrav2 - xE2;
  // Get sqrt(xyz.xyz):
  DVector dist2 = xyz2.cwiseProduct(xyz2).rowwise().sum().cwiseSqrt();
  if ( (dist2.array()>1e4).any())
    throw NonConvergent("distance >10,000 AU in FitterTracklet::iterateTimeDelay()");
  tEmit2 = tObs2 - dist2/SpeedOfLightAU;
  tdbEmit2.array() = tEmit2.array() + f.tdb0;


  // This invalidates the solution
  // Solutions no longer match gravity
  newABG();
  abgIsFit = false;
  
  return;
}

void
FitterTracklet::createTrajectory() {
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
FitterTracklet::calculateGravity() {
  // Calculate the non-inertial terms in trajectory of object
  // using current ABG and gravity settings
  if (grav==Gravity::INERTIAL) {
    xGrav1.setZero();
    xGrav2.setZero();
    return;
  }

  createTrajectory();
  
  DMatrix xyz1 = fullTrajectory->position(tdbEmit1);
  DMatrix xyz2 = fullTrajectory->position(tdbEmit2);
  //cerr << xyz1 << endl;
  //cerr << xyz2 << endl;

  // Convert back to Fitter reference frame and subtract inertial motion
  // Note Frame and Trajectory use 3 x n, we are using n x 3 here.
  xGrav1 = f.fromICRS(xyz1) - abg.getXYZ(tEmit1);
  xGrav2 = f.fromICRS(xyz2) - abg.getXYZ(tEmit2);
  //cerr << "xgrav1" << endl;
  //cerr << xGrav1 << endl;
  //cerr << "xgrav2" << endl;

  //cerr << xGrav1 << endl;


  // Solutions no longer match gravity
  newABG();
  abgIsFit = false;
  return;
}

double
FitterTracklet::getChisq() {
  calculateChisq(false);
  return chisq;
}

void
FitterTracklet::calculateChisq(bool doDerivatives, bool doCovariances) {
  if (tObs1.size() < 2)
    throw std::runtime_error("ERROR: FitterTracklet::calculateChisq called with <2 observations");

  // Is it already done?
  //if (doDerivatives && chisqIsValid && chisqDerivsAreValid) return;
  //if (!doDerivatives && chisqIsValid) return;
  
  // Make sure that model positions are current (and derivatives, if needed)
  calculateOrbit(doDerivatives);
    
  DVector dx1 = thetaX1 - thetaXModel1;
  DVector dy1 = thetaY1 - thetaYModel1;
  DVector dx2 = thetaX2 - thetaXModel2;
  DVector dy2 = thetaY2 - thetaYModel2;

  DVector invCovX1X1 = invCovArr.row(0);
  DVector invCovY1Y1 = invCovArr.row(1);
  DVector invCovX2X2 = invCovArr.row(2);
  DVector invCovY2Y2 = invCovArr.row(3);
  DVector invCovX1Y1 = invCovArr.row(4);
  DVector invCovX1X2 = invCovArr.row(5);
  DVector invCovX1Y2 = invCovArr.row(6);
  DVector invCovY1X2 = invCovArr.row(7);
  DVector invCovY1Y2 = invCovArr.row(8);
  DVector invCovX2Y2 = invCovArr.row(9);
    

  chisq = dx1.transpose() * (invCovX1X1.asDiagonal() * dx1);
  chisq += 2.* dx1.transpose() * (invCovX1Y1.asDiagonal() * dy1);
  chisq += dy1.transpose() * (invCovY1Y1.asDiagonal() * dy1);

  chisq += dx2.transpose() * (invCovX2X2.asDiagonal() * dx2);
  chisq += 2.* dx2.transpose() * (invCovX2Y2.asDiagonal() * dy2);
  chisq += dy2.transpose() * (invCovY2Y2.asDiagonal() * dy2);

  if (doCovariances){
  chisq += 2.* dx1.transpose() * (invCovX1X2.asDiagonal() * dx2);
  chisq += 2.* dx1.transpose() * (invCovX1Y2.asDiagonal() * dy2);

  chisq += 2.* dy1.transpose() * (invCovY1X2.asDiagonal() * dx2);
  chisq += 2.* dy1.transpose() * (invCovY1Y2.asDiagonal() * dy2);
  }



  if (doDerivatives) {
    DMatrix tmp = invCovX1X1.asDiagonal() * dThetaX1dABG;
    b = tmp.transpose() * dx1;
    A = dThetaX1dABG.transpose() * tmp;

    tmp = invCovX1Y1.asDiagonal() * dThetaX1dABG;
    b += tmp.transpose() * dy1;
    A += dThetaY1dABG.transpose() * tmp;

    tmp = invCovX1Y1.asDiagonal() * dThetaY1dABG;
    b += tmp.transpose() * dx1;
    A += dThetaX1dABG.transpose() * tmp;

    tmp = invCovY1Y1.asDiagonal() * dThetaY1dABG;
    b += tmp.transpose() * dy1;
    A += dThetaY1dABG.transpose() * tmp;

    tmp = invCovX2X2.asDiagonal() * dThetaX2dABG;
    b += tmp.transpose() * dx2;
    A += dThetaX2dABG.transpose() * tmp;

    tmp = invCovX2Y2.asDiagonal() * dThetaX2dABG;
    b += tmp.transpose() * dy2;
    A += dThetaY2dABG.transpose() * tmp;

    tmp = invCovX2Y2.asDiagonal() * dThetaY2dABG;
    b += tmp.transpose() * dx2;
    A += dThetaX2dABG.transpose() * tmp;

    tmp = invCovY2Y2.asDiagonal() * dThetaY2dABG;
    b += tmp.transpose() * dy2;
    A += dThetaY2dABG.transpose() * tmp;

    if (doCovariances) {
    tmp = invCovX1X2.asDiagonal() * dThetaX1dABG;
    b += tmp.transpose() * dx2;
    A += dThetaX2dABG.transpose() * tmp;

    tmp = invCovX1X2.asDiagonal() * dThetaX2dABG;
    b += tmp.transpose() * dx1;
    A += dThetaX1dABG.transpose() * tmp;

    tmp = invCovX1Y2.asDiagonal() * dThetaX1dABG;
    b += tmp.transpose() * dy2;
    A += dThetaY2dABG.transpose() * tmp;

    tmp = invCovX1Y2.asDiagonal() * dThetaY2dABG;
    b += tmp.transpose() * dx1;
    A += dThetaX1dABG.transpose() * tmp;

    tmp = invCovY1X2.asDiagonal() * dThetaY1dABG;
    b += tmp.transpose() * dx2;
    A += dThetaX2dABG.transpose() * tmp;

    tmp = invCovY1X2.asDiagonal() * dThetaX2dABG;
    b += tmp.transpose() * dy1;
    A += dThetaY1dABG.transpose() * tmp;

    tmp = invCovY1Y2.asDiagonal() * dThetaY1dABG;
    b += tmp.transpose() * dy2;
    A += dThetaY2dABG.transpose() * tmp;

    tmp = invCovY1Y2.asDiagonal() * dThetaY2dABG;
    b += tmp.transpose() * dy1;
    A += dThetaY1dABG.transpose() * tmp;

  }


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


void FitterTracklet::setSingleOrbit() {
  const int MAX_ITERATIONS=50;
  double lambda = 0.01;
  double olda, oldb, oldg, oldadot, oldbdot, oldgdot;
  int iterations = 0;
  double oldChisq = 0;
  DVector invCovX1X1 = invCovArr.row(0);
  DVector invCovY1Y1 = invCovArr.row(1);
  DVector invCovX2X2 = invCovArr.row(2);
  DVector invCovY2Y2 = invCovArr.row(3);
  DVector invCovMean1 = (invCovY1Y1 + invCovX1X1)/2;
  DVector invCovMean2 = (invCovY2Y2 + invCovX2X2)/2;

  do{
    olda = abg[ABG::A];
    oldb = abg[ABG::B];
    oldg = abg[ABG::G];
    oldadot = abg[ABG::ADOT];
    oldbdot = abg[ABG::BDOT];
    oldgdot = abg[ABG::GDOT];

  calculateOrbit(true);


  DVector dx1 = thetaX1 - thetaXModel1;
  DVector dy1 = thetaY1 - thetaYModel1;
  DVector dx2 = thetaX1 - thetaXModel1;
  DVector dy2 = thetaY1 - thetaYModel1;



  chisq = dx1.transpose() * (invCovMean1.asDiagonal() * dx1);
  chisq += dy1.transpose() * (invCovMean1.asDiagonal() * dy1);
  chisq += dx2.transpose() * (invCovMean2.asDiagonal() * dx2);
  chisq += dy2.transpose() * (invCovMean2.asDiagonal() * dy2);

    DMatrix tmp = invCovMean1.asDiagonal() * dThetaX1dABG;
    b = tmp.transpose() * dx1 ;
    A = dThetaX1dABG.transpose() * tmp;

    tmp = invCovMean1.asDiagonal() * dThetaY1dABG;
    b += tmp.transpose() * dy1;
    A += dThetaY1dABG.transpose() * tmp;

    tmp = invCovMean2.asDiagonal() * dThetaY2dABG;
    b += tmp.transpose() * dy2;
    A += dThetaY2dABG.transpose() * tmp;

    tmp = invCovMean2.asDiagonal()  * dThetaX2dABG;
    b += tmp.transpose() * dx2;
    A += dThetaX2dABG.transpose() * tmp;



  if (gammaPriorSigma > 0.) {
    double w = pow(gammaPriorSigma, -2);
    double dg = (gamma0-abg[ABG::G]);
    chisq += dg * w * dg;
    
      b[ABG::G] += w * dg;
      A(ABG::G,ABG::G) += w;
    
  }
  // Extract parameters other than GDOT, which we assume is last
  DVector blin = b.subVector(0, ABG::GDOT);
  DMatrix Alin = A.subMatrix(0, ABG::GDOT, 0, ABG::GDOT);

    
    Alin(ABG::A, ABG::A) *= (1.+lambda);
    Alin(ABG::B, ABG::B) *= (1.+lambda);
    Alin(ABG::G, ABG::G) *= (1.+lambda);
    Alin(ABG::ADOT, ABG::ADOT) *= (1.+lambda);
    Alin(ABG::BDOT, ABG::BDOT) *= (1.+lambda);
    //A(ABG::GDOT, ABG::GDOT) *= (1+lambda);
  
  // Solve (check degeneracies??)
  auto llt = Alin.llt();
  DVector answer = llt.solve(blin);


  // Transfer answer to instance members
  abg.subVector(0,ABG::GDOT) += answer;
  //if (abg[ABG::G] < 0) abg[ABG::G] *= -1;
  abg[ABG::GDOT] = 0.;

  oldChisq = chisq;

  calculateChisq(false);
  /*
  if (chisq > oldChisq && iterations > 0){
      lambda *= 10;

    abg[ABG::A] = olda;
    abg[ABG::B] = oldb;
    abg[ABG::G] = oldg;
    abg[ABG::ADOT] = oldadot;
    abg[ABG::BDOT] = oldbdot;
    abg[ABG::GDOT] = oldgdot;
      //cerr << "Updated" << endl;

      //cerr << abg << endl;

    }
    else if (iterations > 0)
    {
      lambda /= 10;
    }
  */
    positionsAreValid = false;
    positionDerivsAreValid = false;
    //newABG();
    iterations++;
    //cerr << chisq << " " << oldChisq << endl;
  } while (iterations<MAX_ITERATIONS && abs(chisq-oldChisq) > 0.01);

  newABG(); // Reset results flags
  // But don't consider linear solution to be valid
  abgIsFit = false;

  // Update positions, chisq - no derivatives
  calculateChisq(false, true);

}

void
FitterTracklet::setLinearOrbit() {

  // Get positions/derivatives for current ABG
  calculateChisq(true, false);

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
  calculateChisq(false, false);
}

void
FitterTracklet::abgSanityCheck() const {
  const double MAX_GAMMA = 1.; // Any closer than this is failure
  const double MAX_KE = 10.; // Failure when crude |KE/PE| exceeds this

  if (abs(abg[ABG::G]) > MAX_GAMMA) {
    FormatAndThrow<NonConvergent>()
      << "gamma grows too large: " << abg[ABG::G];
  }
  double ke = (abg[ABG::ADOT] * abg[ABG::ADOT]
	       + abg[ABG::BDOT] * abg[ABG::BDOT]
	       + abg[ABG::GDOT] * abg[ABG::GDOT]) / (2. * GM * pow(abs(abg[ABG::G]),3.));
  if (ke > MAX_KE) {
    FormatAndThrow<NonConvergent>()
      << "|KE/PE| grows too large: " << ke;
  }
}

void
FitterTracklet::newtonFit(double chisqTolerance, bool dump, bool doCovariances) {
  // Assuming that we are coming into this with a good starting point
  iterateTimeDelay();
  calculateGravity();
  calculateChisq(true, doCovariances);
  const int MAX_ITERATIONS = 1000000;

  double lambda = 0.01;
  ABG oldabg;

  // cancel validity of results until re-converged
  abgIsFit = false;
  
  // Do a series of iterations with time delay and gravity fixed
  int iterations = 0;
  double oldChisq, prevchisq;
  do {
    prevchisq = chisq;
    calculateChisq(true, doCovariances);
    // Solve (check degeneracies??)
    // LM change
    A(ABG::A, ABG::A) *= (1.+lambda);
    A(ABG::B, ABG::B) *= (1.+lambda);
    A(ABG::G, ABG::G) *= (1.+lambda);
    A(ABG::ADOT, ABG::ADOT) *= (1.+lambda);
    A(ABG::BDOT, ABG::BDOT) *= (1.+lambda);
    A(ABG::GDOT, ABG::GDOT) *= (1.+lambda);
     
    auto llt = A.llt();
    oldabg = abg;
    abg += llt.solve(b);

    if (dump) {
      cerr << "Iteration " << iterations <<endl;
      auto dd = llt.solve(b);
      cerr << " change: ";
      write6(dd,cerr);
      cerr << endl << " abg:    ";
      write6(abg,cerr);
      cerr << endl << "old:     ";
      write6(oldabg,cerr);
    }

    // Quit if iterations have gone completely awry
    //abgSanityCheck();

    // Recalculate with new ABG
    oldChisq = chisq;
    calculateChisq(false, doCovariances);
    if (chisq >= oldChisq){
      cerr << "changing lambda " << lambda << endl;
      lambda *= 10.;
      abg = oldabg;
      //calculateChisq(true);

    }
    else{
      lambda /= 10.;
    }
    
    cerr << endl << " New chisq: " << chisq << endl;
    cerr << "Delta chisq " << abs(chisq - prevchisq) <<  endl;
    iterations++;
    positionsAreValid = false;
    positionDerivsAreValid = false;

  } while (iterations < MAX_ITERATIONS && abs(chisq-prevchisq) > chisqTolerance);
  /**
#pragma omp critical (io)
  cerr << "--NewtonFit done at iteration " << iterations << endl; 
  ***/
  //if (iterations >= MAX_ITERATIONS)
  //  throw NonConvergent("newtonFit exceeded max iterations");
   
  // Then converge with time delay and gravity recalculated
  iterateTimeDelay();
  calculateGravity();
  calculateChisq(true, doCovariances);
  iterations = 0;
  lambda = 0.01;
  do {
    // Solve (check degeneracies??)

    oldabg = abg;
    A(ABG::A, ABG::A) *= (1.+lambda);
    A(ABG::B, ABG::B) *= (1.+lambda);
    A(ABG::G, ABG::G) *= (1.+lambda);
    A(ABG::ADOT, ABG::ADOT) *= (1.+lambda);
    A(ABG::BDOT, ABG::BDOT) *= (1.+lambda);
    A(ABG::GDOT, ABG::GDOT) *= (1.+lambda);
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
    //abgSanityCheck();

    // Update with new ABG
    iterateTimeDelay();
    calculateGravity();
    positionsAreValid = false;
    oldChisq = chisq;
    calculateChisq(false);

    
    if (chisq > oldChisq){
      lambda *= 10.;
      //cerr << "changing lambda " << lambda << endl;

      abg = oldabg;
      //calculateChisq(true);

    }
    else{
      lambda /= 10.;
    }
    
    if (dump) cerr << endl << " New chisq: " << chisq << endl;
    iterations++;
  } while (iterations < MAX_ITERATIONS && abs(chisq-oldChisq) > chisqTolerance);
  if (iterations >= MAX_ITERATIONS)
    throw NonConvergent("newtonFit exceeded max fine iterations");
  abgIsFit = true;
  return;
}

void
FitterTracklet::printResiduals(std::ostream& os) {
  // Update residuals and chisq to current ABG if needed
  calculateChisq(false);

  stringstuff::StreamSaver ss(os);
  os << "# Residuals: " << endl << std::fixed << std::setprecision(2);
  double chitot = 0;
  double chi;
  DVector dx1 = thetaX1 - thetaXModel1;
  DVector dy1 = thetaY1 - thetaYModel1;
  DVector dx2 = thetaX2 - thetaXModel2;
  DVector dy2 = thetaY2 - thetaYModel2;
  
  os << "# N    T        dx      dy  chisq" << endl;
  os << "#    (days)        (mas)   " << endl;
  for (int i=0; i<dx1.size(); i++) {
    Matrix44 covar, invcov; 

      covar(0,0) = invCovArr(0,i);
      covar(1,1) = invCovArr(1,i);
      covar(2,2) = invCovArr(2,i);
      covar(3,3) = invCovArr(3,i);
      covar(0,1) = covar(1,0) = invCovArr(4,i);
      covar(0,2) = covar(2,0) = invCovArr(5,i);
      covar(0,3) = covar(3,0) = invCovArr(6,i);
      covar(1,2) = covar(2,1) = invCovArr(7,i);
      covar(1,3) = covar(3,1) = invCovArr(8,i);
      covar(2,3) = covar(3,2) = invCovArr(9,i);

      invcov = covar;

      chi = dx1[i] * dx1[i] * invcov(0,0) + dy1[i] * dy1[i] * invcov(1,1) + 2 * dx1[i] * dy1[i] * invcov(0,1);
      chi += dx2[i] * dx2[i] * invcov(2,2) + dy2[i] * dy2[i] * invcov(3,3) + 2 * dx2[i] * dy2[i] * invcov(2,3);
      chi += 2*dx1[i] * dx2[i] * invcov(0,2) + 2*dx1[i] * dy2[i] * invcov(0,3);
      chi += 2 * dy1[i] * dx2[i] * invcov(1,2) + 2 * dy1[i] * dy2[i] * invcov(1,3);


    os << std::setw(3) << i << "  "
       << std::setprecision(2)
       << std::setw(7) << std::showpos << tObs1[i]/DAY << " "
       << std::setw(7) << std::showpos << tObs2[i]/DAY << " "
       << std::setprecision(2)
       << std::setw(7) << 1000.*dx1[i]/ARCSEC << " "
       << std::setw(7) << 1000.*dy1[i]/ARCSEC << " "
       << std::setw(7) << 1000.*dx2[i]/ARCSEC << " "
       << std::setw(7) << 1000.*dy2[i]/ARCSEC << " "
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
FitterTracklet::getElements(bool heliocentric) const {
  Vector3 x0;
  Vector3 v0;
  // Get state in our Frame:
  abg.getState(0.,x0,v0);
  State sRef;
  sRef.x = x0;
  sRef.v = v0;
  sRef.tdb = f.tdb0;
  // Convert to ICRS
  State s = f.toICRS(sRef);
  return heliocentric ? orbits::getElements(s, true, &eph) : orbits::getElements(s);
}

ElementCovariance
FitterTracklet::getElementCovariance(bool heliocentric) const {
  // Derivative of state vector in reference frame:
  Matrix66 dSdABG_frame = abg.getStateDerivatives();

  // Rotate x and v derivatives into ICRS
  Matrix66 dSdABG;
  // Note that Frame is working with Nx3 arrays so we need
  // to put xyz on the 2nd index.
  DMatrix tmp = dSdABG_frame.subMatrix(0,3,0,6).transpose();
  dSdABG.subMatrix(0,3,0,6) = f.toICRS(tmp , true).transpose();
  tmp = dSdABG_frame.subMatrix(3,6,0,6).transpose();
  dSdABG.subMatrix(3,6,0,6) = f.toICRS(tmp, true).transpose();

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
  Matrix66 dEdABG = getElementDerivatives(s, heliocentric, &eph) * dSdABG;
  //**/cerr << "dEdABG: " << endl << dEdABG << endl;
  return Matrix66(dEdABG * A.inverse() * dEdABG.transpose());
}

Elements
FitterTracklet::getElements(double tdb, bool heliocentric) const {
  State s = predictState(tdb - f.tdb0);
  return heliocentric ? orbits::getElements(s, heliocentric, &eph) : orbits::getElements(s) ;
}

ElementCovariance
FitterTracklet::getElementCovariance(double tdb, bool heliocentric) const {
  // Derivative of state vector in reference frame:
  Matrix66 covS;
  State s = predictState(tdb-f.tdb0, &covS);
  // Get element derivatives
  Matrix66 dEdS = getElementDerivatives(s, heliocentric, &eph);
  return Matrix66(dEdS * covS * dEdS.transpose());
}

// Forecast position using current fit.  Cov matrix elements given if filled:
void
FitterTracklet::predict(const DVector& t_obs,    // Time of observations, relative to tdb0
		const DMatrix& earth,  // Observation coordinates, in our frame, Nx3
		DVector* xOut,         // Angular coordinates, in our frame
		DVector* yOut,
		DVector* covXX,        // Covar matrix of coordinates
		DVector* covYY,
		DVector* covXY) const {

  // Resize arrays if necessary
  int nobs = t_obs.rows();
  if (earth.cols()!=3 || earth.rows()!=nobs)
    throw std::runtime_error("Wrong dimensions for earth positions in FitterTracklet::predict");
  if (!xOut || !yOut)
    throw std::runtime_error("Must provide output xOut and yOut arrays for FitterTracklet::predict");

  bool doDerivs = covXX || covXY || covYY;
  if (doDerivs && (!covXX || !covXY || !covYY))
    throw std::runtime_error("FitterTracklet::predict did not get arrays for all three covariance elements");

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

// Forecast position using current fit.  Cov matrix elements given if filled.
// Position uses full orbit but the covariance will be propagated with inertial approx
State
FitterTracklet::predictState(const double tdb,    // Dynamical time, relative to tdb0
		     Matrix66* stateCov) const {

  State out;

  // Final calculation of 3d positions.
  if (fullTrajectory) {
    out.x = fullTrajectory->position(tdb + f.tdb0, &out.v);
  } else {
    Vector3 x,v;
    abg.getState(tdb, x, v);
    // Convert to ICRS
    out.x = f.toICRS(x);
    out.v = f.toICRS(v,true);
  }
  out.tdb = tdb + f.tdb0;

  if (stateCov) {
    // Calculate covariance matrix of state.  Propagate from frame time assuming inertial motion.

    // Get ABG covariance at reference time
    Matrix66 abgCov = A.inverse();

    // Get derivatives d(ICRS state)/d(ABG) at reference time
    // by combining d(ICRS state vector) / d(Ref state vector)  * d(ref state)/d(ABG)
    Matrix66 dICRSdABG = f.dICRSdRef() * abg.getStateDerivatives();
    // Propagate uncertainty to emission time.
    Matrix66 dNowdThen(0.);
    for (int i=0; i<6; i++)
      dNowdThen(i,i) = 1.;
    for (int i=0; i<3; i++)
      dNowdThen(i,i+3) = tdb;

    Matrix66 dNowdABG = dNowdThen * dICRSdABG;

    *stateCov = dNowdABG * abgCov * dNowdABG.transpose();
  }  
  return out ;
}

FitterTracklet*
FitterTracklet::augmentObservation(double tObs1_,    // TDB since reference time
         double tObs2_,
			   double thetaX1_,  // Observed positions
			   double thetaY1_,
         double thetaX2_,  // Observed positions
         double thetaY2_,
			   DVector& cov_,   // Covariance of observed posns
			   const Vector3& xE1_,     // Observatory posn at observations,
         const Vector3& xE2_,     // Observatory posn at observations
			   bool newGravity) const {

  // Create new Fitter with room for one more data point
  auto out = new FitterTracklet(eph, grav);
  out->setFrame(f);
  int oldN = tObs1.size();
  out->resizeArrays(oldN + 1);
  out->tObs1.subVector(0,oldN) = tObs1;
  out->tObs1[oldN] = tObs1_;
  out->tObs2.subVector(0,oldN) = tObs2;
  out->tObs2[oldN] = tObs2_;

  out->thetaX1.subVector(0,oldN) = thetaX1;
  out->thetaX1[oldN] = thetaX1_;
  out->thetaY1.subVector(0,oldN) = thetaY1;
  out->thetaY1[oldN] = thetaY1_;

  out->thetaX2.subVector(0,oldN) = thetaX2;
  out->thetaX2[oldN] = thetaX2_;
  out->thetaY2.subVector(0,oldN) = thetaY2;
  out->thetaY2[oldN] = thetaY2_;


  out->invCovArr.subMatrix(0,oldN, 0, 9) = invCovArr;

  Matrix44 covar, invcov;

  covar(0,0) = cov_(0);
  covar(1,1) = cov_(1);
  covar(2,2) = cov_(2);
  covar(3,3) = cov_(3);
  covar(0,1) = covar(1,0) = cov_(4);
  covar(0,2) = covar(2,0) = cov_(5);
  covar(0,3) = covar(3,0) = cov_(6);
  covar(1,2) = covar(2,1) = cov_(7);
  covar(1,3) = covar(3,1) = cov_(8);
  covar(2,3) = covar(3,2) = cov_(9);
  invcov = covar.inverse();



  out->invCovArr(0,oldN) = invcov(0,0);
  out->invCovArr(1,oldN) = invcov(1,1);
  out->invCovArr(2,oldN) = invcov(2,2);
  out->invCovArr(3,oldN) = invcov(3,3);
  out->invCovArr(4,oldN) = invcov(0,1);
  out->invCovArr(5,oldN) = invcov(0,2);
  out->invCovArr(6,oldN) = invcov(0,3);
  out->invCovArr(7,oldN) = invcov(1,2);
  out->invCovArr(8,oldN) = invcov(1,3);
  out->invCovArr(9,oldN) = invcov(2,3);



  out->xE1.subMatrix(0,oldN,0,3) = xE1;
  out->xE1.row(oldN) = xE1_.transpose();
  
  out->xE2.subMatrix(0,oldN,0,3) = xE2;
  out->xE2.row(oldN) = xE2_.transpose();
  
  out->xGrav1.subMatrix(0,oldN,0,3) = xGrav1;
  out->tEmit1.subVector(0,oldN) = tEmit1;
  out->tdbEmit1.subVector(0,oldN) = tdbEmit1;
  out->xGrav2.subMatrix(0,oldN,0,3) = xGrav2;
  out->tEmit2.subVector(0,oldN) = tEmit2;
  out->tdbEmit2.subVector(0,oldN) = tdbEmit2;
  
  // Get gravity position and time delay for new point from current fit
  {
    if (!fullTrajectory)
      throw std::runtime_error("ERROR: FitterTracklet::augmentObservation called without "
			       "a valid Trajectory");
    Vector3 xFull = f.fromICRS(fullTrajectory->position(tObs1_+f.tdb0).getVector());
    Vector3 xInertial,v;
    abg.getState(tObs1_,xInertial,v);
    Vector3 g = xFull-xInertial;
    out->xGrav1.row(oldN) = g.transpose();

    Vector3 dx = (xFull-xE1_);
    out->tEmit1[oldN] = tObs1_ - sqrt(dx.dot(dx))/SpeedOfLightAU; 
    out->tdbEmit1[oldN] = out->tEmit1[oldN] + f.tdb0;

    xFull = f.fromICRS(fullTrajectory->position(tObs2_+f.tdb0).getVector());
    xInertial,v;
    abg.getState(tObs2_,xInertial,v);
    g = xFull-xInertial;
    out->xGrav2.row(oldN) = g.transpose();

    dx = (xFull-xE2_);
    out->tEmit2[oldN] = tObs2_ - sqrt(dx.dot(dx))/SpeedOfLightAU; 
    out->tdbEmit2[oldN] = out->tEmit2[oldN] + f.tdb0;

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


std::map<Gravity, string> gravityNamesTracklet = {
  {Gravity::INERTIAL, "INERTIAL"},
  {Gravity::BARY, "BARYCENTER"},
  {Gravity::GIANTS, "GIANTS"}
  };
std::map<string, Gravity> gravityTypesTracklet = {
  {"INERTIAL",Gravity::INERTIAL},
  {"BARYCENTER", Gravity::BARY},
  {"GIANTS", Gravity::GIANTS}
};


void
FitterTracklet::save(std::ostream& os, string comment) const {
  if (!abgIsFit) {
    cerr << "WARNING: saving an invalid ABG fit" << endl;
  }
  // print comment if there is one
  if (!comment.empty()) {
    os << "# " << comment << endl;
  }
  // Save config info: gravity type
  os << "# Gravity type:" << endl;
  os << gravityNamesTracklet[grav] << endl;

  // Binding constraint
  os << "# Binding constraint strength " << endl;
  os << bindingConstraintFactor << endl;

  // Gamma priors
  os << "# Gamma prior " << endl;
  os << gamma0 << " " << gammaPriorSigma << endl;
  
  // write frame
  os << "# Coordinate frame " << endl;
  f.write(os);

  // write ABG
  ABG::writeHeader(os);
  getABG().write(os);

  // write covariance
  ABGCovariance cov(A.inverse());
  cov.write(os,5);

  // save observations
  os << "#\n# Observations " << endl;
  os << "#  TDB1  |  ThetaX1   |   ThetaY1  |  TDB2  |  ThetaX2   |   ThetaY2  |  covX1X1 Y1Y1 X2X2 Y2Y2 X1Y1 X1X2 X1Y2 Y1X2 Y1Y2 X2Y2 |  Observatory posn X1 Y1 Z1 |  Observatory posn X2 Y2 Z2" << endl;
  int n = tObs1.size();
  os << n << " observations" << endl;
  stringstuff::StreamSaver ss(os);
  Matrix44 covar;
  Matrix44 invcov;
  for (int i=0; i<n; i++) {
    invcov(0,0) = invCovArr(0,i);
    invcov(1,1) = invCovArr(1,i);
    invcov(2,2) = invCovArr(2,i);
    invcov(3,3) = invCovArr(3,i);
    invcov(0,1) = invcov(1,0) = invCovArr(4,i);
    invcov(0,2) = invcov(2,0) = invCovArr(5,i);
    invcov(0,3) = invcov(3,0) = invCovArr(6,i);
    invcov(1,2) = invcov(2,1) = invCovArr(7,i);
    invcov(1,3) = invcov(3,1) = invCovArr(8,i);
    invcov(2,3) = invcov(3,2) = invCovArr(9,i);
    covar = invcov.inverse();

    os << std::setprecision(9)
       << tObs1[i] << " "
       << thetaX1[i] << " "
       << thetaY1[i] << " "
       << tObs2[i] << " "
       << thetaX2[i] << " "
       << thetaY2[i] << " "
       << std::setprecision(4)
       << cov(0,0) << " " << cov(1,1) << " " << cov(2,2) << " " << cov(3,3) << " "
       << cov(0,1) << " " << cov(0,2) << " " << cov(0,3) << " "
       << cov(1,2) << " " << cov(1,3) << " " << cov(2,3) << " "
       << std::setprecision(9)
       << xE1(i,0) << " "
       << xE1(i,1) << " "
       << xE1(i,2)
       << xE2(i,0) << " "
       << xE2(i,1) << " "
       << xE2(i,2)

       << endl;
  }
}

void
FitterTracklet::restore(std::istream& is) {
  // Get gravity type
  string buffer;
  stringstuff::getlineNoComment(is,buffer);
  {
    string g;
    std::istringstream iss(buffer);
    iss >> g;
    try {
      grav = gravityTypesTracklet.at(g);
    } catch (std::out_of_range& e) {
      throw std::runtime_error("Unknown Gravity type in Fitter restore: " + g);
    }
  }
  
  // Get priors for binding and gamma
  stringstuff::getlineNoComment(is,buffer);
  {
    std::istringstream iss(buffer);
    double d;
    iss >> d;
    setBindingConstraint(d);
  }
  stringstuff::getlineNoComment(is,buffer);
  {
    std::istringstream iss(buffer);
    double d1,d2;
    iss >> d1 >> d2;
    setGammaConstraint(d1,d2);
  }

  // Get frame
  stringstuff::getlineNoComment(is,buffer);
  {
    std::istringstream iss(buffer);
    Frame f_;
    f_.read(iss);
    setFrame(f_);
  }

  // Read ABG and covariance
  ABG abg_;
  ABGCovariance cov_;
  stringstuff::getlineNoComment(is,buffer);
  {
    std::istringstream iss(buffer);
    abg_.read(iss);
  }
  cov_.read(is);

  // Now read the observation data
  int nObs;
  stringstuff::getlineNoComment(is,buffer);
  {
    std::istringstream iss(buffer);
    iss >> nObs;
  }
  // Create arrays
  DVector tObs1_(nObs);
  DVector thetaX1_(nObs);
  DVector thetaY1_(nObs);
  DMatrix xE1_(nObs,3);
  DVector tObs2_(nObs);
  DVector thetaX2_(nObs);
  DVector thetaY2_(nObs);
  DMatrix xE2_(nObs,3);
  DMatrix covFull_(nObs, 10);


  // Read each line
  for (int i=0; i<nObs; i++) {
    stringstuff::getlineNoComment(is, buffer);
    std::istringstream iss(buffer);
    iss >> tObs1_[i]
	>> thetaX1_[i] >> thetaY1_[i]
  >> tObs2_[i]
  >> thetaX2_[i] >> thetaY2_[i] 
	>> covFull_(i,0) >> covFull_(i,1) >> covFull_(i,2)
  >> covFull_(i,3) >> covFull_(i,4) >> covFull_(i,5)
  >> covFull_(i,6) >> covFull_(i,7) >> covFull_(i,8)
  >> covFull_(i,9) 
	>> xE1_(i,0) >> xE1_(i,1) >> xE1_(i,2)
  >> xE2_(i,0) >> xE2_(i,1) >> xE2_(i,2);

  }
  setObservationsInFrame(tObs1_, tObs2_, thetaX1_, thetaY1_, thetaX2_, thetaY2_, covFull_, xE1_, xE2_);

  // Save the solution
  setABG(abg_, cov_);
  
}
