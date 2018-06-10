// Class that manages orbit fitting
#ifndef FITTER_H
#define FITTER_H

#include "Astrometry.h"
#include "OrbitTypes.h"
#include "Ephemeris.h"
#include "Trajectory.h"

// Need to following to make a vector containing structures with fixed-size
// Eigen arrays
#ifdef USE_EIGEN
#include <Eigen/StdVector>
#endif
namespace orbits {

  class Fitter {
    typedef linalg::Vector<double> Vector;
    typedef linalg::Matrix<double> Matrix;
  public:
    Fitter(const Ephemeris& eph_, Gravity grav_=Gravity::GIANTS);
    // Initialize fitter with ephemeris and choice of gravity approximation

    // Ingest a sequence of MPC-style observations from stream
    void readMPCObservations(istream& is);

    // Set reference frame to the one given, put observations into this frame
    void setFrame(const Frame& f_);

    // Choose a reference frame as the location and direction of
    // the chosen observation.  If negative number,
    // nearest observation to mean MJD is used.
    void chooseFrame(int obsNumber = 0);
    
    // Set abg with simple linear fit: inertial orbit, nominal gamma, gdot=0
    void setLinearOrbit(double nominalGamma=-1.);

    // Place a Gaussian prior on gamma with the given mean and sigma.
    // sigg<=0 will disable prior.
    void setGammaConstraint(double g, double sigg) {
      gamma0=g;
      gammaPriorSigma = sigg;
    }

    void newtonFit(double chisqTolerance=0.01);
    
    ABG abg;         // Orbit in abg basis
    linalg::SMatrix<double,6,6> abgInvCov; // Inverse covariance of abg
    double getChisq() const {return chisq;}  // Chisq at current abg
    Matrix getInvCovarABG() const {return A;}  // Inverse covariance of ABG from last fit
    
  private:
    const Ephemeris& eph; // Solar system Ephemeris

#ifdef USE_EIGEN
    vector<Observation, Eigen::aligned_allocator<Observation> > observations;
#else
    vector<Observation> observations;
#endif
    
    void iterateTimeDelay(); // Update tEmit based on light-travel time in current orbit.
    void calculateGravity(); // Calculate non-inertial terms from current ABG

    // Calculate positions and their derivs wrt ABG
    void calculateOrbitDerivatives(); 

    void calculateChisq(); // Calculate chisq at current abg
    void calculateChisqDerivatives(); // Calculate chisq and derivs wrt abg

    Frame f;   // Reference frame for our coordinates

    // Observational information:
    Vector tObs;  // Times of observations, in Julian years since reference time
    Vector tEmit; // Times of light emission, in Julian years since reference time
    Vector tdbEmit; // TDB at time of light emission per observation, years since J2000 TDB 
    Vector thetaX; // Observed angles in our frame
    Vector thetaY; // Observed angles in our frame
    Vector invcovXX; // Inverse cov matrix of observations
    Vector invcovXY; 
    Vector invcovYY;
    Matrix xE;    // Earth positions in our ref frame

    // Model/fitting information:
    Trajectory* fullTrajectory;
    Vector thetaXModel;  // Positions predicted by current ABG
    Vector thetaYModel;  // Positions predicted by current ABG
    Matrix xGrav; // Gravitational component of motion in our frame
    Matrix dThetaXdABG;  // Derivatives of x position wrt abg
    Matrix dThetaYdABG;  // Derivatives of y position wrt abg
    double chisq;  // Chisq at current abg
    Vector b;  // -1/2 d(chisq)/d(abg)
    Matrix A; // 1/2 d^2(chisq)/d(abg)^2

    // Configuration and prior
    Gravity grav;       // Use giant planets, SS barycenter, or no gravity?
    double energyConstraintFactor;
    double gamma0;	// Nominal gamma and uncertainty when using gamma prior
    double gammaPriorSigma;
  };
} // namespace orbits
#endif
