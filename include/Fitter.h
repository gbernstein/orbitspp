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
    Fitter(const Ephemeris& eph_);

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
    
    ABG abg;         // Orbit in abg basis
    linalg::SMatrix<double,6,6> abgInvCov; // Inverse covariance of abg
    
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
    void calculateChisqDerivs(); // Calculate chisq and derivs wrt abg

    Frame f;   // Reference frame for our coordinates

    // Observational information:
    Vector tdb;    // TDB's of observations
    Vector dt;    // TDB's of observations, relative to Frame's reference time
    Vector tEmit;    //  TDB at time of light emission per observation
    Vector thetaX; // Observed angles in our frame
    Vector thetaY; // Observed angles in our frame
    Vector invcovXX; // Inverse cov matrix of observations
    Vector invcovXY; 
    Vector invcovYY;
    Matrix xE;    // Earth positions in our ref frame

    // Model/fitting information:
    Trajectory* inertialTrajectory;
    Trajectory* fullTrajectory;
    Vector thetaXModel;  // Positions predicted by current ABG
    Vector thetaYModel;  // Positions predicted by current ABG
    Matrix xGrav; // Gravitational component of motion in our frame
    Matrix dThetaXdABG;  // Derivatives of x position wrt abg
    Matrix dThetaYdABG;  // Derivatives of y position wrt abg
    double chisq;  // Chisq at current abg
    Vector b;  // First derivs of chisq wrt abg
    Matrix A; // 2nd derivs of chisq wrt abg.

    // Configuration and prior
    bool useGiants;       // Use giant planets or just SS barycenter?
    double energyConstraintFactor;
    double gamma0;	// Nominal gamma and uncertainty when using gamma prior
    double gammaPriorSigma;
  };
} // namespace orbits
#endif
