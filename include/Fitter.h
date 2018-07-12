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
  public:
    Fitter(const Ephemeris& eph_, Gravity grav_=Gravity::GIANTS);
    // Initialize fitter with ephemeris and choice of gravity approximation

    // Give Fitter a previous result (can be used for prediction without Observations)
    void setABG(const ABG& abg_, const ABGCovariance& cov_);
    
    // Add an Observation
    void addObservation(const Observation& obs);

    // Ingest a sequence of MPC-style observations from stream
    void readMPCObservations(istream& is);

    // Number of observations
    int nObservations() const {return observations.size();}

    // Set reference frame to the one given, put observations into this frame
    void setFrame(const Frame& f_);
    const Frame& getFrame() const {return f;}

    // Choose a reference frame as the location and direction of
    // the chosen observation.  If negative number,
    // nearest observation to mean MJD is used.
    void chooseFrame(int obsNumber = 0);
    
    // Set abg with simple linear fit: inertial orbit, gdot=0, 1 Newton iteration
    // from current starting point.
    void setLinearOrbit();

    // Place a Gaussian prior on gamma with the given mean and sigma.
    // sigg<=0 will disable prior.
    void setGammaConstraint(double g, double sigg) {
      gamma0=g;
      gammaPriorSigma = sigg;
    }

    // Place prior of f * (v/v_esc)^2, mild push to zero velocity away from escaping.
    // Parabolic orbit will get chisq penalty = f.
    void setBindingConstraint(double f=1.) {
      bindingConstraintFactor = f;
    }

    void newtonFit(double chisqTolerance=0.01, bool dump=false);

    // Print residuals (in arcsec) and chisq contributions per point
    void printResiduals(std::ostream& os) const;
    
    // Obtain fitting results
    const ABG& getABG(bool invalidOK=false) const {
      if (!abgIsValid && !invalidOK)
	throw std::runtime_error("ERROR: Fitter::getABG is not getting converged result");
      return abg;
    }
    double getChisq() const {return chisq;}  // Chisq at current abg
    ABGCovariance getInvCovarABG() const {
      // Inverse covariance of ABG from last fit
      if (!abgIsValid) throw std::runtime_error("ERROR: Fitter::getInvCovarABG is not getting converged result");
      return A;}  
    Elements getElements() const;
    ElementCovariance getElementCovariance() const;

    const Trajectory& getTrajectory() const {return *fullTrajectory;}

    // Forecast position using current fit.  Cov matrix elements given if filled:
    void predict(const DVector& t_obs,    // Time of observations, relative to tdb0
		 const DMatrix& earth,  // Observation coordinates, in our frame
		 DVector* xOut,         // Angular coordinates, in our frame
		 DVector* yOut,
		 DVector* covXX = nullptr,   // Covar matrix of coordinates
		 DVector* covXY = nullptr,   // (not computed if nullptrs)
		 DVector* covYY = nullptr) const;
		 
    
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
    DVector tObs;  // Times of observations, in Julian years since reference time
    DVector tEmit; // Times of light emission, in Julian years since reference time
    DVector tdbEmit; // TDB at time of light emission per observation, years since J2000 TDB 
    DVector thetaX; // Observed angles in our frame
    DVector thetaY; // Observed angles in our frame
    DVector invcovXX; // Inverse cov matrix of observations
    DVector invcovXY; 
    DVector invcovYY;
    DMatrix xE;    // Earth positions in our ref frame

    // Model/fitting information:
    Trajectory* fullTrajectory;
    DVector thetaXModel;  // Positions predicted by current ABG
    DVector thetaYModel;  // Positions predicted by current ABG
    DMatrix xGrav; // Gravitational component of motion in our frame
    DMatrix dThetaXdABG;  // Derivatives of x position wrt abg
    DMatrix dThetaYdABG;  // Derivatives of y position wrt abg

    // Fitting results
    ABG abg;         // Orbit in abg basis
    double chisq;  // Chisq at current abg
    Vector6 b;  // -1/2 d(chisq)/d(abg)
    Matrix66 A; // 1/2 d^2(chisq)/d(abg)^2

    // Configuration and prior
    Gravity grav;       // Use giant planets, SS barycenter, or no gravity?
    double bindingConstraintFactor;  // Chisq penalty for marginally unbound orbit
    double gamma0;	// Nominal gamma and uncertainty when using gamma prior
    double gammaPriorSigma;

    // Some state indicators
    bool frameIsSet; // The reference frame has been set.
    bool abgIsValid; // The ABG is a solution to observations and its (inv) covariance A is valid.
  };
} // namespace orbits
#endif
