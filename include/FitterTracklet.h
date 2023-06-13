// Class that manages orbit fitting
#ifndef FITTERTRACKLET_H
#define FITTERTRACKLET_H

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

  class FitterTracklet{
  public:
    FitterTracklet(const Ephemeris& eph_, Gravity grav_=Gravity::GIANTS);

    ~FitterTracklet() {
      if (fullTrajectory) delete fullTrajectory;
    }
    
    // Define an exception class for non-convergence of fits
    class NonConvergent: public std::runtime_error {
    public:
      NonConvergent(string msg=""): std::runtime_error("Fit Non-Convergence: " + msg) {}
    };

    // Save and restore observations & ABG fit to/from a stream.
    void save(std::ostream& os, string comment="") const;
    void restore(std::istream& is);
    
    // Give Fitter a previous result (can be used for prediction without Observations)
    void setABG(const ABG& abg_, const ABGCovariance& cov_);
    
    // Add an Observation
    void addObservation(const Observation& obs);

    // Ingest a sequence of MPC-style observations from stream
    void readMPCObservations(istream& is);

    // Add an Observation
    void addTracklet(const Tracklet& track);


    // Replace observations with data that is already in desired Frame
    void setObservationsInFrame(const DVector& tObs1_,    // TDB since reference time
        const DVector& tObs2_,
				const DVector& thetaX1_,  // Observed positions
				const DVector& thetaY1_,
        const DVector& thetaX2_,  // Observed positions
        const DVector& thetaY2_,
				const DMatrix& covFull_,   // Covariance of observed posns
				const DMatrix& xE1_,     // Observatory posn at observations
        const DMatrix& xE2_);     // Observatory posn at observations


    // Number of observations
    int nObservations() const {return frameIsSet ? tObs1.size() : tracks.size();}

    // Set reference frame to the one given, put observations into this frame
    void setFrame(const Frame& f_);
    const Frame& getFrame() const {return f;}

    // Choose a reference frame as the location and direction of
    // the chosen observation.  If negative number,
    // nearest observation to mean MJD is used.
    void chooseFrame(int obsNumber = 0);
    
    void chiSqNoCovariance();

    // Set abg with simple linear fit: inertial orbit, gdot=0, 1 Newton iteration
    // from current starting point.
    void setLinearOrbit();
    void setSingleOrbit(double chisqTolerance=0.01);


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

    void newtonFit(double chisqTolerance=0.01, bool dump=false, bool doCovariances=true);

    // Print residuals (in arcsec) and chisq contributions per point
    void printResiduals(std::ostream& os);
    
    // Obtain fitting results
    const ABG& getABG(bool invalidOK=false) const {
      if (!abgIsFit && !invalidOK)
	throw std::runtime_error("ERROR: Fitter::getABG is not getting converged result");
      return abg;
    }
    double getChisq(); // Chisq at current abg
    int getDOF() const {return 2*2*tObs1.size() - 6;}  // DOF does not count priors
    ABGCovariance getInvCovarABG() const {
      // Inverse covariance of ABG from last fit
      if (!abgIsFit) throw std::runtime_error("ERROR: Fitter::getInvCovarABG is not getting converged result");
      return A;}  
    // Get elements - barycentric by default, heliocentric if helio=true
    Elements getElements(bool helio=false) const;
    ElementCovariance getElementCovariance(bool helio=false) const;
    // Get elements at a specified tdb, given relative to J2000.
    Elements getElements(double tdb, bool helio=false) const;
    ElementCovariance getElementCovariance(double tdb, bool helio=false) const;

    const Trajectory& getTrajectory() const {return *fullTrajectory;}

    // Forecast position using current fit.  Cov matrix elements given if filled:
    void predict(const DVector& t_obs,    // Time of observations, relative to tdb0
		 const DMatrix& earth,  // Observation coordinates, in our frame
		 DVector* xOut,         // Angular coordinates, in our frame
		 DVector* yOut,
		 DVector* covXX = nullptr,   // Covar matrix of coordinates
		 DVector* covYY = nullptr,   // (not computed if nullptrs)
		 DVector* covXY = nullptr) const;
		 
    // Forecast ICRS state vector at arbitrary time, optionally fill covariance too.
    State predictState(const double tdb,    // Dynamical time, relative to tdb0
		       Matrix66* stateCov=nullptr) const;
    
    // Return a new Fitter that has orbit updated using an additional observation
    // The additional data is assumed already in desired Frame.
    // Use updated orbit to recalculate non-inertial motion if newGravity=true,
    // otherwise keep trajectory from previous fit.
    FitterTracklet* augmentObservation(double tObs1_,    // TDB since reference time
             double tObs2_,
			       double thetaX1_,  // Observed position
			       double thetaY1_,
             double thetaX2_,
             double thetaY2_,
             DVector& cov_,
			       const Vector3& xE1_,     // Observatory posn at observation
             const Vector3& xE2_,     // Observatory posn at observation

			       bool newGravity=false) const; 
    EIGEN_NEW


  private:
    const Ephemeris& eph; // Solar system Ephemeris

    vector<Tracklet, Eigen::aligned_allocator<Tracklet> > tracks;
    
    void iterateTimeDelay(); // Update tEmit based on light-travel time in current orbit
    void calculateGravity(); // Calculate non-inertial terms from current ABG
    void createTrajectory(); // Initialize the orbit integrator from current ABG

    // Calculate angular positions and their derivs wrt ABG
    void calculateOrbit(bool doDerivatives=true); 

    // Calculate chisq (and derivatives) at current abg.
    // (Calls calculateOrbit if needed)
    void calculateChisq(bool doDerivatives=true, bool doCovariances=true); 


    // Throw an exception if we have a clearly invalid ABG
    void abgSanityCheck() const;
      
    void resizeArrays(int n); // Set all of arrays below to desired sizes

    Frame f;   // Reference frame for our coordinates

    // Observational information:
    DVector tObs1;  // Times of observations, in Julian years since reference time
    DVector tEmit1; // Times of light emission, in Julian years since reference time
    DVector tdbEmit1; // TDB at time of light emission per observation, years since J2000 TDB 
    DVector thetaX1; // Observed angles in our frame
    DVector thetaY1; // Observed angles in our frame
    DMatrix xE1;    // Earth positions in our ref frame
    
    DVector tObs2;  // Times of observations, in Julian years since reference time
    DVector tEmit2; // Times of light emission, in Julian years since reference time
    DVector tdbEmit2; // TDB at time of light emission per observation, years since J2000 TDB 
    DVector thetaX2; // Observed angles in our frame
    DVector thetaY2; // Observed angles in our frame
    DMatrix xE2;    // Earth positions in our ref frame

    DMatrix invCovArr;


    // Model/fitting information:
    Trajectory* fullTrajectory;
    DVector thetaXModel1;  // Positions predicted by current ABG
    DVector thetaYModel1;  // Positions predicted by current ABG
    DMatrix xGrav1; // Gravitational component of motion in our frame
    DMatrix dThetaX1dABG;  // Derivatives of x position wrt abg
    DMatrix dThetaY1dABG;  // Derivatives of y position wrt abg

    DVector thetaXModel2;  // Positions predicted by current ABG
    DVector thetaYModel2;  // Positions predicted by current ABG
    DMatrix xGrav2; // Gravitational component of motion in our frame
    DMatrix dThetaX2dABG;  // Derivatives of x position wrt abg
    DMatrix dThetaY2dABG;  // Derivatives of y position wrt abg


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
    bool abgIsFit; // The ABG is a solution to the observations
    bool positionsAreValid; // theta[XY]Model correspond to current ABG.
    bool positionDerivsAreValid; // dTheta[XY]dABG correspond to current ABG.
    bool chisqIsValid;   // Chisq corresponds to current ABG.
    bool chisqDerivsAreValid;  // Chisq derivatives (b,A) correspond to current ABG.

    // Reset state flags for new solution
    void newABG() {
      abgIsFit = true;
      positionsAreValid = positionDerivsAreValid = false;
      chisqIsValid = chisqDerivsAreValid = false;
    }

    // Reset state flags for new data
    void newData() {
      newABG();
      abgIsFit = false;
    }

  };
} // namespace orbits
#endif