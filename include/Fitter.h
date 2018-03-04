// Class that manages orbit fitting
#ifndef FITTER_H
#define FITTER_H

#include "Astrometry.h"
#include "OrbitTypes.h"
#include "Ephemeris.h"
#include "Trajectory.h"

namespace orbits {

  class Fitter {
    typedef linalg::Vector<double> Vector;
    typedef linalg::Matrix<double> Matrix;
  public:
    Fitter(const Ephemeris& eph_);

    // Ingest a sequence of MPC-style observations from stream
    void readObservations(istream& is);

    void setFrame(const Frame& f_);

    // Set reference frame to the one given
    void chooseFrame();
    
    // Choose a good reference frame from observations
    void orbitDerivatives();
    void setLinearOrbit();
    
    
  private:
    vector<MPCObservation> observations;
    void setupObservations();  // Calculate projected positions, Earth posns

    Frame f;   // Reference frame for our coordinates
    Vector tdb;    // TDB's of observations
    Vector thetaX; // Observed angles in our frame
    Vector thetaY; // Observed angles in our frame
    Vector invcovXX; // Inverse cov matrix of observations
    Vector invcovXY; 
    Vector invcovYY;
    Matrix xE;    // Earth positions in our ref frame
    Matrix xGrav; // Gravitational component of motion in our frame
    const Ephemeris& eph; // Solar system Ephemeris
    ABG abg;         // Orbit in abg basis
    Vector b;  // First derivs of chisq wrt abg
    Matrix A; // 2nd derivs of chisq wrt abg.
    Trajectory* inertialTrajectory;
    Trajectory* fullTrajectory;
    bool useGiants;       // Use giant planets or just SS barycenter?
    double energyConstraintFactor;
    double gamma0;	// Nominal gamma and uncertainty when using gamma prior
    double gammaPriorSigma;
  };
} // namespace orbits
#endif
