// Orbit integration code

#include "AstronomicalConstants.h"
#include "Ephemeris.h"
#include "Trajectory.h"

using namespace orbits;
using namespace astrometry;

/**
In this class the xfwd/xbwd arrays have at element i
the position at time t=tdb0 +- i*dt, respectively.

Likewise the [va][fwd|bwd] arrays hold v*dt and 0.5*a*dt*dt
at time tdb0 +- i*dt.  Note that the v's are actually calculated
at the midpoints of timestep during leapfrog integration (these
are vnext_[fwd|bwd].  But we save an (approximate) value at the
start of timestep by adding half the acceleration.

With this convention, a quadratic interpolation of the position
at a fraction f of the way from step i*dt to (i+1)*dt is
x + f*v + f^2 * a.

A velocity estimate would be v + 2*f*a.

This interpolation scheme yields continuous x agrees with
leapfrog at the nodes.  The velocity also agrees with leapfrog
at times (i+1/2), but will have a discontinuity at the nodes
of 
step(v) = 0.5*dt* (a(i*dt)-a((i+1)*dt)) = 0.5 * dt^2 * da/dt
and the calculated acceleration will be a step function that
is off by 1/2 time step.  Could take the time to make a 
higher-order interpolator...
**/


Trajectory::Trajectory (const Ephemeris& ephem_,
			const State& s0,
			Gravity grav_,
			double dt_): ephem(ephem_),
				     x0(s0.x.getVector()),
				     v0(s0.v.getVector()),
				     tdb0(s0.tdb),
				     dt(dt_),
				     grav(grav_)
{
  if (grav != INERTIAL) {
    // Put initial positions in caches
    xfwd.push_back(x0);
    xbwd.push_back(x0);
    // Put in the initial half time step for the leap frog
    auto dv = deltaV(x0, tdb0);
    vnext_fwd = v0 + 0.5*dv;
    vnext_bwd = v0 - 0.5*dv;
    vfwd.push_back(v0*dt);
    vbwd.push_back(v0*dt);
    afwd.push_back(dv*0.5*dt);
    abwd.push_back(dv*0.5*dt);
  }
}

void
Trajectory::span(double tdb) const {
  // How many time steps fwd / back must we go?
  double tstep = (tdb - tdb0)/dt;

  if (tstep >= xfwd.size()) {
    // Extend cache forward
    int nfwd = static_cast<int> (ceil(tstep));
    double tdb = tdb0 + dt * xfwd.size();
    for (int i=xfwd.size(); i<=nfwd; i++, tdb+=dt) {
      // Leap frog forward
      xfwd.push_back( xfwd.back() + dt * vnext_fwd);
      auto dv = deltaV(xfwd.back(), tdb);
      afwd.push_back(0.5*dt*dv);
      vfwd.push_back( dt*(vnext_fwd + 0.5*dv));
      vnext_fwd += dv;
    }
  }

  if (tstep < -xbwd.size()) {
    // Extend cache backward
    int nbwd = static_cast<int> (ceil(-tstep));
    double tdb = tdb0 - dt * xbwd.size();
    for (int i=xbwd.size(); i<=nbwd; i++, tdb-=dt) {
      // Leap frog backward
      xbwd.push_back( xbwd.back() - dt * vnext_bwd);
      auto dv = deltaV(xbwd.back(), tdb);
      abwd.push_back(0.5*dt*dv);
      // Save negated velocity since LUT is using (-t)
      vbwd.push_back( (-dt)*(vnext_bwd - 0.5*dv));
      vnext_bwd -= dv;
    }
  }
}
		      
linalg::Matrix<double>
Trajectory::position(const linalg::Vector<double>& tdb,
		     linalg::Matrix<double>* velocity) const {
  linalg::Matrix<double> out(3,tdb.size());
  if (velocity)
    velocity->resize(3,tdb.size());
  if (tdb.size()==0) return out;
  if (grav==INERTIAL) {
    // Inertial motion is just linear algebra
    for (int i=0; i<tdb.size(); i++)
      out.col(i) = x0 + (tdb[i]-tdb0)*v0;
    velocity->colwise() = v0;
  } else {
    // Grow the integration of the orbit
    span(tdb.minCoeff());
    span(tdb.maxCoeff());

    // Interpolate all positions, backward ones first
    linalg::Vector<double> tstep = (tdb.array()-tdb0)/dt;
    auto tabs = tstep.cwiseAbs();
    linalg::Vector<double> istep = tabs.array().floor(); // Integer and fraction parts of |tstep|
    linalg::Vector<double> f = tabs - istep;  
    int i;
    for (i=0 ; i<tstep.size(); i++) {
      int i0 = static_cast<int> (istep[i]); // Index of time step nearer 0
      if (tstep[i]<0.) {
	// Use backward integration
	out.col(i) = xbwd[i0] + f[i]*(vbwd[i0] + f[i]*abwd[i0]);
	if (velocity)
	  // Negate the velocity because LUT is function of (-t)
	  velocity->col(i) = -(vbwd[i0] + (2.*f[i])*abwd[i0]);
      } else {
	// Use forward integration
	out.col(i) = xfwd[i0] + f[i]*(vfwd[i0] + f[i]*afwd[i0]);
	if (velocity)
	  velocity->col(i) = vfwd[i0] + (2.*f[i])*afwd[i0];
      } 
    }
    if (velocity) (*velocity)/=dt;
  }
  return out;
}

astrometry::CartesianICRS
Trajectory::position(double tdb,
		     astrometry::CartesianICRS* velocity) const {
  linalg::Vector<double> tvec(1,tdb);
  if (velocity) {
    linalg::Matrix<double> vmat(3,1);
    auto m = position(tvec, &vmat);
    *velocity = CartesianICRS(vmat(0,0), vmat(1,0), vmat(2,0));
    return CartesianICRS(m(0,0), m(1,0), m(2,0));
  } else {
    auto m = position(tvec);
    return CartesianICRS(m(0,0), m(1,0), m(2,0));
  }
}

astrometry::Vector3
Trajectory::deltaV(const Vector3& x, double tdb) const {
  Vector3 out(0.);
  if (grav==INERTIAL) {
    // no acceleration
  } else if (grav==BARY) {
    // Coordinates are already barycentric.
    double scale = -SolarSystemGM * dt * pow(x.dot(x), -1.5);
    out = scale * x;
  } else if (grav==GIANTS) {
    // Sum over Sun and giants
    static vector<double> gm = {JupiterGM, SaturnGM, UranusGM, NeptuneGM,
				GM + MercuryGM + VenusGM + EarthMoonGM + MarsGM};
    static vector<int> bodies = {orbits::JUPITER, orbits::SATURN, orbits::URANUS,
				 orbits::NEPTUNE, orbits::SUN};
    for (int i=0; i<bodies.size(); i++) {
      astrometry::Vector3 dx = x - ephem.position(bodies[i], tdb).getVector();
      out -= (gm[i] * dt * pow(dx.dot(dx), -1.5)) * dx;
    }
  }
  return out;
}
      
