// Orbit integration code

#include "AstronomicalConstants.h"
#include "Ephemeris.h"
#include "Trajectory.h"


using namespace orbits;
using namespace astrometry;

/**
In this class the xvaLUT array has at element i
the (x,vv,aa) at time t=tdb0 +- i*dt, respectively.

The x holds the coordinate; vv holds v*dt; and aa holds 0.5*a*dt*dt,
so that the interpolation across the time step is
x(tdb0 + (i+f)*dt) = x_i + f*vv_i + f*f*aa_i
v(tdb0 + (i+f)*dt)  * dt = v_i + 2*f*aa_i
where f is the fractional part of the time step.

This interpolation scheme yields continuous x agrees with
leapfrog at the nodes.  The velocity also agrees with leapfrog
at times (i+1/2), but will have a discontinuity at the nodes
of 
step(v) = 0.5*dt* (a(i*dt)-a((i+1)*dt)) = 0.5 * dt^2 * da/dt
and the calculated acceleration will be a step function that
is off by 1/2 time step.  Could take the time to make a 
higher-order interpolator...

The leapfrog wants to work with the velocity at the midpoint of
the interval, not at the interval (in fact the formula above produces
a slight discontinuity in v at the nodes).  Let's call the velocity at
the midpoint of the interval i->i+1 as vmid_i.  The leapfrog is
x -> x + vmid_i*dt
vmid*dt -> vmid_i*dt + a_i*dt*dt = vmid_i*dt + dv_i
where dv = a_i*dt, a_i is the gravity evaluated at time i, position x_i

and this is the way we will do the integration.  But the tabulated
values differ slightly from vmid.  Equating the leapfrog step to the 
interpolation formula for f=1 gives
vmid_i = vv_i + aa_i, or
vv_i = vmid_i - aa_i

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
  if (grav == INERTIAL) {
    return;
  } else {
    // Put in initial values for fwd and bwd leapfrog tables
    auto dv = deltaV(x0, tdb0);
    Matrix33 xva;
    xva.col(0) = x0;
    xva.col(1) = v0*dt;
    xva.col(2) = dv*0.5*dt;
    xvaLUT.extend(0,xva);
  }
}

void
Trajectory::span(double tdb) const {
  // How many time steps fwd / back must we go?
  double tstep = (tdb - tdb0)/dt;
  int s = xvaLUT.iEnd();  // Size might increase before we append
  if (tstep >= s) {
    // Extend cache forward
    vector<Matrix33, Eigen::aligned_allocator<Matrix33>> addFwd;
    int nfwd = static_cast<int> (ceil(tstep));
    double tdb = tdb0 + dt * s;
    auto xva = xvaLUT[s-1];
    // half-step v to the midpoint, which
    // is what we'll use for leapfrogging
    Vector3 vmid = xva.col(1) + xva.col(2);  
    
    for (int i=s; i<=nfwd; i++, tdb+=dt) {
      // Leap frog forward to x_{i+1}
      xva.col(0) += vmid;
      // Calculate and save acceleration
      auto dv = deltaV(xva.col(0),tdb)*dt;
      xva.col(2) = 0.5*dv;
      // Leapfrog the velocity, leaving half of it at the time step.
      xva.col(1) = vmid  + 0.5*dv;
      vmid += dv;
      // But run it back from midpoint to node for the interpolation table
      addFwd.push_back(xva);
    }
    // Extend the cache (this will block other writes)
    xvaLUT.extend(s,addFwd.begin(),addFwd.end());
  }

  // Now check for backward interpolations
  s = xvaLUT.iStart();
  if (tstep < s) {
    // Extend cache backward
    vector<Matrix33, Eigen::aligned_allocator<Matrix33>> addBwd;
    auto xva = xvaLUT[s];
    int nbwd = static_cast<int> (floor(tstep));
    double tdb = tdb0 + dt * (s-1);
    // half-step v to the midpoint, which
    // is what we'll use for leapfrogging
    Vector3 vmid = xva.col(1) - xva.col(2);  
    for (int i=s-1; i>=nbwd; i--, tdb-=dt) {
      // Leap frog backward
      xva.col(0) -= vmid;
      auto dv = deltaV(xva.col(0), tdb)*dt;
      xva.col(2) = 0.5*dv;
      // Leapfrog velocity
      xva.col(1) = vmid - 0.5*dv;
      vmid -= dv;
      // Run velocity back to node for interp table
      addBwd.push_back(xva);
    }
    // Extend the cache (this will block other writes)
    // Note we need to reverse the order since we integrated backwards.
    std::reverse(addBwd.begin(), addBwd.end());
    xvaLUT.extend(nbwd,addBwd.begin(),addBwd.end());
  }
}
		      
DMatrix
Trajectory::position(const DVector& tdb,
		     DMatrix* velocity) const {
  DMatrix out(tdb.size(),3);
  if (velocity) {
    velocity->resize(tdb.size(),3);
  }
  if (tdb.size()==0) return out;
  if (grav==INERTIAL) {
    // Inertial motion is just linear algebra
    for (int i=0; i<tdb.size(); i++)
      out.row(i) = (x0 + (tdb[i]-tdb0)*v0).transpose();
    if (velocity)
      velocity->rowwise() = v0.transpose();
  } else {
    // Grow the integration of the orbit
    span(tdb.minCoeff());
    span(tdb.maxCoeff());

    // Interpolate all positions
    DVector tstep = (tdb.array()-tdb0)/dt;
    DVector istep = tstep.array().floor(); // Integer and fraction parts of |tstep|
    DVector f = tstep - istep;  

    // Note that the SharedLUT takes care of parallel reads,
    // no need to worry about it here.
    // Create vectors used for interpolation
    Vector3 xf, vf;
    Matrix33 xva;
    xf[0] = 1.;
    vf[0] = 0.;
    vf[1] = 1.;
    for (int i=0 ; i<tstep.size(); i++) {
      int i0 = static_cast<int> (istep[i]); 
      xva = xvaLUT[i0];
      xf[1] = f[i];
      xf[2] = f[i]*f[i];
      out.row(i) = (xva * xf).transpose();
      if (velocity) {
	vf[2] = 2.*f[i];
	velocity->row(i) = (xva*vf).transpose();
      }
    }
    if (velocity) (*velocity)/=dt;  // we store v*dt, return just v
  }
  return out;
}

astrometry::CartesianICRS
Trajectory::position(double tdb,
		     astrometry::CartesianICRS* velocity) const {
  DVector tvec(1,tdb);
  if (velocity) {
    DMatrix vmat(1,3);
    auto m = position(tvec, &vmat);
    *velocity = CartesianICRS(vmat(0,0), vmat(0,1), vmat(0,2));
    return CartesianICRS(m(0,0), m(0,1), m(0,2));
  } else {
    auto m = position(tvec);
    return CartesianICRS(m(0,0), m(0,1), m(0,2));
  }
}

Vector3
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
      Vector3 dx = x - ephem.position(bodies[i], tdb).getVector();
      out -= (gm[i] * dt * pow(dx.dot(dx), -1.5)) * dx;
    }
  }
  return out;
}
      
// Return observed astrometric position from observer position/time.
astrometry::SphericalICRS
Trajectory::observe(double tdbObserve,
		    const astrometry::CartesianICRS& observer) const {
  double tEmit = tdbObserve;
  // Get the light-travel time
  CartesianICRS velocity;
  CartesianICRS target = position(tEmit, &velocity);
  target -= observer;
  double v_los = velocity.getVector().dot(target.getVector()) / target.radius();
  tEmit = tdbObserve - target.radius() / (SpeedOfLightAU + v_los);
  
  // Final position:
  target = position(tEmit);
  target -= observer;
  return SphericalICRS(target);
}

DMatrix
Trajectory::observe(const DVector& tdbObserve,
		    const DMatrix& observer) const {
  DVector tEmit = tdbObserve;
  // Get the light-travel time
  DMatrix velocity;
  DMatrix target = position(tEmit, &velocity);
  target -= observer;
  DVector targetRadius = target.rowwise().norm();
  DVector v_los = (velocity.array()*target.array()).rowwise().sum() /
    targetRadius.array();
  tEmit = tdbObserve.array() - targetRadius.array() / (v_los.array() + SpeedOfLightAU);
  
  // Final position:
  target = position(tEmit);
  target -= observer;
  // Return direction cosines
  targetRadius = target.rowwise().norm();
  target.col(0).array() /= targetRadius.array();
  target.col(1).array() /= targetRadius.array();
  target.col(2).array() /= targetRadius.array();
  return target;
}
