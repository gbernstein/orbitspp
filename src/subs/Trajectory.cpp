// Orbit integration code

#include "AstronomicalConstants.h"
#include "Ephemeris.h"
#include "Trajectory.h"

using namespace orbits;
using namespace astrometry;

Trajectory::Trajectory (const Ephemeris& ephem_,
			const State& s0,
			double dt_,
			Gravity grav_): ephem(ephem_),
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
    vfwd.push_back( v0 + 0.5*dv);
    vbwd.push_back( v0 - 0.5*dv);
  }
}

void
Trajectory::span(double tdbMin, double tdbMax) const {
  // How many time steps fwd / back must we go?
  int nbwd = static_cast<int> (ceil((tdb0-tdbMin)/dt));
  int nfwd = static_cast<int> (ceil((tdbMax-tdb0)/dt));

  if (nfwd >= xfwd.size()) {
    // Extend cache forward
    double tdb = tdb0 + dt * xfwd.size();
    for (int i=xfwd.size(); i<=nfwd; i++, tdb+=dt) {
      // Leap frog forward
      xfwd.push_back( xfwd.back() + dt * vfwd.back() );
      vfwd.push_back( vfwd.back() + deltaV(xfwd.back(), tdb));
    }
  }

  if (nbwd >= xbwd.size()) {
    // Extend cache backward
    double tdb = tdb0 - dt * xbwd.size();
    for (int i=xbwd.size(); i<=nbwd; i++, tdb-=dt) {
      // Leap frog backward
      xbwd.push_back( xbwd.back() - dt * vbwd.back() );
      vbwd.push_back( vbwd.back() - deltaV(xbwd.back(), tdb));
    }
  }
}
		      
linalg::Matrix<double>
Trajectory::position(const linalg::Vector<double>& tdb) const {
  linalg::Matrix<double> out(3,tdb.size());
  if (tdb.size()==0) return out;
  if (grav==INERTIAL) {
    // Inertial motion is just linear algebra
    for (int i=0; i<tdb.size(); i++)
      out.col(i) = x0 + (tdb[i]-tdb0)*v0;
  } else {
    // Grow the integration of the orbit
    span(tdb[0], tdb[tdb.size()-1]);

    // Interpolate all positions, backward ones first
    linalg::Vector<double> tsteps = tdb;
    for (int i=0; i<tsteps.size(); i++)
      tsteps[i] = (tsteps[i]-tdb0) / dt;
    
    int i;
    for (i=0 ;tsteps[i]<0 && i<tsteps.size(); i++) {
      int i0 = static_cast<int> (-floor(tsteps[i])); // Index of time before
      out.col(i) = xbwd[i0] + (i0+tsteps[i])*dt*vbwd[i0-1];
    }

    // Now forward ones
    for ( ; i<tsteps.size(); i++) {
      int i0 = static_cast<int> (floor(tsteps[i])); // Index of time before
      out.col(i) = xfwd[i0] + (tsteps[i]-i0)*dt*vfwd[i0];
    }
  }
  return out;
}

astrometry::CartesianICRS
Trajectory::position(double tdb) const {
  linalg::Vector<double> vv(1,tdb);
  auto m = position(vv);
  return CartesianICRS(m(0,0), m(1,0), m(2,0));
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
      
