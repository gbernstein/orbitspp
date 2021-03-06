// Tranformation between state vector and orbital elements.
// Not good for unbound orbits!

#include "OrbitTypes.h"
#include "AstronomicalConstants.h"
#include "Astrometry.h"
#include "Ephemeris.h"

#include "Eigen/Geometry"     // For curl operation

using namespace astrometry;

namespace orbits {

  Elements
  getElements(const State& s,
	      bool heliocentric=false, const Ephemeris* ephem=nullptr) {
    // Return orbital elements for state vector.
    // Use solar mass and Sun if heliocentric=true, else full solar system
    // mass and barycenter.
    // The incoming state is assumed ICRS so we will convert to ecliptic
    // to get ecliptic J2000 elements.
  
    double mass = heliocentric ? GM : SolarSystemGM;
    double tdb = s.tdb;

    // Convert state to ecliptic
    Vector3 R, V;
    if (heliocentric) {
      // Remove solar motion from state
      if (!ephem)
	throw runtime_error("Need an Ephemeris for heliocentric getElements");
      auto sun = ephem->state(orbits::SUN, tdb);
      auto net = s;
      net.x -= sun.x;
      net.v -= sun.v;
      R = CartesianEcliptic(net.x).getVector();
      V = CartesianEcliptic(net.v).getVector();
    } else {
      R = CartesianEcliptic(s.x).getVector();
      V = CartesianEcliptic(s.v).getVector();
    }

    double rMagnitude, velSquare, radVelDotProduct;
    Vector3 eccentricityVector, angularMomentum, ascendingNode;
    double semimajor, eccentricity, inclination, longitudeOfAscendingNode;
    double semiLatusRectum;
    double hMagnitude, ascendingNodeMagnitude; /* magnitude of angular momentum */
    double ascEccDotProduct, argumentOfPerifocus;
    double xBar, yBar;
    double cosE, sinE, E1, E2;
    /* E1 and E2 are used to decide the quadrant of Eccentric Anomaly */
    double meanAnomaly, meanMotion, timeOfPerifocalPassage;

    rMagnitude = sqrt(R.dot(R));
    velSquare = V.dot(V);
    radVelDotProduct = R.dot(V);
    eccentricityVector = (velSquare/mass - 1/rMagnitude)*R
      - (radVelDotProduct/mass)*V;

    /* Angular mom is cross product of rad and vel */
    angularMomentum[0] = R[1]*V[2] - R[2]*V[1];
    angularMomentum[1] = R[2]*V[0] - R[0]*V[2];
    angularMomentum[2] = R[0]*V[1] - R[1]*V[0];

    /* Ascending Node vector is cross product of k-unit vector and angular momentum
       vector. k = [0 0 1] and h = [hx, hy, hz] */
    ascendingNode[0] = -angularMomentum[1];
    ascendingNode[1] = angularMomentum[0];
    ascendingNode[2] = 0.0;

    semimajor = 1/(2/rMagnitude - velSquare/mass);
    eccentricity = sqrt(eccentricityVector.dot(eccentricityVector));
    semiLatusRectum = angularMomentum.dot(angularMomentum)/mass;

    /* p = h-square by mu */
    hMagnitude = sqrt( angularMomentum.dot(angularMomentum));
    inclination = acos(angularMomentum[2]/hMagnitude); /* in radians here */
    ascendingNodeMagnitude = sqrt(ascendingNode.dot(ascendingNode));
    longitudeOfAscendingNode = acos(ascendingNode[0]/ascendingNodeMagnitude);

    /* Capital Omega in radians here */
    if (ascendingNode[1] < 0) longitudeOfAscendingNode = 
				2*PI - longitudeOfAscendingNode;
    /* ???could use atan2 here?? */
    ascEccDotProduct = ascendingNode.dot(eccentricityVector);
    argumentOfPerifocus = acos(ascEccDotProduct/
			       (ascendingNodeMagnitude*eccentricity)); 
    /* Small omega in radians here */
    if (eccentricityVector[2] < 0) argumentOfPerifocus = 
				     2*PI - argumentOfPerifocus;
    xBar = (semiLatusRectum - rMagnitude)/eccentricity;
    yBar = radVelDotProduct*sqrt(semiLatusRectum/mass)/eccentricity;

    if (eccentricity < 1) {
      /* Elliptical solution*/

      cosE = (xBar/semimajor) + eccentricity;
      sinE = yBar/(semimajor*sqrt(1-eccentricity*eccentricity));
      /* where semimajor*sqrt(1-eccentricity*eccentricity) is semiminor */
      double eccentricAnomaly = atan2(sinE,cosE);

      meanAnomaly = eccentricAnomaly - eccentricity*sinE; /* radians */
      meanMotion = sqrt(mass/(pow(semimajor,3.)));
      timeOfPerifocalPassage = tdb - meanAnomaly/meanMotion;
      /* This comes from M=n(t-T) where t is epoch time and T is time of perifocal
	 passage, in days */
    } else {
      /* hyperbolic solution */
      /* several ways to get the hyperbolic anomaly.  Try this one which has
	 no sign ambiguity:
	 tanh(H/2) = sqrt( (e-1)/(e+1))) tan (trueAnomaly/2)
         where tan(trueAnomaly) = yBar / xBar.

	 Or the much simpler r = a(1 - e*cosh(H)) which has sign ambiguity resolved by
	 sign of r dot v, but also is nearly degenerate near perihelion.  So maybe the
	 first one is better since we care more about numerical precision near peri than
	 when asymptoting.

	 Also have tan(trueAnomaly/2) = yBar / (rMagnitude + xBar)
      */
      double tmp = sqrt( (eccentricity-1)/(eccentricity+1) ) * yBar / (rMagnitude+xBar);
      double hyperbolicAnomaly = 2 * atanh(tmp);
      meanAnomaly = eccentricity * sinh(hyperbolicAnomaly) - hyperbolicAnomaly;
      meanMotion = sqrt(mass/(pow(-semimajor,3.)));
      timeOfPerifocalPassage = tdb - meanAnomaly/meanMotion;     
    }			

    Elements out;
    out[Elements::A] = semimajor;
    out[Elements::E] = eccentricity; 
    out[Elements::I] = inclination;
    out[Elements::LAN] = longitudeOfAscendingNode;
    out[Elements::AOP] = argumentOfPerifocus;
    out[Elements::TOP] = timeOfPerifocalPassage;

    return out; 
  } 

  State
  getState(const Elements& el, double tdb,
	   bool heliocentric=false, const Ephemeris* ephem=nullptr) {
    Vector3  r0, v0, r1, v1, r2, v2;
    double eccentricAnomaly;
    double meanAnomaly, meanMotion;
    double c, s, dt;

    double mass = heliocentric ? GM : SolarSystemGM;
    double e = el[Elements::E];
    double a = el[Elements::A];
    const double TOLERANCE = 0.001 * TIMESEC; // Tolerance is motion is 1 msec
  
    if (e == 1.)
      throw std::runtime_error("elements_to_xv not coded for parabolic orbits");
    if ( a*(1-e) < 0.)
      throw std::runtime_error("elements_to_xv obtained negative perihelion");

    if (e<1) {
      // Elliptical orbit
      meanMotion = pow(a,-1.5) * sqrt(mass);
      const double tol = TOLERANCE * meanMotion;
      meanAnomaly = (tdb - el[Elements::TOP]) * meanMotion;
      /* Put mean anomaly near 0 */
      meanAnomaly -= floor(meanAnomaly / TPI) * TPI;
      /* get the eccentric Anomaly from the mean anomaly using Newton iterations */
      const int JMAX=40;
      double ea = meanAnomaly; // starting guess for eccentricAnomaly
      double dea;  // error in ea
      for (int j=0; j<JMAX+1; j++) {
	if (j>=JMAX)
	  throw std::runtime_error("eccentricAnomaly solution took too long");
	dea = meanAnomaly + e*sin(ea) - ea;
	if (abs(dea) < tol)
	  break;  // Converged.
	ea += dea / (1 - e*cos(ea));  // Newton correction.
      }
	
      /*Coordinates and velocity in the system aligned with orbit: */
      c = cos(ea);
      s = sin(ea);
      double b = a * sqrt(1-e*e);
      r0[0] = a * (c - e);
      r0[1] = b * s;
      dt = ( 1 - e*c) / meanMotion;  // This is dt/dEccentricAnomaly
      v0[0] = -a * s / dt;
      v0[1] = b * c  / dt;
    } else {
      // Hyperbolic orbit
      meanMotion = pow(-a,-1.5) * sqrt(mass);
      const double tol = TOLERANCE * meanMotion;
      meanAnomaly = (tdb - el[Elements::TOP]) * meanMotion;
      // Solve for hyperbolic anomaly with Newton iterations
      const int JMAX=40;
      double h = log(1+abs(meanAnomaly)); // starting guess at hyperbolic anom
      double dh;  // error in hyperbolicAnomaly
      if (meanAnomaly < 0)
	h *= -1.;
      
      for (int j=0; j<JMAX+1; j++) {
	if (j>=JMAX)
	  throw std::runtime_error("hyperbolicAnomaly solution took too long");
	dh = meanAnomaly - e*sinh(h) + h;
	if (abs(dh) < tol)
	  break;  // Converged.
	h += dh / (e*cosh(h)-1);  // Newton correction.
      }
      /*Coordinates and velocity in the system aligned with orbit: */
      c = cosh(h);
      s = sinh(h);
      double b = -a * sqrt(e*e-1);
      r0[0] = a * (c - e);
      r0[1] = b * s;
      dt = (e*c-1) / meanMotion;  // This is dt/dHyperbolicAnomaly
      v0[0] = a * s / dt;
      v0[1] = b * c  / dt;
	  
    } // End hyperbolic/elliptical split     
    
    /* Rotate about z to put perihelion at arg of peri */
    c = cos(el[Elements::AOP]);  s = sin(el[Elements::AOP]);
    r1[0] = r0[0]*c - r0[1]*s;
    r1[1] = r0[0]*s + r0[1]*c;
    v1[0] = v0[0]*c - v0[1]*s;
    v1[1] = v0[0]*s + v0[1]*c;

    /* Rotate about x axis to incline orbit */
    c = cos(el[Elements::I]);  s = sin(el[Elements::I]);
    r2[0] = r1[0];
    r2[1] = r1[1]*c;
    r2[2] = r1[1]*s;
    v2[0] = v1[0];
    v2[1] = v1[1]*c;
    v2[2] = v1[1]*s;

    /* Rotate about z axis to bring node to longitude */
    c = cos(el[Elements::LAN]);  s = sin(el[Elements::LAN]);

    CartesianEcliptic x(r2[0]*c - r2[1]*s, r2[0]*s + r2[1]*c, r2[2]);
    CartesianEcliptic v(v2[0]*c - v2[1]*s, v2[0]*s + v2[1]*c, v2[2]);

    State out;
    out.x = CartesianICRS(x);
    out.v = CartesianICRS(v);
    out.tdb = tdb;

    if (heliocentric) {
      if (!ephem)
	throw runtime_error("Need an Ephemeris for heliocentric getState");
      auto sun = ephem->state(orbits::SUN, tdb);
      out.x += sun.x.getVector();
      out.v += sun.v.getVector();
    }
  
    return out;
  }

  Elements
  heliocentric2barycentric(const Elements& el, double tdb, const Ephemeris& ephem) {
    State s = getState(el, tdb, true, &ephem);
    return getElements(s);
  }

  Elements
  barycentric2heliocentric(const Elements& el, double tdb, const Ephemeris& ephem) {
    State s = getState(el, tdb);
    return getElements(s, true, &ephem);
  }



  // For calculation of orbit derivatives wrt state vector, define these units:
  class ScalarDerivative {
  public:
    typedef linalg::SMatrix<double,1,6> Deriv;
    EIGEN_NEW

    ScalarDerivative(double s_, const Deriv& ds_): s(s_), ds(ds_) {}
    ScalarDerivative(): s(0.), ds(0.) {}

    // A scalar and its derivative wrt state
    double s;
    Deriv ds;  // Row vector

    // Unaries:
    ScalarDerivative operator-() const {
      return ScalarDerivative(-s, -ds);
    }
    ScalarDerivative inverse() const {
      return ScalarDerivative(1./s, -ds / (s*s));
    }

    // Binaries with scalar:
    ScalarDerivative operator*(double rhs) const {
      return ScalarDerivative(s*rhs, ds*rhs);
    }
    ScalarDerivative operator/(double rhs) const {
      return ScalarDerivative(s/rhs, ds/rhs);
    }
    ScalarDerivative operator+(double rhs) const {
      return ScalarDerivative(s+rhs, ds);
    }
    ScalarDerivative operator-(double rhs) const {
      return ScalarDerivative(s-rhs, ds);
    }
    // Binaries with ScalarDerivative:
    ScalarDerivative operator*(const ScalarDerivative& rhs) const {
      return ScalarDerivative(s*rhs.s, ds*rhs.s + s*rhs.ds);
    }
    ScalarDerivative operator/(const ScalarDerivative& rhs) const {
      return ScalarDerivative(s/rhs.s, (ds*rhs.s - s*rhs.ds)/(rhs.s*rhs.s));
    }
    ScalarDerivative operator+(const ScalarDerivative& rhs) const {
      return ScalarDerivative(s+rhs.s, ds + rhs.ds);
    }
    ScalarDerivative operator-(const ScalarDerivative& rhs) const {
      return ScalarDerivative(s-rhs.s, ds - rhs.ds);
    }
  };

  class VectorDerivative {
  public:
    // A vector and its derivative wrt state
    typedef     linalg::SMatrix<double,3,6> Deriv;
    EIGEN_NEW
    
    Vector3 v;
    Deriv dv;

    VectorDerivative(const Vector3& v_, const Deriv& dv_): v(v_), dv(dv_) {}
    VectorDerivative(): v(0.), dv(0.) {}
    
    ScalarDerivative operator[](int i) const {
      // Get one component of vector
      return ScalarDerivative(v[i], dv.row(i));
    }

    // Unaries:
    VectorDerivative operator-() const {
      return VectorDerivative(-v, -dv);
    }
    ScalarDerivative magnitude() const {
      // Returns |v| and its derivative
      ScalarDerivative out;
      out.s = sqrt(v.dot(v));
      out.ds = v.transpose() * dv / out.s;
      return out;
    }

    // Binary w/scalar
    VectorDerivative operator*(double rhs) const {
      return VectorDerivative(v*rhs, dv*rhs);
    }
    VectorDerivative operator/(double rhs) const {
      return VectorDerivative(v/rhs, dv/rhs);
    }

    // Binary w/ScalarDerivative
    VectorDerivative operator*(const ScalarDerivative& rhs) const {
      return VectorDerivative(v * rhs.s,
			      dv * rhs.s + v * rhs.ds);  // Last term is outer product
    }

    // Binary w/VectorDerivative
    VectorDerivative operator+(const VectorDerivative& rhs) const {
      return VectorDerivative(v + rhs.v, dv + rhs.dv);
    }
    VectorDerivative operator-(const VectorDerivative& rhs) const {
      return VectorDerivative(v - rhs.v, dv - rhs.dv);
    }
    ScalarDerivative dot(const VectorDerivative& rhs) const {
      return ScalarDerivative(v.transpose() * rhs.v,
			      v.transpose() * rhs.dv + rhs.v.transpose() * dv);
    }

    VectorDerivative cross(const VectorDerivative& rhs) const {
    VectorDerivative out;
    out.v = v.cross(rhs.v);
    for (int i=0; i<6; i++)
      out.dv.col(i) = v.cross(rhs.dv.col(i)) + dv.col(i).cross(rhs.v);
    return out;
    }

  };

  ScalarDerivative
  atan2(const ScalarDerivative& y, const ScalarDerivative& x) {
    return ScalarDerivative( std::atan2(y.s,x.s),
			     (x.s * y.ds - x.ds * y.s) / (x.s*x.s+y.s*y.s));
  }
  ScalarDerivative
  atanh(const ScalarDerivative& x) {
    return ScalarDerivative( std::atanh(x.s),
			     x.ds / (1-x.s*x.s));
  }
  ScalarDerivative
  acos(const ScalarDerivative& x) {
    return ScalarDerivative( std::acos(x.s),
			     x.ds * (-pow(1-x.s*x.s,-0.5)));
  }
  ScalarDerivative
  sin(const ScalarDerivative& x) {
    return ScalarDerivative( std::sin(x.s),
			     x.ds * std::cos(x.s));
  }
  ScalarDerivative
  sinh(const ScalarDerivative& x) {
    return ScalarDerivative( std::sinh(x.s),
			     x.ds * std::cosh(x.s));
  }
  ScalarDerivative
  pow(const ScalarDerivative& x, double p) {
    return ScalarDerivative( std::pow(x.s,p),
			     x.ds * (p * std::pow(x.s, p-1.)));
  }
  
  Matrix66
  getElementDerivatives(const State& s,
			bool heliocentric=false, const Ephemeris* ephem=nullptr) {

    double mass = heliocentric ? GM : SolarSystemGM;
    double tdb = s.tdb;

    // Convert state to ecliptic
    Vector3 R, V;
    VectorDerivative x;
    VectorDerivative v;
    CartesianEcliptic ecliptic;
    Matrix33 icrs2ecliptic;

    if (heliocentric) {
      // Remove solar motion from state
      if (!ephem)
	throw runtime_error("Need an Ephemeris for heliocentric getElements");
      auto sun = ephem->state(orbits::SUN, tdb);
      auto net = s;
      net.x -= sun.x;
      net.v -= sun.v;

      // Rotate into ecliptic system
      ecliptic.convertFrom(net.x, icrs2ecliptic);
      x.v = ecliptic.getVector();
      x.dv.subMatrix(0,3,0,3) = icrs2ecliptic;

      ecliptic.convertFrom(net.v, icrs2ecliptic);
      v.v = ecliptic.getVector();
      v.dv.subMatrix(0,3,3,6) = icrs2ecliptic;
    } else {
      // Rotate into ecliptic system
      ecliptic.convertFrom(s.x, icrs2ecliptic);
      x.v = ecliptic.getVector();
      x.dv.subMatrix(0,3,0,3) = icrs2ecliptic;

      ecliptic.convertFrom(s.v, icrs2ecliptic);
      v.v = ecliptic.getVector();
      v.dv.subMatrix(0,3,3,6) = icrs2ecliptic;
    }

    Matrix66 out;
#ifdef CHECK
    Elements elements;
#endif

    // The z unit vector:
    VectorDerivative z;
    z.v.setZero();
    z.v[2] = 1.;
    z.dv.setZero();

    auto r = x.magnitude();
    auto rInverse = r.inverse();

    // Angular momentum vector:
    auto h = x.cross(v);

    // Eccentricity vector:
    auto evec = v.cross(h) / mass - x * rInverse;

    // Ascending node vector:
    auto ascendingNode = z.cross(h);

    // Semi-latus rectum
    auto p = h.dot(h) / mass;

    // Ready to start assigning some orbital elements
    // Semi-major axis
    ScalarDerivative stmp =  rInverse * 2. - v.dot(v) / mass;
    auto a = stmp.inverse();
    out.row(Elements::A) = a.ds;
#ifdef CHECK
    elements[Elements::A] = a.s;
#endif

    // Eccentricity
    auto e = evec.magnitude();
    out.row(Elements::E) = e.ds;
#ifdef CHECK
    elements[Elements::E] = e.s;
#endif

    // Inclination:
    stmp = h[2] * h.magnitude().inverse();
    auto inc = acos(stmp);
    out.row(Elements::I) = inc.ds;
#ifdef CHECK
    elements[Elements::I] = inc.s;
#endif

    // Ascending node:
    stmp = atan2( ascendingNode[1], ascendingNode[0]);
    out.row(Elements::LAN) = stmp.ds;
#ifdef CHECK
    elements[Elements::LAN] = stmp.s;
#endif

    // Argument of Perihelion:
    stmp = ascendingNode.dot(evec) / (ascendingNode.magnitude() * e);
    ScalarDerivative aop = evec.v[2]<0. ? -acos(stmp) : acos(stmp);
    out.row(Elements::AOP) = aop.ds;
#ifdef CHECK
    elements(Elements::AOP) = aop.s;
#endif

    // Now the manipulations to get mean anomaly
    auto xBar = (p-r)/e;
    auto yBar = x.dot(v) * pow(p/mass,0.5) / e;
    if (e.s < 1.) {
      // Elliptical solution
      auto cE = xBar + a*e;
      auto sE = yBar * pow(-evec.dot(evec)+1.,-0.5);
      auto eccentricAnomaly = atan2(sE, cE);
      auto meanAnomaly = eccentricAnomaly - e * sin(eccentricAnomaly);

      // And the final element, time of perhelion passage:
      auto top = -meanAnomaly * pow(a,1.5) / sqrt(mass) + tdb;
      out.row(Elements::TOP) = top.ds;
    } else {
      // Hyperbolic solution
      stmp = pow((e-1.) / (e+1.),0.5) * yBar / (x.magnitude()+xBar);
      auto hyperbolicAnomaly = atanh(stmp) * 2.;
      auto meanAnomaly = e * sinh(hyperbolicAnomaly) - hyperbolicAnomaly;
      // And the final element, time of perhelion passage:
      auto top = -meanAnomaly * pow(-a,1.5) / sqrt(mass) + tdb;
      out.row(Elements::TOP) = top.ds;
    }      
    return out;
  }
#ifdef CHECK
  elements(Elements::TOP) = top.s;

  /**/cerr << "New vs old:" << endl;
  auto oldel = getElements(s,heliocentric);
  for (int i=0; i<6; i++)
    cerr << i << " " << elements[i] << " -- " << oldel[i] << endl;
#endif

} // end namespace
