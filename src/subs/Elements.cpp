// Tranformation between state vector and orbital elements.
// Not good for unbound orbits!

#include "OrbitTypes.h"
#include "AstronomicalConstants.h"
#include "Astrometry.h"

using namespace astrometry;

namespace orbits {

Elements
getElements(const State& s, bool heliocentric) {
  // Return orbital elements for state vector.
  // Use solar mass if heliocentric=true, else full solar system mass
  // The incoming state is assumed ICRS so we will convert to ecliptic
  // to get ecliptic J2000 elements.
  
  double mass = heliocentric ? GM : SolarSystemGM;
  double tdb = s.tdb;

  Vector3 R = CartesianEcliptic(s.x).getVector();
  Vector3 V = CartesianEcliptic(s.v).getVector();
  
  double rMagnitude, velSquare, radVelDotProduct;
  Vector3 eccentricityVector, angularMomentum, ascendingNode;
  double semimajor, eccentricity, inclination, longitudeOfAscendingNode;
  double semiLatusRectum;
  double hMagnitude, ascendingNodeMagnitude; /* magnitude of angular momentum */
  double ascEccDotProduct, argumentOfPerifocus;
  double xBar, yBar;
  double cosE, sinE, E1, E2, eccentricAnomaly;
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

  /* From here, we assume that the motion is elliptical */

  cosE = (xBar/semimajor) + eccentricity;
  sinE = yBar/(semimajor*sqrt(1-eccentricity*eccentricity));
  /* where semimajor*sqrt(1-eccentricity*eccentricity) is semiminor */
  eccentricAnomaly = atan2(sinE,cosE);

  meanAnomaly = eccentricAnomaly - eccentricity*sinE; /* radians */
  meanMotion = sqrt(mass/(pow(semimajor,3.)));
  timeOfPerifocalPassage = tdb - meanAnomaly/meanMotion;
  /* This comes from M=n(t-T) where t is epoch time and T is time of perifocal
     passage, in days */
				
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
getState(const Elements& el, double tdb, bool heliocentric) {
  Vector3  r0, v0, r1, v1, r2, v2;
  double eccentricAnomaly;
  double meanAnomaly;
  double c, s, dt;

  double mass = heliocentric ? GM : SolarSystemGM;
  double e = el[Elements::E];
  double a = el[Elements::A];
  
  if (e >= 1. || el[Elements::A] <=0.)
    throw std::runtime_error("elements_to_xv only for closed orbits now");

  /* get the eccentric Anomaly from the mean anomaly */
  meanAnomaly = (tdb - el[Elements::TOP]) * pow(a,-1.5) * sqrt(mass);
  /* Put it near 0 */
  meanAnomaly -= floor(meanAnomaly / TPI) * TPI;
  /* Use bisection to find solution (adapt Numerical Recipes)*/  
  {
    const int JMAX=40;
    const double TOLERANCE = 0.1*SECOND;
    double f, fmid, x1, x2, dx, xmid, rtb;
    int j;
    x1 = meanAnomaly - e;
    x2 = meanAnomaly + e;
    f   = x1 - e * sin(x1) - meanAnomaly;
    fmid= x2 - e * sin(x2) - meanAnomaly;
    if (f*fmid > 0.0) 
      FormatAndThrow<std::runtime_error>() << "Error, eccentricAnomaly root not bracketed; f,fmid = "
					   << f << " " << fmid;

    rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
    for (j=1;j<=JMAX;j++) {
      xmid=rtb+(dx *= 0.5);
      fmid= xmid - e * sin(xmid) - meanAnomaly;
      if (fmid <= 0.0) rtb=xmid;
      if (fabs(dx) < TOLERANCE || fmid == 0.0) break;
    }
    if (j>=JMAX)
      throw std::runtime_error("eccentricAnomaly took too long");
    meanAnomaly = rtb;
  }

  /*Coordinates and velocity in the system aligned with orbit: */
  c = cos(meanAnomaly);
  s = sin(meanAnomaly);
  r0[0] = a * (c - e);
  r0[1] = a * s * sqrt(1-e*e);
  dt = sqrt(pow(a,3.)/mass) * ( 1 - e*c);
  v0[0] = -a * s / dt;
  v0[1] = a * c * sqrt(1-e*e) / dt;

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
  
  return out;
}
} // end namespace
