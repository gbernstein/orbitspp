// Numerical tests of analytic derivatives of conversions between
// ABG, state vectors, and orbital elements
#include "OrbitTypes.h"
#include "Elements.h"
#include "AstronomicalConstants.h"

using namespace orbits;

const double tdb0 = 1.5;  // Epoch for orbit

Vector6
rotateState(const Frame& f,
	    const ABG& abg) {
  Vector3 x;
  Vector3 v;
  abg.getState(0., x, v);
  Vector6 out;
  out.subVector(0,3) = f.toICRS(x);
  out.subVector(3,6) = f.toICRS(v, true);
  return out;
}

Elements
elementfunc(const Frame& f,
	    const ABG& abg) {
  Vector3 x;
  Vector3 v;
  abg.getState(0., x, v);
  State s;
  s.x = astrometry::CartesianICRS(f.toICRS(x));
  s.v = astrometry::CartesianICRS(f.toICRS(v,true));
  s.tdb = 0.;
  return getElements(s);
}


int main(int argc,
	 char *argv[])
{
  // Set up a reference frame and an ABG
  //  Coordinate origin somewhere 1 AU from Sun
  astrometry::CartesianICRS origin( 0.8, 0.3, -sqrt(1. - 0.8*0.8 - 0.3*0.3));
  //  Pick a direction
  astrometry::SphericalICRS axis(30*DEGREE, -8*DEGREE);
  astrometry::Orientation orient(axis);
  orient.alignToEcliptic();
  Frame f(origin, orient, tdb0);

  // An initial orbit:
  double gamma = 1. / 50.;
  double vesc = sqrt(2. * GM * pow(gamma,3.));
  ABG abg;
  abg[ABG::A] = -0.5*DEGREE;
  abg[ABG::B] = -0.1*DEGREE;
  abg[ABG::G] = gamma;
  abg[ABG::ADOT] = +0.5*vesc;
  abg[ABG::BDOT] = -0.1*vesc;
  abg[ABG::GDOT] = +0.2*vesc;

  // Define small steps for numerical derivatives
  ABG dabg;
  dabg[ABG::A] = 20*ARCSEC;
  dabg[ABG::B] = 20*ARCSEC;
  dabg[ABG::G] = 0.01 * gamma;
  dabg[ABG::ADOT] = 0.01*vesc;
  dabg[ABG::BDOT] = 0.01*vesc;
  dabg[ABG::GDOT] = 0.01*vesc;
  
  /**/ dabg /= 4.;
  // First check dS/dABG
  Matrix66 dSdABG_frame = abg.getStateDerivatives();
  /**/cerr << "dSdABG_frame: " << endl << dSdABG_frame << endl;

  Matrix66 numeric;
  for (int i=0; i<6; i++) {
    ABG tmp = abg;
    tmp[i] += dabg[i];
    astrometry::Vector3 xp,vp,xm,vm;
    tmp.getState(0.,xp,vp);
    tmp = abg;
    tmp[i] -= dabg[i];
    tmp.getState(0.,xm,vm);
    Vector6 ds(0.);
    ds.subVector(0,3) = (xp - xm) / (2*dabg[i]);
    ds.subVector(3,6) = (vp - vm) / (2*dabg[i]);
    numeric.col(i) = ds;
  }
  /**/cerr << "--Numeric: " << endl << numeric << endl;
  
 
  // Then include rotation to ICRS
  Matrix66 dSdABG;
  
  // Analytic; note that Frame takes Nx3 matrices so put xyz second
  DMatrix tmp = dSdABG_frame.subMatrix(0,3,0,6).transpose();
  dSdABG.subMatrix(0,3,0,6) = f.toICRS(tmp , true).transpose();
  tmp = dSdABG_frame.subMatrix(3,6,0,6).transpose();
  dSdABG.subMatrix(3,6,0,6) = f.toICRS(tmp, true).transpose();
  /**/cerr << "dSdABG: " << endl << dSdABG << endl;
  
  // Numeric version:
  Vector6 sp,sm;
  for (int i=0; i<6; i++) {
    ABG tmp = abg;
    tmp[i] += dabg[i];
    sp = rotateState(f,tmp);
    tmp[i] -= 2.*dabg[i];
    sm = rotateState(f,tmp);
    numeric.col(i) = (sp - sm) / (2*dabg[i]);
  }
  /**/cerr << "--Numeric: " << endl << numeric << endl;

  // Then include orbital element derivatives
  Vector3 x0;
  Vector3 v0;
  abg.getState(0.,x0,v0);
  State s;
  s.x = astrometry::CartesianICRS(f.toICRS(x0));
  s.v = astrometry::CartesianICRS(f.toICRS(v0,true));
  s.tdb = f.tdb0;
  /**/cerr << "Elements: " << getElements(s) << endl;
  /**/cerr << " or " << elementfunc(f,abg) << endl;
  // Get element derivatives
  Matrix66 dEdABG = getElementDerivatives(s) * dSdABG;
  /**/cerr << "dEdABG: " << endl << dEdABG << endl;
 
  // Numeric:

  Elements ep, em;
  for (int i=0; i<6; i++) {
    ABG tmp = abg;
    tmp[i] += dabg[i];
    ep = elementfunc(f,tmp);
    tmp[i] -= 2.*dabg[i];
    em = elementfunc(f,tmp);
    numeric.col(i) = (ep - em) / (2*dabg[i]);
  }
  /**/cerr << "--Numeric: " << endl << numeric << endl;

  // A cross check, using a= (2/r - v^2/GM)^{-1}
  Elements el = getElements(s);
  Vector3 x = s.x.getVector();
  double r = sqrt(x.dot(x));
  Vector6 dads;
  dads.subVector(0,3) = 2. * pow(el[Elements::A],2.) * pow(r,-3.) * s.x.getVector();
  dads.subVector(3,6) = (2. * pow(el[Elements::A],2.) / GM)* s.v.getVector();
  //**/cerr << "dads: " << endl << dads << endl;
  /**/cerr << "dadABG (should match first row above): "  << endl<< dads.transpose() * dSdABG << endl;
  exit(0);
}
