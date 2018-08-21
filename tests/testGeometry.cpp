// Test PlaneGeometry
#include "PlaneGeometry.h"
#include "AstronomicalConstants.h"
#include <iostream>

using namespace std;
using namespace orbits;

int
main(int argc,
     char *argv[]) {
  double a=12., b=2.;
  double PA=120.*DEGREE;

  Point ctr(-1., 0.);
  Matrix22 cov;
  cov(0,0) = 0.5* (a*a+b*b + (a*a-b*b)*cos(2*PA));
  cov(1,1) = 0.5* (a*a+b*b - (a*a-b*b)*cos(2*PA));
  cov(0,1) = cov(1,0) = 0.5*(a*a-b*b)*sin(2*PA);

  Ellipse e(ctr,cov);
  AffineTransformation aff = e.getCircularizer();
  
  // Walk around border of 2x ellipse
  for (double theta=0; theta<TPI; theta+=PI/10.) {
    double xx = a * cos(theta);
    double yy = b * sin(theta);
    Point p( xx*cos(PA) - yy*sin(PA), xx*sin(PA)+yy*cos(PA));
    p *= 2;
    p += ctr;
    Point u = aff(p);
    cout << p[0] << "," << p[1] << " " << e.chisq(p) << " " << u.norm() << endl;
  }

  cout << "Enter x0,y0,x1,y1 segments-> ";
  double x0,y0,x1,y1;
  while (cin >> x0 >> y0 >> x1 >> y1) 
    cout << Segment(Point(x0,y0),Point(x1,y1)).minDistanceSq(Point()) << "\n->";

  exit(0);
}
