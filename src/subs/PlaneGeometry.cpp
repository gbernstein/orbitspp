#include "PlaneGeometry.h"

using namespace orbits;

double
Segment::minDistanceSq(const Point& p) const {
  Point dp = Vector2(p-p1);
  // Point projection onto segment:
  double along = dp[1]*perpUnit[0]-dp[0]*perpUnit[1];
  if (along<=0.) {
    // closest approach at p1
    return dp.squaredNorm();
  } else if (along >= (p2-p1).norm()) {
    // closest approach at p2
    return (p-p2).squaredNorm();
  } else {
    return dp.dot(perpUnit) * dp.dot(perpUnit);
  }
}

ConvexPolygon::ConvexPolygon(const std::vector<Point>& vertices) {
  for (int i=0; i<vertices.size()-1; i++) 
    edges.push_back(Segment(vertices[i],vertices[i+1]));
  // Close polygon
  edges.push_back(Segment(vertices.back(), vertices.front()));
}

bool
ConvexPolygon::isConvex() const {
  for (int i=0; i<edges.size()-1; i++) 
    if (edges[i].perpDistance(edges[i+1].p2)<0.)
      return false;
  return edges.back().perpDistance(edges.front().p2)>=0.;
}

double
ConvexPolygon::area() const {
  // Use shoelaces formulae, relative to first point to avoid roundoff
  double x0 = edges.front().p1[0];
  double y0 = edges.front().p1[1];
  double out=0.;
  for (auto& e : edges)
    out += (e.p1[0]-x0)*(e.p2[1]-y0)-(e.p1[1]-y0)*(e.p2[0]-x0);
  return -0.5*out;
}


bool
ConvexPolygon::inside(const Point& p) const {
  // Check whether one or many points are interior to polygon
  for (auto& e : edges)
    if (!e.onRight(p)) return false;
  return true;
}
BVector
ConvexPolygon::inside(const DMatrix& pts) const {
  BVector out = edges.front().onRight(pts);
  for (int i=1; i<edges.size(); i++)
    out.array() *= edges[i].onRight(pts).array();
  return out;
}

// Is the point within requested distance of interior?
bool
ConvexPolygon::isWithin(const Point& p, double distance) const {
  bool interior = true;
  double dsq = distance*distance;
  for (auto& edge : edges) {
    if (edge.minDistanceSq(p)<=dsq) return true;
    interior &= edge.onRight(p);
  }
  return interior;
}

// Affine transformations
Segment
AffineTransformation::operator()(const Segment& s) const {
  return Segment( (*this)(s.p1), (*this)(s.p2));
}
ConvexPolygon
AffineTransformation::operator()(const ConvexPolygon& p) const {
  ConvexPolygon out;
  for (auto& edge : p.edges)
    out.edges.push_back( (*this)(edge) );
  return out;
}

Ellipse
AffineTransformation::operator()(const Circle& c) const {
  Matrix22 ccov(0.);
  ccov(0,0) = ccov(1,1) = c.rsq;  // Diagonal covariance matrix
  Matrix22 tcov = m.transpose() * ccov * m;
  return Ellipse( (*this)(c.ctr), tcov);
}

Ellipse
AffineTransformation::operator()(const Ellipse& e) const {
  Matrix22 tcov = m.transpose() * e.cov * m;
  return Ellipse( (*this)(e.ctr), tcov);
}

void
Ellipse::setup() {
  double invdet = 1./(cov(0,0)*cov(1,1)-cov(0,1)*cov(1,0));
  invCov(0,0) = cov(1,1)*invdet;
  invCov(1,1) = cov(0,0)*invdet;
  invCov(0,1) = invCov(1,0) = -cov(0,1)*invdet;
  double t = 0.5*(cov(0,0) + cov(1,1));  // Half of trace
  double e1 = 0.5*(cov(0,0) - cov(1,1));
  double e2 = cov(1,0);
  double e = hypot(e1,e2);
  if (e==0.) {
    a = b = sqrt(t);
    circularizer(1,1) = circularizer(0,0) = 1./a;
  } else {
    a = sqrt(t+e);
    b = sqrt(t-e);
    double c = sqrt(0.5*(1.+e1/e));
    double s = copysign(sqrt(0.5*(1-e1/e)),cov(1,1));
    circularizer(0,0) = c/a;
    circularizer(0,1) = -s/a;
    circularizer(1,0) = s/b;
    circularizer(1,1) = c/b;
  }
}

bool
Ellipse::intersects(const ConvexPolygon& poly) const {
  // Transform polygon and ellipse into unit circles:
  return getCircularizer()(poly).isWithin(Point(0.,0.),1.);
}

int
Ellipse::roughIntersects(const Ellipse& rhs) const {
  double dsq = ctr.distanceSq(rhs.ctr);
  if (dsq < (b + rhs.b)*(b + rhs.b)) {
    // Definite intersection
    return 1;
  } else if (dsq > (a + rhs.a)*(a + rhs.a)) {
    // Definite miss
    return -1;
  } else {
    // maybe
    return 0;
  }
}

int
Ellipse::roughIntersects(const Circle& rhs) const {
  double dsq = ctr.distanceSq(rhs.ctr);
  if (dsq < (b + rhs.radius)*(b + rhs.radius)) {
    // Definite intersection
    return 1;
  } else if (dsq > (a + rhs.radius)*(a + rhs.radius)) {
    // Definite miss
    return -1;
  } else {
    // maybe
    return 0;
  }
}
