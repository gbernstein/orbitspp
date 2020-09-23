// Classes to use for geometric calculations

#ifndef PLANEGEOMETRY_H
#define PLANEGEOMETRY_H

#include <vector>
#include "LinearAlgebra.h"
#include "OrbitTypes.h"

namespace orbits {

  class Point: public Vector2 {
  public:
    Point(double x=0., double y=0.) {
      (*this)[0] = x; (*this)[1] = y;
    }
    Point(const Vector2& rhs): Vector2(rhs) {}
    
    double distanceSq(const Point& rhs) const {
      return (*this - rhs).squaredNorm();
    }
    DVector distanceSq(const DMatrix& rhs) const {
      // distance-sq to a set of points in Nx2 matrix
      DMatrix dx = rhs.rowwise() - this->transpose();
      return dx.rowwise().squaredNorm();
    }
  };

  class Segment {
    // Class describing a directed line segment that goes from p1 to p2
  public:
  EIGEN_NEW
  Segment(const Point& p1_, const Point& p2_): p1(p1_), p2(p2_) {
      Vector2 slope = p2-p1;
      slope /= slope.norm();
      perpUnit[0] = slope[1];
      perpUnit[1] = -slope[0];
    }
    // Perp distance to point from line (not just segment) (signed, pos to right)
    double perpDistance(const Point& p) const {
      return (p-p1).dot(perpUnit);
    }
    // Distance to a list of points (Nx2 matrix)
    DVector perpDistance(const DMatrix& p) const {
      DMatrix dx = p.rowwise() - p1.transpose();
      return dx * perpUnit; //dx.col(0)*perpUnit[0]+dx.col(1)*perpUnit[1];
    }

    // Is the point to the right of the line?
    bool onRight(const Point& p) const {return perpDistance(p)>0.;}  
    // same for Nx2 array of points
    BVector onRight(const DMatrix& pts) const {return perpDistance(pts).array()>0.;}
    
    // Return minimum (squared) distance from any point on Segment:
    double minDistanceSq(const Point& p) const;
    Point p1;
    Point p2;
  private:
    Vector2 perpUnit; // Unit vector perpendicular to ray (to right)
  };

  class ConvexPolygon {
    // Convex polygon defined by clockwise list of vertices.  There is no automatic check
    // for convexity or clockwise.
    friend class AffineTransformation;
  public:
    ConvexPolygon() = default;
    ConvexPolygon(const std::vector<Point>& vertices);

    // Force check of convexity - each edge's endpoint must be to
    // right of previous edge (or colinear)
    bool isConvex() const;
    double area() const;
    
    // Check whether one or many points are interior to polygon
    bool inside(const Point& p) const;
    BVector inside(const DMatrix& pts) const;

    // Is the point within requested distance of interior?
    bool isWithin(const Point& p, double distance) const;
    
    int size() const {return edges.size();}
    Segment edge(int i) const {return edges[i];}

    EIGEN_NEW
  private:
    vector<Segment, Eigen::aligned_allocator<Segment>> edges;
  };
    
  class Circle {
    friend class AffineTransformation;
    friend class Ellipse;
  public:
    Circle(const Point& ctr_, double radius_): ctr(ctr_), radius(radius_), rsq(radius_*radius_) {}
    Point ctr;
    double radius;
    
    // Is point interior to circle?
    bool inside(const Point& p) const {
      return ctr.distanceSq(p) < rsq;
    }
    // booleans of whether Nx2 point array is interior
    BVector inside(const DMatrix& p) const {
      return ctr.distanceSq(p).array() < rsq;
    }

    // Whether circles share interior
    bool intersects(const Circle& rhs) const {
      return ctr.distanceSq(rhs.ctr) < (radius+rhs.radius)*(radius+rhs.radius);
    }
    
    // Whether any finite part of segment is inside the circle
    bool intersects(const Segment& rhs) const {
      return rhs.minDistanceSq(ctr) < rsq;
    }
    
    // True if finite area in common.
    bool intersects(const ConvexPolygon& rhs) const {
      return rhs.isWithin(ctr,radius);
    }
    
    // ?? Rest of Circle not implemented yet
    // Intersection tests that return locations of intersection.
    // in which case determine intersection points ?? always 2?
    bool intersects(const Segment& rhs, Point& p1, Point& p2) const;
    bool intersects(const Circle& rhs, Point& p1, Point& p2) const;
    
    // Return area of intersection with polygon.
    double intersectionArea(const ConvexPolygon& rhs) const;

    EIGEN_NEW

  private:
    double rsq;
  };

  class Ellipse;  // forward declaration
  
  class AffineTransformation {
    // Map of the form (x,y) -> M * (x-x0, y-y0)
  public:
    AffineTransformation(const Point& xy0_, const Matrix22& m_): xy0(xy0_), m(m_) {}
    // Operate on single or  set of points as Nx2 matrix
    Point operator()(const Point& p) const {
      return Vector2(m * (p-xy0));
    }
    DMatrix operator()(const DMatrix& pts) const {
      return (pts.rowwise()-xy0.transpose()) * m.transpose();
    }

    // Operate on other shapes
    Segment operator()(const Segment& s) const;
    ConvexPolygon operator()(const ConvexPolygon& p) const;
    Ellipse operator()(const Circle& c) const;
    Ellipse operator()(const Ellipse& e) const;
    EIGEN_NEW
  private:
    Point xy0;
    Matrix22 m;
  };
  
  class Ellipse {
  public:
    // Class representing elliptical region
    Ellipse(const Point& ctr_, const Matrix22& cov_): ctr(ctr_),
						      cov(cov_),
						      circularizer(0.) {setup();}
    Ellipse(const Circle& c): ctr(c.ctr), cov(0.) {
      cov(0,0) = cov(1,1) = c.rsq;
    }
    Point ctr;
    Matrix22 cov;
    
    double area() const {return sqrt(PI*a*b);}

    // Is point interior to ellipse?
    bool inside(const Point& rhs) const {return chisq(rhs)<1.;}
    // for Nx2 array of points:
    BVector inside (const DMatrix& rhs) const {return chisq(rhs).array()<1.;}

    // Return quadratic form for ellipse:
    double chisq(const Point& rhs) const {
      return (rhs-ctr).transpose() * invCov * (rhs-ctr);
    }
    // Return quadratic form for Nx2 array
    DVector chisq(const DMatrix& rhs) const {
      DMatrix dx = rhs.rowwise() - ctr.transpose();
      DMatrix tmp = dx.array() * (dx * invCov).array();
      return tmp.rowwise().sum();
    }

    int roughIntersects(const Circle& rhs) const;
    int roughIntersects(const Ellipse& rhs) const;
    // Return -1, 0, 1 as circle definitely doesn't, might, or definitely does
    // intersect the circle, based on major/minor axes.

    void getSemiMajorMinor(double& a, double& b) const ;

    // Returns a transformation that renders Ellipse into unit circle.
    AffineTransformation getCircularizer() const {
      return AffineTransformation(ctr,circularizer);
    }

    bool intersects(const ConvexPolygon& poly) const;
    // True if finite area in common.
    double intersectionArea(const ConvexPolygon& rhs) const;
    // Return area of intersection with polygon. ?? not implemented yet.
    EIGEN_NEW

  private:
    Matrix22 invCov;
    double a,b; // Major, minor axes.
    Matrix22 circularizer;
    void setup();
  };

} // end namespace

#endif  // PLANEGEOMETRY_H
