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

  class Ray {
    // Class describing a directed line that goes from p1 to p2
  public:
  Ray(const Point& p1_, const Point& p2_): p1(p1_), p2(p2_) {
      Vector2 slope = p2-p1;
      slope /= slope.norm();
      perpUnit[0] = slope[1];
      perpUnit[1] = -slope[0];
    }
    // Perp distance to point from line (signed, pos to right)
    double distance(const Point& p) const {
      return (p-p1).dot(perpUnit);
    }
    // Distance to a list of points (Nx2 matrix)
    DVector distance(const DMatrix& p) const {
      DMatrix dx = p.rowwise() - p1.transpose();
      return dx * perpUnit; //dx.col(0)*perpUnit[0]+dx.col(1)*perpUnit[1];
    }

    // Is the point to the right of the ray?
    bool onRight(const Point& p) const {return distance(p)>0.;}  
    // same for Nx2 array of points
    BVector onRight(const DMatrix& pts) const {return distance(pts).array()>0.;}
    
    Point p1;
    Point p2;
  private:
    Vector2 perpUnit; // Unit vector perpendicular to ray (to right)
    EIGEN_NEW
  };

  class ConvexPolygon {
    // Convex polygon defined by clockwise list of vertices.  There is no automatic check
    // for convexity or clockwise.
  public:
    ConvexPolygon(const std::vector<Point>& vertices) {
      for (int i=0; i<vertices.size()-1; i++) 
	edges.push_back(Ray(vertices[i],vertices[i+1]));
      // Close polygon
      edges.push_back(Ray(vertices.back(), vertices.front()));
    }
    // Force check of convexity - each edge's endpoint must be to
    // right of previous edge (or colinear)
    bool isConvex() const {
      for (int i=0; i<edges.size()-1; i++) 
	if (edges[i].distance(edges[i+1].p2)<0.)
	  return false;
      return edges.back().distance(edges.front().p2)>=0.;
    }
    double area() const {
      // Use shoelaces formulae, relative to first point to avoid roundoff
      double x0 = edges.front().p1[0];
      double y0 = edges.front().p1[1];
      double out=0.;
      for (auto& e : edges)
	out += (e.p1[0]-x0)*(e.p2[1]-y0)-(e.p1[1]-y0)*(e.p2[0]-x0);
      return -0.5*out;
    }
    
    // Check whether one or many points are interior to polygon
    bool inside(const Point& p) const {
      for (auto& e : edges)
	if (!e.onRight(p)) return false;
      return true;
    }
    BVector inside(const DMatrix& pts) const {
      BVector out = edges.front().onRight(pts);
      for (int i=1; i<edges.size(); i++)
	out.array() *= edges[i].onRight(pts).array();
      return out;
    }
    EIGEN_NEW
    // ??? Eigen vector business ??
  private:
    vector<Ray> edges;
  };
    
  class Circle {
  public:
    Circle(const Point& ctr_, double radius_): ctr(ctr_), radius(radius_), rsq(radius_*radius_) {}
    Point ctr;
    double radius;
    
    bool inside(const Point& p) const {
      // Is point interior to circle?
      return ctr.distanceSq(p) < rsq;
    }
    BVector inside(const DMatrix& p) const {
      // booleans of whether Nx2 point array is interior
      return ctr.distanceSq(p).array() < rsq;
    }
    bool inside(const Ray& rhs) const;
    // Whether line crosses inside the circle
    bool inside(const Ray& rhs, Point& p1, Point& p2) const;
    // Return true if line crosses inside the circle,
    // in which case determine intersection points

    bool intersects(const Circle& rhs) const;
    // Whether circles share interior
    bool crosses(const Ray& rhs, Point& p1, Point& p2) const;
    // Return true if circle boundaries cross.
    // In which case determine intersection points
    
    bool intersects(const ConvexPolygon& rhs) const;
    // True if finite area in common.
    double commonArea(const ConvexPolygon& rhs) const;
    // Return area of intersection with polygon.

  private:
    double rsq;
    EIGEN_NEW
  };

  class AffineTransformation; // Forward declaration
  
  class Ellipse {
    // Class representing elliptical region
    Ellipse(const Point& ctr_, const Matrix22& cov_): ctr(ctr_), cov(cov_) {setup();}
    Point ctr;
    Matrix22 cov;
    
    double area() const {return sqrt(detCov);}

    bool inside(const Point& rhs) const; // Is point interior to ellipse?
    BVector inside (const DMatrix& rhs) const; // Nx2 array of points

    double chisq(const Point& rhs) const; // Return quadratic form
    DVector chisq(const DMatrix& rhs) const; // Return quadratic form for Nx2 array

    int roughIntersects(const Circle& rhs) const;
    int roughIntersects(const Ellipse& rhs) const;
    // Return -1, 0, 1 as circle definitely doesn't, might, or definitely does
    // intersect the circle, based on major/minor axes.

    void getSemiMajorMinor(double& a, double& b) const ;
    AffineTransformation getCircularizer() const;
    // Returns a transformation that renders Ellipse into unit circle.

    bool intersect(const ConvexPolygon& rhs) const;
    // True if finite area in common.
    double commonArea(const ConvexPolygon& rhs) const;
    // Return area of intersection with polygon.

  private:
    void setup(); // Calculate various things...
    double detCov; // Determinant of covariance
    Matrix22 invCov;
    // ?? Eigenvalues etc
    
    EIGEN_NEW
  };

  class AffineTransformation {
    // Map of the form (x,y) -> M * (x-x0, y-y0)
  public:
    AffineTransformation(const Point& xy0_, const Matrix22& m_): xy0(xy0_), m(m_) {}
    // Operate on single or  set of points as Nx2 matrix
    Point operator()(const Point& p) const {
      Vector2 out = m * (p-xy0);
      return out;
    }
    DMatrix operator()(const DMatrix& pts) const {
      return (pts.rowwise()-xy0.transpose()) * m.transpose();
    }

    // Operate on other shapes
    Ray operator()(const Ray& r) const;
    Ellipse operator()(const Circle& c) const;
    Ellipse operator()(const Ellipse& c) const;
  private:
    Point xy0;
    Matrix22 m;
    EIGEN_NEW
  };
    
} // end namespace

#endif  // PLANEGEOMETRY_H
