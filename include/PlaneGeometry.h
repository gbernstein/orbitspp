// Classes to use for geometric calculations

#ifndef PLANEGEOMETRY_H
#define PLANEGEOMETRY_H

#include "LinearAlgebra.h"


namespace orbits {

  typedef linalg::SVector<double,2> Vector2;
  typedef linalg::SMatrix<double,2,2> Matrix22;
  typedef linalg::Matrix<double> DMatrix;
  typedef linalg::Vector<double> DVector;
  typedef linalg::Vector<bool> BVector;

  class Point: public Vector2 {
  public:
    Point(double x=0., double y=0.) {
      (*this)[0] = x; (*this)[1] = y;
    }
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
    Ray(const Point& p1_, const Point& p2_);
    bool onRight(const Point& p) const;  // Is the point to the right of the right?
    double distance(const Point& p) const; // Perp distance to point from line (signed???)
    DVector distance(const DMatrix& p) const; // Distance to a list of points (Nx2 matrix)
    
  private:
    Point p1;
    Point p2;
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
    
  private:
    double rsq;
  };

  class Ellipse {
    // Class representing elliptical region
    Ellipse(const Point& ctr_, const Matrix22& cov_): ctr(ctr_), cov(cov_) {}
    Point ctr;
    Matrix22 cov;
    
    double area() const;

    bool inside(const Point& rhs) const; // Is point interior to ellipse?
    BVector inside (const DMatrix& rhs) const; // Nx2 array of points

    double chisq(const Point& rhs) const; // Return quadratic form
    DVector chisq(const DMatrix& rhs) const; // Return quadratic form for Nx2 array

    int intersect(const Circle& rhs) const;
    int intersect(const Ellipse& rhs) const;
    // Return -1, 0, 1 as circle definitely doesn't, might, or definitely does
    // intersect the circle, based on major/minor axes.

    void getSemiMajorMinor(double& a, double& b) const ;
  };
    
} // end namespace

#endif  // PLANEGEOMETRY_H
