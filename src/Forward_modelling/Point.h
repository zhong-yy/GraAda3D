
#ifndef _POINT_H
#define _POINT_H

// C++ includes
#include <cassert>
#include <cmath>
#include <iostream>
// Point operator*(double a, const Point& p);
// Point cross(const Point& v1, const Point& v2);
// Point unitCross(const Point& v1, const Point& v2);
/**
 * @brief 3D point
 *
 */
class Point {
 public:
  Point(const double x = 0., const double y = 0., const double z = 0.);
  Point(const Point& p);
  virtual ~Point();

  friend std::ostream& operator<<(std::ostream& os, const Point& point);

  // getter and setter
  double* get_xyz() { return _coords; }
  void set_xyz(double x, double y, double z);

  /**
   * @brief get x coordinate
   */
  double x() const { return _coords[0]; }

  /**
   * @brief get y coordinate
   */
  double y() const { return _coords[1]; }

  /**
   * @brief get z coordinate
   */
  double z() const { return _coords[2]; }
  double& x() { return _coords[0]; }
  double& y() { return _coords[1]; }
  double& z() { return _coords[2]; }

  // Overload operators.
  Point& operator=(const Point& p);
  double operator()(const unsigned int i) const;
  double& operator()(const unsigned int i);
  Point operator+(const Point& v) const;
  Point operator-(const Point& v) const;
  Point operator*(const double a) const;
  friend Point operator*(double a, const Point& p);
  Point operator/(double a) const;  // divided by a number
  double operator*(const Point& v) const;
  bool operator==(const Point& v) const;
  bool operator<(const Point& v) const;
  void operator=(const double a);

  // compound assignment
  Point& operator+=(const Point& v);
  Point& operator*=(double a);
  Point& operator/=(double a);

  Point operator-() const;
  void reverse();
  // Math functions
  // cross product
  friend Point cross(const Point& v1, const Point& v2);
  friend Point unitCross(const Point& v1, const Point& v2);
  // Point cross(const Point& v) const;

  Point unit() const;
  Point getUnit() const;  // get the unit vector
  void setUnit();

  double size() const;
  void zero();

  /**
   * @brief get the projection onto the plane defined by p1,p2 and p3
   * @param p1 
   * @param p2 
   * @param p3 
   * @return Point 
   */
  Point projection(const Point& p1, const Point& p2, const Point& p3) const;

  /**
   * @brief get the projection point on to the straight line passing through points p1 and p2
   * 
   * @param p1 
   * @param p2 
   * @return Point 
   */
  Point projection_line(const Point& p1, const Point& p2) const;

  double distance(const Point& p1) const;

 protected:
  double _coords[3];
  const double TOL;
};

#endif  // _POINT_H
