#include "Point.h"
Point::Point(const double x, const double y, const double z) : TOL(1.0e-10)
{
  _coords[0] = x;
  _coords[1] = y;
  _coords[2] = z;
}

Point::Point(const Point &p) : TOL(1.0e-10)
{
  for (unsigned int i = 0; i < 3; i++)
    _coords[i] = p._coords[i];
}

Point::~Point()
{
  // no space is allocated by the new operator.
  // so,there is no delete [].
}

Point &Point::operator=(const Point &p)
{
  _coords[0] = p._coords[0];
  _coords[1] = p._coords[1];
  _coords[2] = p._coords[2];

  return (*this);
}

Point Point::operator-() const
{
  return Point(-_coords[0], -_coords[1], -_coords[2]);
}

void Point::reverse()
{
  _coords[0] = -_coords[0];
  _coords[1] = -_coords[1];
  _coords[2] = -_coords[2];
}

void Point::set_xyz(double x, double y, double z)
{
  this->_coords[0] = x;
  this->_coords[1] = y;
  this->_coords[2] = z;
}

double Point::operator()(const unsigned int i) const
{
  assert(i < 3);
  return _coords[i];
}

double &Point::operator()(const unsigned int i)
{
  assert(i < 3);
  return _coords[i];
}

Point Point::operator+(const Point &p) const
{
  return Point(_coords[0] + p._coords[0], _coords[1] + p._coords[1],
               _coords[2] + p._coords[2]);
}

Point Point::operator-(const Point &p) const
{
  return Point(_coords[0] - p._coords[0], _coords[1] - p._coords[1],
               _coords[2] - p._coords[2]);
}

Point Point::operator*(const double factor) const
{
  return Point(_coords[0] * factor, _coords[1] * factor, _coords[2] * factor);
}

double Point::operator*(const Point &p) const
{
  return (_coords[0] * p(0) + _coords[1] * p(1) + _coords[2] * p(2));
}

Point Point::operator/(double a) const
{
  return Point(_coords[0] / a, _coords[1] / a, _coords[2] / a);
}

Point &Point::operator+=(const Point &v)
{
  this->_coords[0] += v._coords[0];
  this->_coords[1] += v._coords[1];
  this->_coords[2] += v._coords[2];
  return *this;
}

Point &Point::operator*=(double a)
{
  this->_coords[0] *= a;
  this->_coords[1] *= a;
  this->_coords[2] *= a;
  return *this;
}
Point &Point::operator/=(double a)
{
  this->_coords[0] = this->_coords[0] / a;
  this->_coords[1] = this->_coords[1] / a;
  this->_coords[2] = this->_coords[2] / a;
  return *this;
}

bool Point::operator==(const Point &rhs) const
{
  return ((std::abs(_coords[0] - rhs._coords[0]) +
           std::abs(_coords[1] - rhs._coords[1]) +
           std::abs(_coords[2] - rhs._coords[2])) < 3 * TOL);
}

bool Point::operator<(const Point &rhs) const
{
  // First we assume (this)<rhs true
  if (*this == rhs)
    return false;
  if ((*this)(0) < rhs(0))
    return true; //  <
  else if ((*this)(0) > rhs(0))
    return false; //  >
  else if (std::abs((*this)(0) - rhs(0)) < TOL)
  { // vx=rhsx
    if ((*this)(1) < rhs(1))
      return true;
    else if ((*this)(1) > rhs(1))
      return false;
    else if (std::abs((*this)(1) - rhs(1)) < TOL)
    { // vy=rhsy
      if ((*this)(2) < rhs(2))
        return true;
      else if ((*this)(2) > rhs(2))
        return false;
    }
  }
  return false;
}

void Point::operator=(const double a)
{
  _coords[0] = a;
  _coords[1] = a;
  _coords[2] = a;
}

// Point Point::cross(const Point& p) const {
//   return Point(_coords[1] * p._coords[2] - _coords[2] * p._coords[1],
//                -_coords[0] * p._coords[2] + _coords[2] * p._coords[0],
//                _coords[0] * p._coords[1] - _coords[1] * p._coords[0]);
// }

// v1 cross v2
Point cross(const Point &v1, const Point &v2)
{
  double i, j, k;
  double x1 = v1.x(), x2 = v2.x();
  double y1 = v1.y(), y2 = v2.y();
  double z1 = v1.z(), z2 = v2.z();
  i = y1 * z2 - y2 * z1;
  j = x2 * z1 - x1 * z2;
  k = x1 * y2 - x2 * y1;
  return Point(i, j, k);
}

// return unit vector of cross product
Point unitCross(const Point &v1, const Point &v2)
{
  double i, j, k;
  double x1 = v1.x(), x2 = v2.x();
  double y1 = v1.y(), y2 = v2.y();
  double z1 = v1.z(), z2 = v2.z();
  i = y1 * z2 - y2 * z1;
  j = x2 * z1 - x1 * z2;
  k = x1 * y2 - x2 * y1;
  Point cro(i, j, k);
  return cro.unit();
}

Point Point::unit() const
{
  const double length = size();
  return Point(_coords[0] / length, _coords[1] / length, _coords[2] / length);
}
Point Point::getUnit() const
{
  double value = sqrt(_coords[0] * _coords[0] + _coords[1] * _coords[1] +
                      _coords[2] * _coords[2]);
  return Point(_coords[0] / value, _coords[1] / value, _coords[2] / value);
}
void Point::setUnit()
{
  double length = size();
  for (int i = 0; i < 3; i++)
  {
    _coords[i] /= length;
  }
  return;
}

double Point::size() const
{
  double value = _coords[0] * _coords[0] + _coords[1] * _coords[1] +
                 _coords[2] * _coords[2];
  return std::sqrt(value);
}

void Point::zero()
{
  _coords[0] = 0.;
  _coords[1] = 0.;
  _coords[2] = 0.;
}
// projection onto plane pasing through p1,p2 and p3
Point Point::projection(const Point &p1,
                        const Point &p2,
                        const Point &p3) const
{
  // normal vector of the plane
  Point n = cross(p3 - p1, p2 - p1);
  // plane equation:Ax+By+Cz+D=0;
  double A = n.x(), B = n.y(), C = n.z();
  double D = -(A * p1.x() + B * p1.y() + C * p1.z());
  double t = -(D + A * x() + B * y() + C * z()) / (A * A + B * B + C * C);
  return Point(A * t + x(), B * t + y(), C * t + z());
}

Point Point::projection_line(const Point &p1, const Point &p2) const
{
  Point t = (p2 - p1).getUnit();
  Point prj = p1 + (((*this) - p1) * t) * t;
  return prj;
}
double Point::distance(const Point &p1) const
{
  return sqrt((p1._coords[0] - _coords[0]) * (p1._coords[0] - _coords[0]) +
              (p1._coords[1] - _coords[1]) * (p1._coords[1] - _coords[1]) +
              (p1._coords[2] - _coords[2]) * (p1._coords[2] - _coords[2]));
}

std::ostream &operator<<(std::ostream &os, const Point &point)
{
  os << point(0) << "\t" << point(1) << "\t" << point(2);
  return os;
}

/*********************friends**************************************/
Point operator*(double a, const Point &p)
{
  return Point(a * p._coords[0], a * p._coords[1], a * p._coords[2]);
}
