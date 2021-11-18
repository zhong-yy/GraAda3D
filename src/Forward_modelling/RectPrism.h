#ifndef RECTPRISM_H
#define RECTPRISM_H
#include <cassert>
#include <string>
// #include "Point.h"
using namespace std;
#include "gs.h"
/**
 * @brief Rectangular prism
 *
 */
class RectPrism {
 public:
  RectPrism(double x[2], double y[2], double z[2]);

  RectPrism(double x0, double y0, double z0, double x1, double y1, double z1);

  RectPrism();
  ~RectPrism(){};

  double get_volumn() const {
    double dx, dy, dz;
    this->get_size(dx, dy, dz);
    return (dx * dy * dz);
  }
  void get_size(double& dx, double& dy, double& dz) const;

  void get_center(double& xcenter, double& ycenter, double& zcenter) const;

  void get_xlim(double& x0, double& x1)const;
  void get_ylim(double& y0, double& y1)const;
  void get_zlim(double& z0, double& z1)const;

 public:
  // int id;
  // double dx, dy, dz;
  // Point minPoint;

  double _x[2];
  double _y[2];
  double _z[2];
};

#endif
