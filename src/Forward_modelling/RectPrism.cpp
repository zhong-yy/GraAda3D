#include "RectPrism.h"
#include <string>
RectPrism::RectPrism(double x0,
                     double y0,
                     double z0,
                     double x1,
                     double y1,
                     double z1) {
  this->_x[0] = x0;
  this->_x[1] = x1;
  this->_y[0] = y0;
  this->_y[1] = y1;
  this->_z[0] = z0;
  this->_z[1] = z1;
}

RectPrism::RectPrism(double x[2], double y[2], double z[2]) {
  for (int i = 0; i < 2; i++) {
    this->_x[i] = x[i];
    this->_y[i] = y[i];
    this->_z[i] = z[i];
  }
}

RectPrism::RectPrism() {
  this->_x[0] = 0;
  this->_x[1] = 0;
  this->_y[0] = 0;
  this->_y[1] = 0;
  this->_z[0] = 0;
  this->_z[1] = 0;
}

void RectPrism::get_size(double& dx, double& dy, double& dz) const {
  dx = abs(_x[1] - _x[0]);
  dy = abs(_y[1] - _y[0]);
  dz = abs(_z[1] - _z[0]);
}

void RectPrism::get_center(double& xc, double& yc, double& zc) const {
  xc = 0.5 * (_x[1] + _x[0]);
  yc = 0.5 * (_y[1] + _y[0]);
  zc = 0.5 * (_z[1] + _z[0]);
}
