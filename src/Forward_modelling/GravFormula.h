#ifndef GRAVFORMULA_H
#define GRAVFORMULA_H
#include <bitset>
#include "Integral.h"
#include "RectPrism.h"
/**
 * @brief A group of formulae to compute gravitational potential, acceleration
 * and gradient tensor.
 * 
 * Length: m
 * Potential V: m^2/s^2
 * Gravity anomaly: mGal (1e-5 m/s^2)
 * Gravity gradient tensor: E (1e-9 1/s^2)
 */
class GravFormula {
 public:
  void field_caused_by_single_prism(const Point& p,
                                    const RectPrism& rect,
                                    double rho,
                                    vector<double>& field,
                                    std::bitset<10> flag);

  double V(const RectPrism& rect, const Point& p, double rho);


  double gx(const RectPrism& rect, const Point& p, double rho);


  double gy(const RectPrism& rect, const Point& p, double rho);

  double gz(const RectPrism& rect, const Point& p, double rho);


  double Txx(const RectPrism& rect, const Point& p, double rho);


  double Txy(const RectPrism& rect, const Point& p, double rho);


  double Txz(const RectPrism& rect, const Point& p, double rho);


  double Tyy(const RectPrism& rect, const Point& p, double rho);


  double Tyz(const RectPrism& rect, const Point& p, double rho);


  double Tzz(const RectPrism& rect, const Point& p, double rho);
};
#endif