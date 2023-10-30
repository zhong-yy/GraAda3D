#include "GravFormula.h"
void GravFormula::field_caused_by_single_prism(const Point &p,
                                               const RectPrism &rect,
                                               double rho,
                                               vector<double> &field,
                                               std::bitset<10> flag)
{
  std::vector<unsigned int> fields_needed;
  const unsigned count = flag.count();
  for (int i = 0; i < flag.size(); i++)
  {
    if (flag[i])
    {
      fields_needed.push_back(i);
    }
  }
  double (GravFormula::*f[10])(const RectPrism &rect, const Point &p,
                               double rho) = {
      &GravFormula::V, &GravFormula::gz, &GravFormula::gx,
      &GravFormula::gy, &GravFormula::Tzz, &GravFormula::Txz,
      &GravFormula::Tyz, &GravFormula::Txx, &GravFormula::Txy,
      &GravFormula::Tyy};

  field.clear();
  field.resize(10);
  for (int i = 0; i < 10; i++)
  {
    field[i] = 0.0;
  }

  for (unsigned i = 0; i < count; i++)
  {
    unsigned index = fields_needed[i];
    field[index] += (this->*f[index])(rect, p, rho);
  }
}

double GravFormula::V(const RectPrism &rect, const Point &p, double rho)
{
  double result = 0;
  Point p0(rect._x[0], rect._y[0], rect._z[0]);
  Point p1(rect._x[1], rect._y[0], rect._z[0]);
  Point p2(rect._x[1], rect._y[1], rect._z[0]);
  Point p3(rect._x[0], rect._y[1], rect._z[0]);
  Point p4(rect._x[0], rect._y[0], rect._z[1]);
  Point p5(rect._x[1], rect._y[0], rect._z[1]);
  Point p6(rect._x[1], rect._y[1], rect._z[1]);
  Point p7(rect._x[0], rect._y[1], rect._z[1]);

  vector<Point *> face1{&p0, &p4, &p7, &p3};
  vector<Point *> face2{&p1, &p2, &p6, &p5};

  vector<Point *> face3{&p0, &p1, &p5, &p4};
  vector<Point *> face4{&p2, &p3, &p7, &p6};

  vector<Point *> face5{&p3, &p2, &p1, &p0};
  vector<Point *> face6{&p4, &p5, &p6, &p7};

  double h1 = -(p0.x() - p.x());
  double h2 = p1.x() - p.x();
  double h3 = -(p0.y() - p.y());
  double h4 = p2.y() - p.y();
  double h5 = -(p0.z() - p.z());
  double h6 = p6.z() - p.z();

  result = 0.5 * GS::G0 * rho *
           (h1 * Integral::Surface_R_1(p, face1, 4) +
            h2 * Integral::Surface_R_1(p, face2, 4) +
            h3 * Integral::Surface_R_1(p, face3, 4) +
            h4 * Integral::Surface_R_1(p, face4, 4) +
            h5 * Integral::Surface_R_1(p, face5, 4) +
            h6 * Integral::Surface_R_1(p, face6, 4));
  return result;
}
double GravFormula::gx(const RectPrism &rect, const Point &p, double rho)
{
  double gx = 0.0;
  Point p0(rect._x[0], rect._y[0], rect._z[0]);
  Point p1(rect._x[1], rect._y[0], rect._z[0]);
  Point p2(rect._x[1], rect._y[1], rect._z[0]);
  Point p3(rect._x[0], rect._y[1], rect._z[0]);
  Point p4(rect._x[0], rect._y[0], rect._z[1]);
  Point p5(rect._x[1], rect._y[0], rect._z[1]);
  Point p6(rect._x[1], rect._y[1], rect._z[1]);
  Point p7(rect._x[0], rect._y[1], rect._z[1]);

  vector<Point *> face1{&p0, &p4, &p7, &p3};
  vector<Point *> face2{&p1, &p2, &p6, &p5};

  double integral1 = Integral::Surface_R_1(p, face1, 4);
  double integral2 = Integral::Surface_R_1(p, face2, 4);

  gx = -GS::SI2mGal * GS::G0 * rho * (-integral1 + integral2);
  return gx;
}

double GravFormula::gy(const RectPrism &rect, const Point &p, double rho)
{
  double gy = 0.0;

  Point p0(rect._x[0], rect._y[0], rect._z[0]);
  Point p1(rect._x[1], rect._y[0], rect._z[0]);
  Point p2(rect._x[1], rect._y[1], rect._z[0]);
  Point p3(rect._x[0], rect._y[1], rect._z[0]);
  Point p4(rect._x[0], rect._y[0], rect._z[1]);
  Point p5(rect._x[1], rect._y[0], rect._z[1]);
  Point p6(rect._x[1], rect._y[1], rect._z[1]);
  Point p7(rect._x[0], rect._y[1], rect._z[1]);

  vector<Point *> face1{&p0, &p1, &p5, &p4};
  vector<Point *> face2{&p2, &p3, &p7, &p6};

  double integral1 = Integral::Surface_R_1(p, face1, 4);
  double integral2 = Integral::Surface_R_1(p, face2, 4);

  gy = -GS::SI2mGal * GS::G0 * rho * (-integral1 + integral2);
  return gy;
}

double GravFormula::gz(const RectPrism &rect, const Point &p, double rho)
{
  double gz = 0.0;

  Point p0(rect._x[0], rect._y[0], rect._z[0]);
  Point p1(rect._x[1], rect._y[0], rect._z[0]);
  Point p2(rect._x[1], rect._y[1], rect._z[0]);
  Point p3(rect._x[0], rect._y[1], rect._z[0]);
  Point p4(rect._x[0], rect._y[0], rect._z[1]);
  Point p5(rect._x[1], rect._y[0], rect._z[1]);
  Point p6(rect._x[1], rect._y[1], rect._z[1]);
  Point p7(rect._x[0], rect._y[1], rect._z[1]);

  vector<Point *> face1{&p3, &p2, &p1, &p0};
  vector<Point *> face2{&p4, &p5, &p6, &p7};

  double integral1 = Integral::Surface_R_1(p, face1, 4);
  double integral2 = Integral::Surface_R_1(p, face2, 4);

  gz = -GS::SI2mGal * GS::G0 * rho * (-integral1 + integral2);
  return gz;
}

double GravFormula::Txx(const RectPrism &rect, const Point &p, double rho)
{
  double txx;

  Point p0(rect._x[0], rect._y[0], rect._z[0]);
  Point p1(rect._x[1], rect._y[0], rect._z[0]);
  Point p2(rect._x[1], rect._y[1], rect._z[0]);
  Point p3(rect._x[0], rect._y[1], rect._z[0]);
  Point p4(rect._x[0], rect._y[0], rect._z[1]);
  Point p5(rect._x[1], rect._y[0], rect._z[1]);
  Point p6(rect._x[1], rect._y[1], rect._z[1]);
  Point p7(rect._x[0], rect._y[1], rect._z[1]);

  vector<Point *> face1{&p0, &p4, &p7, &p3};
  vector<Point *> face2{&p1, &p2, &p6, &p5};
  Point temp1 = Integral::Surface_rR_3(p, face1, 4);
  Point temp2 = Integral::Surface_rR_3(p, face2, 4);
  txx = -GS::SI2Eotvos * GS::G0 * rho * (-temp1.x() + temp2.x());

  return txx;
}
double GravFormula::Txy(const RectPrism &rect, const Point &p, double rho)
{
  double txy;

  Point p0(rect._x[0], rect._y[0], rect._z[0]);
  Point p1(rect._x[1], rect._y[0], rect._z[0]);
  Point p2(rect._x[1], rect._y[1], rect._z[0]);
  Point p3(rect._x[0], rect._y[1], rect._z[0]);
  Point p4(rect._x[0], rect._y[0], rect._z[1]);
  Point p5(rect._x[1], rect._y[0], rect._z[1]);
  Point p6(rect._x[1], rect._y[1], rect._z[1]);
  Point p7(rect._x[0], rect._y[1], rect._z[1]);

  vector<Point *> face1{&p0, &p4, &p7, &p3};
  vector<Point *> face2{&p1, &p2, &p6, &p5};
  Point temp1 = Integral::Surface_rR_3(p, face1, 4);
  Point temp2 = Integral::Surface_rR_3(p, face2, 4);
  txy = -GS::SI2Eotvos * GS::G0 * rho * (-temp1.y() + temp2.y());

  return txy;
}
double GravFormula::Txz(const RectPrism &rect, const Point &p, double rho)
{
  double txz;

  Point p0(rect._x[0], rect._y[0], rect._z[0]);
  Point p1(rect._x[1], rect._y[0], rect._z[0]);
  Point p2(rect._x[1], rect._y[1], rect._z[0]);
  Point p3(rect._x[0], rect._y[1], rect._z[0]);
  Point p4(rect._x[0], rect._y[0], rect._z[1]);
  Point p5(rect._x[1], rect._y[0], rect._z[1]);
  Point p6(rect._x[1], rect._y[1], rect._z[1]);
  Point p7(rect._x[0], rect._y[1], rect._z[1]);

  vector<Point *> face1{&p0, &p4, &p7, &p3};
  vector<Point *> face2{&p1, &p2, &p6, &p5};
  Point temp1 = Integral::Surface_rR_3(p, face1, 4);
  Point temp2 = Integral::Surface_rR_3(p, face2, 4);
  txz = -GS::SI2Eotvos * GS::G0 * rho * (-temp1.z() + temp2.z());

  return txz;
}
double GravFormula::Tyy(const RectPrism &rect, const Point &p, double rho)
{
  double tyy;

  Point p0(rect._x[0], rect._y[0], rect._z[0]);
  Point p1(rect._x[1], rect._y[0], rect._z[0]);
  Point p2(rect._x[1], rect._y[1], rect._z[0]);
  Point p3(rect._x[0], rect._y[1], rect._z[0]);
  Point p4(rect._x[0], rect._y[0], rect._z[1]);
  Point p5(rect._x[1], rect._y[0], rect._z[1]);
  Point p6(rect._x[1], rect._y[1], rect._z[1]);
  Point p7(rect._x[0], rect._y[1], rect._z[1]);

  vector<Point *> face1{&p0, &p1, &p5, &p4};
  vector<Point *> face2{&p2, &p3, &p7, &p6};

  Point temp1 = Integral::Surface_rR_3(p, face1, 4);
  Point temp2 = Integral::Surface_rR_3(p, face2, 4);
  tyy = -GS::SI2Eotvos * GS::G0 * rho * (-temp1.y() + temp2.y());

  return tyy;
}
double GravFormula::Tyz(const RectPrism &rect, const Point &p, double rho)
{
  double tyz;

  Point p0(rect._x[0], rect._y[0], rect._z[0]);
  Point p1(rect._x[1], rect._y[0], rect._z[0]);
  Point p2(rect._x[1], rect._y[1], rect._z[0]);
  Point p3(rect._x[0], rect._y[1], rect._z[0]);
  Point p4(rect._x[0], rect._y[0], rect._z[1]);
  Point p5(rect._x[1], rect._y[0], rect._z[1]);
  Point p6(rect._x[1], rect._y[1], rect._z[1]);
  Point p7(rect._x[0], rect._y[1], rect._z[1]);

  vector<Point *> face1{&p0, &p1, &p5, &p4};
  vector<Point *> face2{&p2, &p3, &p7, &p6};

  Point temp1 = Integral::Surface_rR_3(p, face1, 4);
  Point temp2 = Integral::Surface_rR_3(p, face2, 4);
  tyz = -GS::SI2Eotvos * GS::G0 * rho * (-temp1.z() + temp2.z());

  return tyz;
}

double GravFormula::Tzz(const RectPrism &rect, const Point &p, double rho)
{
  double tzz;

  Point p0(rect._x[0], rect._y[0], rect._z[0]);
  Point p1(rect._x[1], rect._y[0], rect._z[0]);
  Point p2(rect._x[1], rect._y[1], rect._z[0]);
  Point p3(rect._x[0], rect._y[1], rect._z[0]);
  Point p4(rect._x[0], rect._y[0], rect._z[1]);
  Point p5(rect._x[1], rect._y[0], rect._z[1]);
  Point p6(rect._x[1], rect._y[1], rect._z[1]);
  Point p7(rect._x[0], rect._y[1], rect._z[1]);

  vector<Point *> face1{&p3, &p2, &p1, &p0};
  vector<Point *> face2{&p4, &p5, &p6, &p7};

  Point temp1 = Integral::Surface_rR_3(p, face1, 4);
  Point temp2 = Integral::Surface_rR_3(p, face2, 4);
  tzz = -GS::SI2Eotvos * GS::G0 * rho * (-temp1.z() + temp2.z());

  return tzz;
}
