#ifndef _CELL
#define _CELL
#include "Face.h"
#include "RectPrism.h"

#include <Eigen/Dense>  // Linear Algebra Lib.
#include <Eigen/Sparse> // Sparse Lib
#include <Eigen/StdVector>
#include <algorithm>
#include <cassert>
#include <string>
#include <iostream>
using namespace std;
using namespace Eigen;

class Face;
class Cell : public RectPrism
{
public:
  // ~Cell();
  Cell(double x0,
       double y0,
       double z0,
       double x1,
       double y1,
       double z1,
       int level_ = 0,
       int num_para = 1,
       bool isleaf_ = true);
  bool get_ordering_forward() { return this->ordering_forward; }
  void set_ordering_forward(bool forward)
  {
    this->ordering_forward = forward;
  }
  void set_id(int _id) { this->id = _id; }
  int get_id() { return id; }
  void set_parameter(double new_value, int i = 0)
  {
    this->parameters[i] = new_value;
    if (!isleaf)
    {
      for (int j = 0; j < 8; j++)
      {
        child_cells[j]->set_parameter(new_value, i);
      }
    }
  }
  double get_parameter(int i = 0) { return this->parameters[i]; }
  void set_external_faces(Face *f1, Face *f2, unsigned int normal_dirction);

  void set_internal_faces(Face *f[4], unsigned int normal_dirction);

  int get_level() { return level; }

  //   protected:
  bool ordering_forward;
  int level; // first level is 0
  // Rect_prism *cell;
  Cell *child_cells[8];

  Face *external_faces_z[2]; // external_faces_r
  Face *external_faces_x[2]; // external_faces_theta
  Face *external_faces_y[2]; // external_faces_phi

  Face *internal_faces_z[4]; // r
  Face *internal_faces_x[4]; // theta
  Face *internal_faces_y[4]; // phi

  int id;
  bool isleaf;
  vector<double> parameters;
  // double density;
};
#endif
