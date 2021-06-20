#ifndef _MESH
#define _MESH

#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>

#include <set>
#include <unordered_map>
#include <vector>

#include "Cell.h"
#include "Point.h"

#ifdef USE_NETCDF
#include <netcdf>  //The library netcdf-cxx is required
using namespace netCDF;
using namespace netCDF::exceptions;
#define NC_ERR 2
#endif

class Mesh {
 public:
  Mesh();
  Mesh(const Mesh&);
  Mesh& operator=(const Mesh&);
  ~Mesh();

  void clear_all();

  void set_n_parameter(int n) {
    this->n_parameters = n;
    for (int i = 0; i < leaf_cells.size(); i++) {
      leaf_cells[i]->parameters.resize(n);
      for (int j = 0; j < n; j++) {
        leaf_cells[i]->set_parameter(0, j);
      }
    }
  }

  // Generate a regular mesh
  void generate_regular_mesh(
      double
          x_model[2],  // the dimension of the area of interest in x direction
      int n_x,         // the number of cells in x direction
      double
          y_model[2],  // the dimension of the area of interest in y direction
      int n_y,         // the number of cells in y direction
      double
          z_model[2],  // the dimension of the area of interest in z direction
      int n_z,         // the number of cells in z direction
      int num_para =
          1);  // the number of types of parameters will be stored in this model

  void generate_regular_mesh(double x_model[2],
                             int n_x,
                             double y_model[2],
                             int n_y,
                             VectorXd& z_points0,
                             int num_para = 1);

  void generate_regular_mesh(
      VectorXd& x_points0,  // Coordinates of grid points in the x direction
      VectorXd& y_points0,
      VectorXd& z_points0,
      int num_para = 1);

  void generate_regular_mesh_with_padding(
      double x_model[2],  // dimension of the area of interest in x direction
      int n_x,            // number of cells in x direction
      double y_model[2],  // dimension of the area of interest in y direction
      int n_y,            // number of cells in y direction
      double z_model[2],  // dimension of the area of interest in z direction
      int n_z,            // number of cells in z direction
      int n_pad_x = 7,    // number of padding cells on +x and -x sides
      double pad_stretch_x =
          1.5,  // increasing factor for sizes of padding cells in x direction
      int n_pad_y = 7,  // number of padding cells on +y and -y sides
      double pad_stretch_y =
          1.5,  // increasing factor for sizes of padding cells in y direction
      int n_pad_z = 3,  // number of padding cells on +z and -z sides
      double pad_stretch_z =
          1.2,  // increasing factor for sizes of padding cells  in z direction
      int num_para = 1);

  void set_parameter_in_a_region(double x[2],
                                 double y[2],
                                 double z[2],
                                 double para_value,
                                 int i_th = 0);
  void set_parameter_in_a_region(double x0,
                                 double x1,
                                 double y0,
                                 double y1,
                                 double z0,
                                 double z1,
                                 double para_value,
                                 int i_th = 0);

  void set_block_parameter(unsigned int i_min,
                           unsigned int i_max,
                           unsigned int j_min,
                           unsigned int j_max,
                           unsigned int k_min,
                           unsigned int k_max,
                           double para_value,
                           int i_th = 0);

  void get_minimum_size(double& dx, double& dy, double& dz, int& lev) {
    assert(cells[cells.size() - 1].size() > 0);
    int num = cells[cells.size() - 1].size();
    Cell* c = cells[cells.size() - 1][0];
    lev = c->get_level();
    c->get_size(dx, dy, dz);
    for (int i = 1; i < num; i++) {
      c = cells[cells.size() - 1][i];
      double dx0, dy0, dz0;
      c->get_size(dx0, dy0, dz0);
      int lev0 = c->get_level();
      if ((dx0 * dy0 * dz0) < (dx * dy * dz)) {
        dx = dx0;
        dy = dy0;
        dz = dz0;
        lev = lev0;
      }
    }
  }

  void set_ith_cell_parameter(double para_value, int id, int i) {
    leaf_cells[id]->set_parameter(para_value, i);
  };

  // return a key-value map, whose key is the id of the cell that is to be
  // subdivided, and the value is the pointer to the cell
  map<unsigned int, Cell*> refinement(int i);
  map<unsigned int, Cell*> refinement(Cell* c);

  void out_model_vtk(string filename,
                     int n = 1,
                     vector<string> parameter_name = vector<string>(1,
                                                                    "density"));
  void out_model_vtk_points(string filename,
                     int n = 1,
                     vector<string> parameter_name = vector<string>(1,
                                                                    "density"));
#ifdef USE_NETCDF
  int out_model_netcdf(string filename,
                       int ith_para = 0,
                       string VAL_NAME = "density",
                       string VAL_UNITS = "kg/m3");
#endif

  void fill_data(int offset_i,
                 int offset_j,
                 int offset_k,
                 double*** data,
                 Cell* c,
                 int max_level,
                 int ith_para);

  int n_elems() const { return leaf_cells.size(); }

  void rearrange_id();

  bool great_equal(long double left, long double right);

  void get_model_parameter_from_mesh(VectorXd& m, int ith = 0);

  void sort(int level);

  RectPrism& get_elem(unsigned int i);

  Cell* get_element_level_0(unsigned int i_lat,
                            unsigned int j_lon,
                            unsigned int k_r);
  int get_nz_level_0() const { return nz; }
  int get_nx_level_0() const { return nx; }
  int get_ny_level_0() const { return ny; }
  //   protected:
  // vector<Cell *> root_cells;
  // vector<Face *> root_faces;

  vector<Cell*> leaf_cells;
  vector<Face*> leaf_faces;
  int num_leaf_cells;
  int num_leaf_faces;

  vector<vector<Cell*>> cells;  // cells of different levels
  vector<vector<Face*>> faces;  // faces of different levels

  vector<int> num_cell;
  vector<int> num_face;

  int n_parameters;

  int nz, nx, ny;  // number of cells along x, y, z directions for level 0
  // vector<double> r_intervals;
  VectorXd z_points;
  VectorXd x_points;
  VectorXd y_points;
  double x_lim[2], y_lim[2], z_lim[2];
};
#endif