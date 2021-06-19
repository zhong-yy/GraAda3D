#include <chrono>
#include <fstream>
#include <random>
// #include<function>

#include "GaussNewtonInversion.h"
#include "timer.h"
int main() {
  double x_lim[2] = {0, 2000};
  double y_lim[2] = {0, 2000};
  double z_lim[2] = {0, 1000};

  Mesh mesh;
  // mesh.generate_regular_mesh(x_lim, 40, y_lim, 40, z_lim, 20);
  int n_pad_x=9;
  int n_pad_y=9;
  int n_pad_z=9;
  double pad_stretch_h=1.25;
  double pad_stretch_v=1.25;
  mesh.generate_regular_mesh_with_padding(x_lim,40,y_lim,40,z_lim,20,n_pad_x,pad_stretch_h,n_pad_y,pad_stretch_h,n_pad_z,pad_stretch_v);

  double block_z[2] = {200, 500};
  double block_x[2] = {600, 900};
  double block_y[2] = {550, 850};
  mesh.set_parameter_in_a_region(block_x, block_y, block_z, 500);

  double block_z2[2] = {200, 500};
  double block_x2[2] = {1200, 1500};
  double block_y2[2] = {550, 850};
  mesh.set_parameter_in_a_region(block_x2, block_y2, block_z2, -500);

  double block_z3[2] = {200, 500};
  double block_x3[2] = {900, 1500};
  double block_y3[2] = {1300, 1500};
  mesh.set_parameter_in_a_region(block_x3, block_y3, block_z3, 500);

  // mesh.set_block_density(12, 27, 12, 27, 6, 15, 500);//
  // mesh.set_block_density(6, 13, 6, 13, 6, 15, 500); //

  mesh.out_model_vtk("test_model_with_padding.vtk");
#ifdef USE_NETCDF
  mesh.out_model_netcdf("test_model_with_padding.nc");
#endif

  return 0;
}
