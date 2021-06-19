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
  mesh.generate_regular_mesh(x_lim, 40, y_lim, 40, z_lim, 20);

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

  mesh.out_model_vtk("test_model.vtk");
#ifdef USE_NETCDF
  mesh.out_model_netcdf("test_model.nc");
#endif

  Mesh mesh2;
  mesh2.generate_regular_mesh(x_lim, 40, y_lim, 40, z_lim, 20);
  mesh2.set_parameter_in_a_region(block_x3, block_y3, block_z3, 200);
  mesh2.out_model_vtk("apriori_model.vtk");
#ifdef USE_NETCDF
  mesh2.out_model_netcdf("apriori_model.nc");
#endif

  ofstream out_cross_gradient_constraint("crg_model");
  for (int k = 0; k < mesh2.get_nz_level_0(); k++) {
    for (int j = 0; j < mesh2.get_ny_level_0(); j++) {
      for (int i = 0; i < mesh2.get_nx_level_0(); i++) {
        Cell* c = mesh2.get_element_level_0(i, j, k);
        double xc, yc, zc;
        c->get_center(xc, yc, zc);
        out_cross_gradient_constraint
            << setw(15) << fixed << setprecision(7) << left << xc << setw(15)
            << setprecision(7) << left << yc << setw(15) << setprecision(3)
            << left << zc << scientific << setprecision(15) << left
            << c->get_parameter() << endl;
      }
    }
  }

  function<double(double)> relation = [](double x) -> double {
    return (2.0 * x + 50.0);
  };  //
  ofstream out_petrophysical_constraint("ref_model");
  for (int k = 0; k < mesh2.get_nz_level_0(); k++) {
    for (int j = 0; j < mesh2.get_ny_level_0(); j++) {
      for (int i = 0; i < mesh2.get_nx_level_0(); i++) {
        Cell* c = mesh2.get_element_level_0(i, j, k);
        double xc, yc, zc;
        c->get_center(xc, yc, zc);
        out_petrophysical_constraint
            << setw(15) << fixed << setprecision(7) << left << xc << setw(15)
            << setprecision(7) << left << yc << setw(15) << setprecision(3)
            << left << zc << scientific << setprecision(15) << left
            << relation(c->get_parameter()) << endl;
      }
    }
  }

  /******************Foward modelling*****************/
  VectorXd rho;
  mesh.get_model_parameter_from_mesh(rho, 0);
  Observation ob;

  int nx = 41;
  int ny = 41;
  ob.generate(nx, x_lim[0], x_lim[1], ny, y_lim[0], y_lim[1], -0.1);
  ofstream os("sites");
  os << ob;

  Fwd forward(mesh, ob, Compute_g_z | Compute_T_zz|Compute_T_zx|Compute_T_zy);

  Timer timer;
  timer.start();
  cout << "Generating synthetic data..." << endl;
  forward.compute_G();
  const Eigen::MatrixXd& G = forward.get_G();
  Eigen::VectorXd d_obs = G * rho;
  cout << "Calculation of synthetic data completed" << endl;
  timer.stop();
  cout << "Time: " << timer.getElapsedTimeInSec() << " s" << endl;

  
  VectorXd noise;
  noise.resize(d_obs.rows());
  double equipment_noise = 0.0;
  for (int i = 0; i < d_obs.rows(); i++) {
    static normal_distribution<double> normal_dist(0, 1);
    static default_random_engine e(time(0));
    noise(i) = (0.02 * fabs(d_obs(i)) + equipment_noise) * normal_dist(e);
    d_obs(i) = d_obs(i) + noise(i);
  }
  // cout << d_obs.rows() << endl;

  GaussNewtonInversion write_data(mesh, ob, Compute_g_z | Compute_T_zz|Compute_T_zx|Compute_T_zy);

  write_data.set_dobs(d_obs);
  write_data.output_obs_data("dobs");

  return 0;
}
