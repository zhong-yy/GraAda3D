#ifndef _FWD
#define _FWD
#include <gsl/gsl_sort.h>
#include <gsl/gsl_wavelet.h>

#include <bitset>

#include "GravFormula.h"
#include "Mesh.h"
#include "Observation.h"
#include "gs.h"
using namespace std;

class Fwd
{
public:
  Fwd();
  Fwd(const Mesh &mesh_, const Observation &ob_,
      unsigned long long field_flag_ = Compute_g_z);
  Fwd(const Mesh &mesh_, const Observation &ob_, bitset<10> field_flag);
  virtual ~Fwd();

  VectorXd compute_gobs(const VectorXd &rho)
  {
    VectorXd d = this->G * rho;
    return d;
  }

  void set_field_flag(unsigned long long field_flag1)
  {
    field_flag = field_flag1;
  }

  void compute_G();
  void compute_G_wavelet();

  void set_use_wavelet(bool use_wavelet0);
  // bind mesh of tesseroids
  void set_mesh(const Mesh &mesh0);

  // bind computation points
  void set_observation(const Observation &ob0);

  const Eigen::MatrixXd &get_G() { return G; }

  unsigned int get_n_fields() { return field_flag.count(); }

  Mesh &get_mesh() { return mesh; }
  void set_compression_threshold(double eps)
  {
    assert((std::fabs(eps) < 1e-15 || eps > 0) && (eps < 1));
    this->compression_threshold = eps;
  }

protected:
  Eigen::MatrixXd G;
  Mesh mesh;
  Observation ob;
  int Nm;    // number of parameters
  int N_obs; // number of data points
  int Nd;    // number of observations

  bitset<10> field_flag;

  // wavelet support
  bool use_wavelet;
  int n_dy_for_wavelet;
  // relative threshold for wavelet compression, default 0.005
  double compression_threshold;

  vector<vector<int>> comp_col_ids;
  vector<vector<double>> comp_G_coeffs;

  int inverse_wavelet_transform_vec(vector<double> &data) const;
  int compress_vec(vector<double> &data, vector<int> &ids,
                   vector<double> &coeffs, double relative_threshold);
  int wavelet_transform_vec(vector<double> &data) const;

  void display_info_fields() const;
  /**
   * 0 V
   * 1 g_z
   * 2 g_x
   * 3 g_y
   * 4 T_zz
   * 5 T_xz
   * 6 T_yz
   * 7 T_xx
   * 8 T_xy
   * 9 T_yy
   */
};
#endif
