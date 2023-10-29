#include "Fwd.h"

Fwd::Fwd() : field_flag(Compute_g_z)
{
  // this->mesh = NULL;
  // this->ob = NULL;
  Nm = 0;
  N_obs = 0;
  int n_fields = field_flag.count();
  Nd = N_obs * n_fields;
  use_wavelet = false;
  compression_threshold = 0;
}

Fwd::~Fwd() {}

Fwd::Fwd(const Mesh &mesh_,
         const Observation &ob_,
         unsigned long long field_flag_)
    : field_flag(field_flag_)
{
  mesh = mesh_;
  ob = ob_;
  Nm = mesh.n_elems();
  N_obs = ob.get_n_obs();
  int n_fields = field_flag.count();
  Nd = N_obs * n_fields;
  use_wavelet = false;
  compression_threshold = 0;
  display_info_fields();
}

Fwd::Fwd(const Mesh &mesh_, const Observation &ob_, bitset<10> field_flag_)
    : field_flag(field_flag_)
{
  mesh = mesh_;
  ob = ob_;

  Nm = mesh.n_elems();
  N_obs = ob.get_n_obs();
  int n_fields = field_flag.count();
  Nd = N_obs * n_fields;
  use_wavelet = false;
  compression_threshold = 0;
  display_info_fields();
}

void Fwd::compute_G()
{
  int n_fields = field_flag.count();
  assert(Nm > 0);
  assert(N_obs > 0);

  G.resize(n_fields * N_obs, Nm);

  vector<int> basic_field_index;
  for (int i = 0; i < 10; i++)
  {
    // bool status = flag[i];
    if (field_flag[i])
    {
      basic_field_index.push_back(i);
    }
  }

  for (int i = 0; i < N_obs; i++)
  {
#pragma omp parallel for
    for (int j = 0; j < Nm; j++)
    {
      GravFormula gra;
      vector<double> field;
      gra.field_caused_by_single_prism(ob(i), mesh.get_elem(j), 1.0, field,
                                       field_flag);
      for (int k = 0; k < n_fields; k++)
      {
        G(i + k * N_obs, j) = field[basic_field_index[k]];
      }
    }
  }
}

void Fwd::set_mesh(const Mesh &mesh0)
{
  this->mesh = mesh0;
  Nm = mesh.n_elems();
}

void Fwd::set_observation(const Observation &ob0)
{
  this->ob = ob0;
  N_obs = ob.get_n_obs();
}

void Fwd::display_info_fields() const
{
  int n = field_flag.count();
  vector<string> strs = {"V", "g_z", "g_x", "g_y",
                         "T_zz", "T_xz", "T_yz", "T_xx",
                         "T_xy", "T_yy"};
  vector<int> to_be_computed;
  for (int i = 0; i < strs.size(); i++)
  {
    if (field_flag[i])
    {
      to_be_computed.push_back(i);
    }
  }

  cout << n
       << (n > 1 ? " gravity field components are"
                 : " gravity field component is")
       << " used. " << (n > 1 ? "They are: \t" : "It is: \t");
  for (int i = 0; i < n; i++)
  {
    cout << strs[to_be_computed[i]] << ((n > 1 && i == n - 2) ? " and " : ", ");
  }
  cout << endl;
}

int Fwd::inverse_wavelet_transform_vec(vector<double> &data) const
{
  int n = data.size();
  int n_dy;

  double log_value = std::log2(1.0 * n);
  if (std::fabs(log_value - int(log_value)) < 1e-15)
  {
    n_dy = n;
  }
  else
  {
    n_dy = std::pow(2, ceil(log2(1.0 * n)));
    vector<double> data_copy = data;
    data.resize(n_dy);
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
      data[i] = data_copy[i];
    }
    // padding
#pragma omp parallel for
    for (int i = n; i < n_dy; i++)
    {
      data[i] = 0;
    }
  }

  gsl_wavelet *w;
  gsl_wavelet_workspace *work;
  w = gsl_wavelet_alloc(gsl_wavelet_daubechies, 4);
  work = gsl_wavelet_workspace_alloc(n_dy);

  double *p_data = data.data();
  gsl_wavelet_transform_inverse(w, p_data, 1, n_dy, work);
  gsl_wavelet_free(w);
  gsl_wavelet_workspace_free(work);
  return n_dy;
}

int Fwd::compress_vec(vector<double> &data, vector<int> &ids,
                      vector<double> &coeffs, double relative_threshold)
{
  int n = data.size();
  int n_dy;

  double log_value = std::log2(1.0 * n);
  if (std::fabs(log_value - int(log_value)) < 1e-15)
  {
    n_dy = n;
  }
  else
  {
    n_dy = std::pow(2, ceil(log2(1.0 * n)));
    vector<double> data_copy = data;
    data.resize(n_dy);
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
      data[i] = data_copy[i];
    }
    // padding
#pragma omp parallel for
    for (int i = n; i < n_dy; i++)
    {
      data[i] = 0;
    }
  }

  gsl_wavelet *w;
  gsl_wavelet_workspace *work;
  w = gsl_wavelet_alloc(gsl_wavelet_daubechies, 4);
  work = gsl_wavelet_workspace_alloc(n_dy);

  double *p_data = data.data();
  gsl_wavelet_transform_forward(w, p_data, 1, n_dy, work);

  double max_abscoeff;
  // auto comp = [](double a, double b) { return std::fabs(a) < std::fabs(b);
  // }; max_abscoeff = std::fabs(*(max_element(data.begin(), data.end(),
  // comp)));
  max_abscoeff = std::fabs(data[0]);
  for (int i = 1; i < data.size(); i++)
  {
    if (max_abscoeff < std::fabs(data[i]))
    {
      max_abscoeff = std::fabs(data[i]);
    }
  }

  ids.clear();
  coeffs.clear();
  // int n_zeros = 0;
  double abs_thre = (relative_threshold * max_abscoeff);
  for (int i = 0; i < n_dy; i++)
  {
    if ((std::fabs(data[i]) > abs_thre) ||
        (std::fabs(std::fabs(data[i]) - abs_thre) < 1e-15))
    {
      ids.push_back(i);
      coeffs.push_back(data[i]);
      //        data[i] = 0;
      //        n_zeros++;
    }
  }
  gsl_wavelet_free(w);
  gsl_wavelet_workspace_free(work);
  // printf("Max abscoeff=%f\n", max_abscoeff);
  // printf("Number of padding elements: %d\n", n_dy - n);
  // printf("number of zero coefficients=%d\n", n_zeros);
  // printf("number of non-zero coefficients=%d\n", n_dy - n_zeros);
  // printf("Compression ratio=%f\n", n * 1.0 / (1.0 * n_dy - 1.0 * n_zeros));
  return n_dy;
}

int Fwd::wavelet_transform_vec(vector<double> &data) const
{
  int n = data.size();
  int n_dy;

  double log_value = std::log2(1.0 * n);
  if (std::fabs(log_value - int(log_value)) < 1e-15)
  {
    n_dy = n;
  }
  else
  {
    n_dy = std::pow(2, ceil(log2(1.0 * n)));
    vector<double> data_copy = data;
    data.resize(n_dy);
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
      data[i] = data_copy[i];
    }
    // padding
#pragma omp parallel for
    for (int i = n; i < n_dy; i++)
    {
      data[i] = 0;
    }
  }

  gsl_wavelet *w;
  gsl_wavelet_workspace *work;
  w = gsl_wavelet_alloc(gsl_wavelet_daubechies, 4);
  work = gsl_wavelet_workspace_alloc(n_dy);

  double *p_data = data.data();
  gsl_wavelet_transform_forward(w, p_data, 1, n_dy, work);

  gsl_wavelet_free(w);
  gsl_wavelet_workspace_free(work);
  // printf("Max abscoeff=%f\n", max_abscoeff);
  // printf("Number of padding elements: %d\n", n_dy - n);
  // printf("number of zero coefficients=%d\n", n_zeros);
  // printf("number of non-zero coefficients=%d\n", n_dy - n_zeros);
  // printf("Compression ratio=%f\n", n * 1.0 / (1.0 * n_dy - 1.0 * n_zeros));
  return n_dy;
}

void Fwd::set_use_wavelet(bool use_wavelet0)
{
  this->use_wavelet = use_wavelet0;
}