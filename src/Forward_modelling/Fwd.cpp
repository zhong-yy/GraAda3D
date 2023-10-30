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

Fwd::Fwd(const Mesh &mesh_, const Observation &ob_,
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
void Fwd::GT_vec_mul(const VectorXd &vec, VectorXd &product) const
{
  int n0 = vec.size();
  // assert(n0 == Nd);
  product.resize(Nm);
  if (use_wavelet)
  {
    vector<double> product_wavelet(n_dy_for_wavelet, 0);
    // #pragma omp parallel for
    for (int i = 0; i < Nd; i++)
    {
      int nnz_of_current_row = 0;
      nnz_of_current_row = comp_col_ids[i].size();
#pragma omp parallel for
      for (int j = 0; j < nnz_of_current_row; j++)
      {
        int id = comp_col_ids[i][j];
        // #pragma omp atomic
        product_wavelet[id] += comp_G_coeffs[i][j] * vec(i);
      }
    }
    inverse_wavelet_transform_vec(product_wavelet);
#pragma omp parallel for
    for (int i = 0; i < Nm; i++)
    {
      int original_id_of_i_re_oredered_cell = mesh.get_reordered_id(i);
      product(original_id_of_i_re_oredered_cell) = product_wavelet[i];
    }
  }
  else
  {
    product = G.transpose() * vec;
  }
}
void Fwd::G_vec_mul(const VectorXd &vec, VectorXd &product) const
{
  int n0 = vec.size();
  // assert(n0 == Nm);
  product.resize(Nd);
  if (use_wavelet)
  {
    vector<double> wavelet_coeffs(n0);
#pragma omp parallel for
    for (int i = 0; i < n0; i++)
    {
      wavelet_coeffs[i] = vec[mesh.get_reordered_id(i)];
    }
    int n_coeff = wavelet_transform_vec(wavelet_coeffs);
    // assert(n_coeff == n_dy_for_wavelet);

#pragma omp parallel for
    for (int i = 0; i < Nd; ++i)
    {
      int n_non_zero = 0;
      n_non_zero = comp_col_ids[i].size();
      product(i) = 0.0;
      for (int j = 0; j < n_non_zero; ++j)
      {
        int col_id = comp_col_ids[i][j];
        product(i) += comp_G_coeffs[i][j] * wavelet_coeffs[col_id];
      }
    }
  }
  else
  {
    product = G * vec;
  }
}
void Fwd::compute_G_wavelet()
{
  int n_fields = field_flag.count();
  assert(Nm > 0);
  assert(N_obs > 0);

  comp_col_ids.clear();
  comp_col_ids.resize(n_fields * N_obs);
  comp_G_coeffs.clear();
  comp_G_coeffs.resize(n_fields * N_obs);

  G.resize(0, 0); // In the wavelet mode, G matrix is not longer used. Use G_vec_mul() and GT_vec_mul() instead.

  vector<int> basic_field_index;
  for (int i = 0; i < 10; i++)
  {
    if (field_flag[i])
    {
      basic_field_index.push_back(i);
    }
  }
  cout << "Relative threshold for wavelet compression: "
       << this->compression_threshold << endl;
  vector<vector<double>> sensitivity_data(basic_field_index.size());
  for (int ii = 0; ii < sensitivity_data.size(); ii++)
  {
    sensitivity_data[ii].resize(Nm);
    for (int jj = 0; jj < Nm; jj++)
    {
      sensitivity_data[ii][jj] = 0.0;
    }
  }
#pragma omp parallel for firstprivate(sensitivity_data)
  for (int i = 0; i < N_obs; i++)
  {
    for (int j = 0; j < Nm; j++)
    { // here j denotes the reordered index
      GravFormula gra;
      vector<double> field;
      int ro_id = mesh.get_reordered_id(j);
      gra.field_caused_by_single_prism(ob(i), mesh.get_elem(ro_id), 1.0, field,
                                       field_flag);
      assert(basic_field_index.size() == n_fields);

      for (int k = 0; k < n_fields; k++)
      {
        sensitivity_data[k][j] = field[basic_field_index[k]];
        // G(i + k * N_obs, j) = field[basic_field_index[k]];
      }
    }

    for (int k = 0; k < n_fields; ++k)
    {
      this->n_dy_for_wavelet = compress_vec(
          sensitivity_data[k], comp_col_ids[i + k * N_obs],
          comp_G_coeffs[i + k * N_obs], this->compression_threshold);
      assert(n_dy_for_wavelet == sensitivity_data[k].size());
      // comp_col_ids[i + k * N_obs].clear();
      // comp_G_coeffs[i + k * N_obs].clear();
      // for (int j = 0; j < n_dy_for_wavelet; ++j) {
      //    if (std::fabs(sensitivity_data[k][j]) > 1e-15) {
      //        n_non_zero++;
      //        comp_col_ids[i + k * N_obs].push_back(j);
      //        comp_G_coeffs[i + k * N_obs].push_back(
      //            sensitivity_data[k][j]);
      //    }
      //}
    }
  }
  int n_non_zero = 0;
#pragma omp parallel for reduction(+ : n_non_zero)
  for (int i = 0; i < comp_col_ids.size(); ++i)
  {
    n_non_zero += comp_col_ids[i].size();
  }
  cout << "n_dy=" << n_dy_for_wavelet << endl;
  double compression_ratio = (n_fields * N_obs * Nm) / (1.0 * n_non_zero);
  cout << "The number of elements in the sensitivity matrix: "
       << (n_fields * N_obs * Nm)
       << ", the number of non-zero elements: " << n_non_zero << endl;
  cout << "Compression ratio: " << compression_ratio << endl;
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
  vector<string> strs = {"V", "g_z", "g_x", "g_y", "T_zz",
                         "T_xz", "T_yz", "T_xx", "T_xy", "T_yy"};
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