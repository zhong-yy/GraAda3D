#include "Fwd.h"

Fwd::Fwd() : field_flag(Compute_g_z) {
  // this->mesh = NULL;
  // this->ob = NULL;
  Nm = 0;
  N_obs = 0;
  int n_fields = field_flag.count();
  Nd = N_obs * n_fields;
}

Fwd::~Fwd() {}

Fwd::Fwd(const Mesh& mesh_,
         const Observation& ob_,
         unsigned long long field_flag_)
    : field_flag(field_flag_) {
  mesh = mesh_;
  ob = ob_;
  Nm = mesh.n_elems();
  N_obs = ob.get_n_obs();
  int n_fields = field_flag.count();
  Nd = N_obs * n_fields;
  display_info_fields();
}

Fwd::Fwd(const Mesh& mesh_, const Observation& ob_, bitset<10> field_flag_)
    : field_flag(field_flag_) {
  mesh = mesh_;
  ob = ob_;

  Nm = mesh.n_elems();
  N_obs = ob.get_n_obs();
  int n_fields = field_flag.count();
  Nd = N_obs * n_fields;
  display_info_fields();
}

void Fwd::compute_G() {
  int n_fields = field_flag.count();
  assert(Nm > 0);
  assert(N_obs > 0);

  G.resize(n_fields * N_obs, Nm);

  vector<int> basic_field_index;
  for (int i = 0; i < 10; i++) {
    // bool status = flag[i];
    if (field_flag[i]) {
      basic_field_index.push_back(i);
    }
  }

  for (int i = 0; i < N_obs; i++) {
#pragma omp parallel for
    for (int j = 0; j < Nm; j++) {
      GravFormula gra;
      vector<double> field;
      gra.field_caused_by_single_prism(ob(i), mesh.get_elem(j), 1.0, field,
                                field_flag);
      for (int k = 0; k < n_fields; k++) {
        G(i + k * N_obs, j) = field[basic_field_index[k]];
      }
    }
  }
}

void Fwd::set_mesh(const Mesh& mesh0) {
  this->mesh = mesh0;
  Nm = mesh.n_elems();
}

void Fwd::set_observation(const Observation& ob0) {
  this->ob = ob0;
  N_obs = ob.get_n_obs();
}

void Fwd::display_info_fields() const {
  int n = field_flag.count();
  vector<string> strs = {"V",          "g_z",      "g_x", "g_y",
                         "T_zz",       "T_xz", "T_yz",  "T_xx",
                         "T_xy", "T_yy"};
  vector<int> to_be_computed;
  for (int i = 0; i < strs.size(); i++) {
    if (field_flag[i]) {
      to_be_computed.push_back(i);
    }
  }

  cout << n
       << (n > 1 ? " gravity field components are"
                 : " gravity field component is")
       << " used. " << (n > 1 ? "They are: \t" : "It is: \t");
  for (int i = 0; i < n; i++) {
    cout << strs[to_be_computed[i]] << ((n > 1 && i == n - 2) ? " and " : ", ");
  }
  cout << endl;
}