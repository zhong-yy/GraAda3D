#include "InversionBase.h"
InversionBase::InversionBase()
    : Fwd(),
      a_s(1e-3),
      a_z(1e0),
      a_x(1e0),
      a_y(1e0),
      a_crg(1e0),
      target_misfit(1),
      max_lambda(1e6),
      depth_weighting_factor(2),
      n_lambda(20),
      lambda_decreasing_rate(0.5),
      interpolator_m0(NULL),
      interpolator_m0s(NULL),
      use_cross_gradient_constraint(false),
      use_petrophysical_constraint(false) {
  mesh.set_n_parameter(5);
}

InversionBase::InversionBase(const Mesh& mesh_,
                             const Observation& ob_,
                             unsigned long long field_flag_)
    : Fwd(mesh_, ob_, field_flag_),
      target_misfit(1),
      max_lambda(1e6),
      depth_weighting_factor(2),
      n_lambda(20),
      lambda_decreasing_rate(0.5),
      interpolator_m0(NULL),
      interpolator_m0s(NULL),
      use_cross_gradient_constraint(false),
      use_petrophysical_constraint(false) {
  int n_field = field_flag.count();
  this->init_matrices();

  // initialise the reference model as zero
  m0 = VectorXd::Constant(Nm, 1, 0);
  m_min = VectorXd::Constant(Nm, 1, -1e6);
  m_max = VectorXd::Constant(Nm, 1, 1e6);
  m0_s = VectorXd::Constant(Nm, 1, 0);
  // cout<<"A"<<endl;
  Wd.resize(Nd, Nd);
  Wd.setIdentity();

  // cout << "here" << endl;
  m.resize(Nm);

  m_ini.resize(Nm);
  m_ini = VectorXd::Constant(Nm, 1, 0);

  mesh.set_n_parameter(5);

  this->set_depth_weighting(this->depth_weighting_factor);

  constraint_z[0] = -1e8;
  constraint_z[1] = 1e8;
  constraint_x[0] = -1e8;
  constraint_x[1] = 1e8;
  constraint_y[0] = -1e8;
  constraint_y[1] = 1e8;
  // this->update_S_crg();
}

InversionBase::~InversionBase() {
  if (this->interpolator_m0 != NULL) {
    delete interpolator_m0;
  }
  if (this->interpolator_m0s != NULL) {
    delete interpolator_m0s;
  }
}

void InversionBase::set_mesh(const Mesh& mesh0) {
  this->mesh = mesh0;
  this->Nm = this->mesh.n_elems();
  mesh.set_n_parameter(5);
  m0 = VectorXd::Constant(Nm, 1, 0);
  cout << m0.rows() << endl;
  m.resize(Nm);
  m_min = VectorXd::Constant(Nm, 1, -1e6);
  m_max = VectorXd::Constant(Nm, 1, 1e6);
  m_ini = VectorXd::Constant(Nm, 1, 0);
  m0_s = VectorXd::Constant(Nm, 1, 0);

  this->init_matrices();
}

void InversionBase::set_S() {
  S_s.resize(Nm, Nm);
  S_s.reserve(Nm);
  S_s.setZero();

  S_x.resize(Nm, Nm);
  S_x.reserve(Nm);
  S_x.setZero();

  S_y.resize(Nm, Nm);
  S_y.reserve(Nm);
  S_y.setZero();

  S_z.resize(Nm, Nm);
  S_z.reserve(Nm);
  S_z.setZero();

  S_crg.resize(Nm, Nm);
  S_crg.reserve(Nm);
  S_crg.setZero();

  double dz, dx, dy;
  double z0, x0, y0;

  for (int i = 0; i < Nm; i++) {
    Cell* c = mesh.leaf_cells[i];
    c->get_size(dx, dy, dz);
    c->get_center(x0, y0, z0);

    S_s.coeffRef(i, i) += 1;
    S_x.coeffRef(i, i) += 1;
    S_y.coeffRef(i, i) += 1;
    S_z.coeffRef(i, i) += 1;
    S_crg.coeffRef(i, i) += 1;
    // S_s.coeffRef(i, i) += 1;
    // S_x.coeffRef(i, i) += r0 * dx;
    // S_y.coeffRef(i, i) += r0 * sin(x0) * dy;
    // S_z.coeffRef(i, i) += dr;

    // min_side=(min_side>dr)?dr:min_side;
    // min_side=(min_side>(r0 * dx))?(r0 * dx):min_side;
    // min_side=(min_side>(r0 * sin(x0) * dy))?(r0 * sin(x0) *
    // dy):min_side;
  }
  // S_x=S_x/min_side;
  // S_y=S_y/min_side;
  // S_z=S_z/min_side;
}

void InversionBase::init_matrices() {
  // This is necesseary because cell number is changed after the mesh is refines
  this->Nm = this->mesh.n_elems();
  V.resize(Nm, Nm);
  V.reserve(Nm);
  V.setZero();
  this->set_Vmatrix();
  V.makeCompressed();
  V.data().squeeze();

  Z.resize(Nm, Nm);
  Z.reserve(Nm);
  Z.setZero();
  this->set_depth_weighting(this->depth_weighting_factor);
  Z.makeCompressed();

  this->set_S();

  D_s.resize(Nm, Nm);
  D_s.reserve(Nm * 2);

  D_x1.resize(Nm, Nm);
  D_x1.reserve(Nm * 2);

  D_y1.resize(Nm, Nm);
  D_y1.reserve(Nm * 2);

  D_z1.resize(Nm, Nm);
  D_z1.reserve(Nm * 2);

  this->set_difference_matrix();

  D_s.makeCompressed();
  D_s.data().squeeze();
  D_x1.makeCompressed();
  D_x1.data().squeeze();
  D_y1.makeCompressed();
  D_y1.data().squeeze();
  D_z1.makeCompressed();
  D_z1.data().squeeze();
}

void InversionBase::set_difference_matrix() {
  D_s.setZero();
  D_x1.setZero();
  D_y1.setZero();
  D_z1.setZero();

  for (int i = 0; i < Nm; i++) {
    Cell* c = mesh.leaf_cells[i];
    int id0 = c->get_id();
    double dz0, dx0, dy0, dz1, dx1, dy1;
    double z0c, x0c, y0c;
    double z1c, x1c, y1c;
    c->get_center(x0c, y0c, z0c);
    c->get_size(dx0, dy0, dz0);

    D_s.coeffRef(id0, id0) += 1.0;

    Face* f = c->external_faces_x[1];
    if (f->isleaf) {
      Cell* neigh = f->neigh_cells[1];
      if (f->neigh_cells[0] != NULL && f->neigh_cells[1] != NULL) {
        // if()
        neigh->get_size(dx1, dy1, dz1);
        int id1 = neigh->get_id();
        assert(id0 != id1);
        D_x1.coeffRef(id0, id0) += -1.0 / (0.5 * (dx0 + dx1));
        D_x1.coeffRef(id0, id1) += 1.0 / (0.5 * (dx0 + dx1));
      } else {
        f = c->external_faces_x[0];
        if (f->isleaf) {
          neigh = f->neigh_cells[0];
          // if (neigh == NULL)
          // {
          //     neigh = f->neigh_cells[1];
          // }
          neigh->get_size(dx1, dy1, dz1);
          int id1 = neigh->get_id();
          assert(id0 != id1);
          D_x1.coeffRef(id0, id0) += 1.0 / (0.5 * (dx0 + dx1));
          D_x1.coeffRef(id0, id1) += -1.0 / (0.5 * (dx0 + dx1));
        } else {
          Cell* neigh = f->child_faces[0]->neigh_cells[0];
          neigh->get_size(dx1, dy1, dz1);
          D_x1.coeffRef(id0, id0) += 1.0 / (0.5 * (dx0 + dx1));
          for (int k = 0; k < 4; k++) {
            Cell* neigh = f->child_faces[k]->neigh_cells[0];
            int id1 = neigh->get_id();
            assert(id0 != id1);
            D_x1.coeffRef(id0, id1) += -1.0 / (0.5 * (dx0 + dx1)) * 0.25;
          }
        }
      }
    } else {
      Cell* neigh = f->child_faces[0]->neigh_cells[1];
      neigh->get_size(dx1, dy1, dz1);
      D_x1.coeffRef(id0, id0) += -1.0 / (0.5 * (dx0 + dx1));
      for (int k = 0; k < 4; k++) {
        Cell* neigh = f->child_faces[k]->neigh_cells[1];
        int id1 = neigh->get_id();
        assert(id0 != id1);
        D_x1.coeffRef(id0, id1) += 1.0 / (0.5 * (dx0 + dx1)) * 0.25;
      }
    }
    // y

    f = c->external_faces_y[1];
    if (f->isleaf) {
      Cell* neigh = f->neigh_cells[1];
      if (f->neigh_cells[0] != NULL && f->neigh_cells[1] != NULL) {
        int id1 = neigh->get_id();
        neigh->get_size(dx1, dy1, dz1);
        assert(id0 != id1);
        D_y1.coeffRef(id0, id0) += -1.0 / (0.5 * (dy0 + dy1));
        D_y1.coeffRef(id0, id1) += 1.0 / (0.5 * (dy0 + dy1));
      } else {
        f = c->external_faces_y[0];
        if (f->isleaf) {
          neigh = f->neigh_cells[0];
          // if (neigh == NULL)
          // {
          //     neigh = c->external_faces_y[0]->neigh_cells[1];
          // }
          neigh->get_size(dx1, dy1, dz1);
          int id1 = neigh->get_id();
          assert(id0 != id1);
          D_y1.coeffRef(id0, id0) += 1.0 / (0.5 * (dy0 + dy1));
          D_y1.coeffRef(id0, id1) += -1.0 / (0.5 * (dy0 + dy1));
        } else {
          Cell* neigh = f->child_faces[0]->neigh_cells[0];
          neigh->get_size(dx1, dy1, dz1);
          D_y1.coeffRef(id0, id0) += 1.0 / (0.5 * (dy0 + dy1));
          for (int k = 0; k < 4; k++) {
            Cell* neigh = f->child_faces[k]->neigh_cells[0];
            int id1 = neigh->get_id();
            assert(id0 != id1);
            D_y1.coeffRef(id0, id1) += -1.0 / (0.5 * (dy0 + dy1)) * 0.25;
          }
        }
      }
    } else {
      Cell* neigh = f->child_faces[0]->neigh_cells[1];
      neigh->get_size(dx1, dy1, dz1);
      D_y1.coeffRef(id0, id0) += -1.0 / (0.5 * (dy0 + dy1));
      for (int k = 0; k < 4; k++) {
        Cell* neigh = f->child_faces[k]->neigh_cells[1];
        int id1 = neigh->get_id();
        assert(id0 != id1);
        D_y1.coeffRef(id0, id1) += 1.0 / (0.5 * (dy0 + dy1)) * 0.25;
      }
    }

    // z
    f = c->external_faces_z[1];
    if (f->isleaf) {
      Cell* neigh = f->neigh_cells[1];
      if (f->neigh_cells[0] != NULL && f->neigh_cells[1] != NULL) {
        int id1 = neigh->get_id();
        neigh->get_size(dx1, dy1, dz1);
        if (id0 == id1) {
          cout << "c" << endl;
          cout << "x direction" << endl;
          for (int j = 0; j < 2; j++) {
            c->external_faces_x[j]->display();
          }
          cout << "y direction" << endl;
          for (int j = 0; j < 2; j++) {
            c->external_faces_y[j]->display();
          }
          cout << "z direction" << endl;
          for (int j = 0; j < 2; j++) {
            c->external_faces_z[j]->display();
          }
          cout << "-----------------------------" << endl;
          cout << "neigh 0" << endl;
          cout << "x direction" << endl;
          for (int j = 0; j < 2; j++) {
            f->neigh_cells[0]->external_faces_x[j]->display();
          }
          cout << "y direction" << endl;
          for (int j = 0; j < 2; j++) {
            f->neigh_cells[0]->external_faces_y[j]->display();
          }
          cout << "z direction" << endl;
          for (int j = 0; j < 2; j++) {
            f->neigh_cells[0]->external_faces_z[j]->display();
          }
          cout << "-----------------------------" << endl;
          cout << "neigh 1" << endl;
          cout << "x direction" << endl;
          for (int j = 0; j < 2; j++) {
            f->neigh_cells[1]->external_faces_x[j]->display();
          }
          cout << "y direction" << endl;
          for (int j = 0; j < 2; j++) {
            f->neigh_cells[1]->external_faces_y[j]->display();
          }
          cout << "z direction" << endl;
          for (int j = 0; j < 2; j++) {
            f->neigh_cells[1]->external_faces_z[j]->display();
          }
        }
        assert(id0 != id1);

        D_z1.coeffRef(id0, id0) += -1.0 / (0.5 * (dz0 + dz1));
        D_z1.coeffRef(id0, id1) += 1.0 / (0.5 * (dz0 + dz1));
      } else if (c->external_faces_z[0]->neigh_cells[0] != NULL) {
        f = c->external_faces_z[0];
        if (f->isleaf) {
          neigh = f->neigh_cells[0];
          neigh->get_size(dx1, dy1, dz1);
          int id1 = neigh->get_id();
          assert(id0 != id1);
          D_z1.coeffRef(id0, id0) += 1.0 / (0.5 * (dz0 + dz1));
          D_z1.coeffRef(id0, id1) += -1.0 / (0.5 * (dz0 + dz1));
        } else {
          Cell* neigh = f->child_faces[0]->neigh_cells[0];
          neigh->get_size(dx1, dy1, dz1);
          D_z1.coeffRef(id0, id0) += 1.0 / (0.5 * (dz0 + dz1));
          for (int k = 0; k < 4; k++) {
            Cell* neigh = f->child_faces[k]->neigh_cells[0];
            neigh->get_size(dx1, dy1, dz1);
            int id1 = neigh->get_id();
            assert(id0 != id1);
            D_z1.coeffRef(id0, id1) += -1.0 / (0.5 * (dz0 + dz1)) * 0.25;
          }
        }
      }
    } else {
      Cell* neigh = f->child_faces[0]->neigh_cells[1];
      neigh->get_size(dx1, dy1, dz1);
      D_z1.coeffRef(id0, id0) += -1.0 / (0.5 * (dz0 + dz1));
      for (int k = 0; k < 4; k++) {
        Cell* neigh = f->child_faces[k]->neigh_cells[1];
        neigh->get_size(dx1, dy1, dz1);
        int id1 = neigh->get_id();
        assert(id0 != id1);
        D_z1.coeffRef(id0, id1) += 1.0 / (0.5 * (dz0 + dz1)) * 0.25;
      }
    }
  }
}

void InversionBase::show_differece_matrix(unsigned int direction) {
  if (direction == NORTH_SOUTH) {
    cout << D_x1 << endl;
  } else if (direction == WEST_EAST) {
    cout << D_y1 << endl;
  } else if (direction == UP_DOWN) {
    cout << D_z1 << endl;
  }
}

void InversionBase::set_Vmatrix() {
  V.setZero();
  // V.setIdentity();
  assert(N_obs != 0);
  assert(Nm != 0);
  for (int i = 0; i < mesh.n_elems(); i++) {
    const RectPrism& t = mesh.get_elem(i);
    double v = t.get_volumn();
    V.coeffRef(i, i) += sqrt(v);
    // V.coeffRef(i, i) += 1;
  }
  // cout<<V<<endl;
}

void InversionBase::set_Tmatrix() {
  T_z.resize(Nm, Nm);
  T_z.setZero();
  T_x.resize(Nm, Nm);
  T_x.setZero();
  T_y.resize(Nm, Nm);
  T_y.setZero();
  VectorXd Dzs_v = D_z1 * m0_s;
  VectorXd Dxs_v = D_x1 * m0_s;
  VectorXd Dys_v = D_y1 * m0_s;

  SMatrix Dzs(Nm, Nm);
  SMatrix Dxs(Nm, Nm);
  SMatrix Dys(Nm, Nm);
  Dzs.setZero();
  Dxs.setZero();
  Dys.setZero();
  for (int i = 0; i < Nm; i++) {
    Dzs.coeffRef(i, i) = Dzs_v(i);
    Dxs.coeffRef(i, i) = Dxs_v(i);
    Dys.coeffRef(i, i) = Dys_v(i);
  }
  T_z = V * (Dys * D_x1 - Dxs * D_y1);
  T_x = V * (Dzs * D_y1 - Dys * D_z1);
  T_y = V * (Dxs * D_z1 - Dzs * D_x1);
}

void InversionBase::set_depth_weighting(double beta, int flag) {
  this->depth_weighting_factor = beta;
  Z.setZero();
  assert(N_obs != 0);  // Nd>=N_obs>=0
  assert(Nm != 0);
  if (flag == 0) {
    double zp = 0;
    for (int i = 0; i < N_obs; i++) {
      Point p = ob(i);
      zp += p.z();
    }
    zp = zp / N_obs;
    for (int i = 0; i < Nm; i++) {
      const RectPrism& t = mesh.get_elem(i);
      double z = 0.5 * (t._z[0] + t._z[1]);
      double weight = 1.0 / pow(std::abs(zp - z), beta * 0.5);
      Z.coeffRef(i, i) += weight;
    }
  }
  if (flag == 1) {
    for (int i = 0; i < Nm; i++) {
      double t = (G.col(i)).norm() / Nd;
      Z.coeffRef(i, i) += sqrt(t);
    }
  }
}

void InversionBase::set_dobs(const VectorXd& d) {
  int n = field_flag.count();
  assert(Nd == d.size());

  this->dobs = d;
  // cout<<dobs.rows()<<endl;
  Wd.resize(Nd, Nd);
  Wd.setZero();
  for (int i = 0; i < Nd; i++) {
    Wd.coeffRef(i, i) += 1.0;
  }
}

void InversionBase::set_dobs(const VectorXd& d,
                             double relative_error,
                             double a) {
  assert(relative_error < 0.7);
  int n = field_flag.count();
  assert(Nd == d.size());

  this->dobs = d;
  Wd.resize(Nd, Nd);
  Wd.setZero();
  for (int i = 0; i < Nd; i++) {
    Wd.coeffRef(i, i) += 1.0 / (relative_error * std::fabs(dobs(i)) + a);
  }
}
void InversionBase::set_dobs(const VectorXd& d,
                             vector<double> relative_error,
                             vector<double> a) {
  int n_fields = field_flag.count();
  assert(relative_error.size() == n_fields && a.size() == n_fields);
  assert(Nd == d.size());
  this->dobs = d;
  Wd.resize(Nd, Nd);
  Wd.setZero();
  for (int i = 0; i < n_fields; i++) {
    for (int j = 0; j < N_obs; j++) {
      Wd.coeffRef(j + i * N_obs, j + i * N_obs) +=
          1.0 / (relative_error[i] * std::fabs(dobs(j + i * N_obs)) + a[i]);
    }
  }
}

void InversionBase::set_Wd(const VectorXd& sigma) {
  assert(sigma.rows() == Wd.cols());
  Wd.setZero();
  for (int i = 0; i < Nd; i++) {
    Wd.coeffRef(i, i) += 1.0 / sigma(i);
  }
}

void InversionBase::set_observation(const Observation& ob0) {
  this->ob = ob0;
  N_obs = ob.get_n_obs();
  int n_fields = field_flag.count();
  Nd = N_obs * n_fields;

  Wd.resize(Nd, Nd);
  Wd.reserve(Nd);
  Wd.setIdentity();
}

void InversionBase::set_weights_of_objectives(double a_s,
                                              double a_z,
                                              double a_x,
                                              double a_y,
                                              double a_crg) {
  this->a_s = a_s;
  this->a_z = a_z;
  this->a_x = a_x;
  this->a_y = a_y;
  this->a_crg = a_crg;
}

void InversionBase::set_reference_model(VectorXd& m_ref) {
  this->m0 = m_ref;
}

void InversionBase::set_geometry_reference_model(VectorXd& m_ref) {
  this->m0_s = m_ref;
}
void InversionBase::set_m(VectorXd& m_) {
  this->m = m_;
}

void InversionBase::set_m_ini(VectorXd& m_ini_) {
  this->m_ini = m_ini_;
}

void InversionBase::set_min_max(VectorXd& m_min_, VectorXd& m_max_) {
  this->m_min = m_min_;
  this->m_max = m_max_;
}

void InversionBase::set_density_to_mesh() {
  assert(m.size() == Nm);
  for (int i = 0; i < Nm; i++) {
    Cell* t = mesh.leaf_cells[i];
    t->set_parameter(m(i), 0);
    // t.set_density(log10(Z.coeffRef(i,i)));
  }
}

void InversionBase::set_reference_model_to_mesh() {
  assert(m0.size() == Nm);
  for (int i = 0; i < Nm; i++) {
    Cell* t = mesh.leaf_cells[i];
    assert(t->parameters.size() > 1);
    t->set_parameter(m0(i), 1);
    t->set_parameter(m0_s(i), 2);
    // t.set_density(log10(Z.coeffRef(i,i)));
  }
}

void InversionBase::set_min_max_to_mesh() {
  assert(m.size() == Nm);
  for (int i = 0; i < Nm; i++) {
    Cell* t = mesh.leaf_cells[i];
    t->set_parameter(m_min(i), 3);
    t->set_parameter(m_max(i), 4);
    // t.set_density(log10(Z.coeffRef(i,i)));
  }
}

void InversionBase::set_petrophysics_constraint(
    VectorXd& m_ref_other_para,
    function<double(double)> relation) {
  int ns = m_ref_other_para.size();
  assert(ns == Nm);
  for (int i = 0; i < Nm; i++) {
    // if (fabs(m_ref_other_para(i)) < 1e-9)
    // {
    //     this->m0(i) = 0;
    // }
    // else
    // {
    this->m0(i) = relation(m_ref_other_para(i));
    // }
  }
}

VectorXd InversionBase::get_predicted_field() {
  VectorXd d_pre = (this->G) * (this->m);
  return d_pre;
}

void InversionBase::update_S_crg() {
  S_crg.setZero();
  for (int i_CELL = 0; i_CELL < Nm; i_CELL++) {
    double zc, xc, yc;
    mesh.leaf_cells[i_CELL]->get_center(xc, yc, zc);

    if ((zc < constraint_z[1] || std::abs(zc - constraint_z[1]) < 1e-10) &&
        (zc > constraint_z[0] || std::abs(zc - constraint_z[0]) < 1e-10) &&
        (xc < constraint_x[1] || std::abs(xc - constraint_x[1]) < 1e-10) &&
        (xc > constraint_x[0] || std::abs(xc - constraint_x[0]) < 1e-10) &&
        (yc < constraint_y[1] || std::abs(yc - constraint_y[1]) < 1e-10) &&
        (yc > constraint_y[0] || std::abs(yc - constraint_y[0]) < 1e-10)) {
      S_crg.coeffRef(i_CELL, i_CELL) += 1;
    }
  }
}

void InversionBase::output_predicted_data(string out_name) {
  VectorXd d_pre = this->get_predicted_field();
  this->out_data(d_pre, out_name);
}

void InversionBase::output_obs_data(string out_name) {
  // cout<<this->dobs.rows()<<endl;
  this->out_data(this->dobs, out_name);
}

void InversionBase::out_data(const VectorXd& d, string out_name) {
  vector<string> strs = {"V",    "g_z",  "g_x",  "g_y",  "T_zz",
                         "T_xz", "T_yz", "T_xx", "T_xy", "T_yy"};
  vector<unsigned int> field_label;
  for (unsigned int i = 0; i < strs.size(); i++) {
    if (field_flag[i]) {
      field_label.push_back(i);
    }
  }

  int n_com = field_flag.count();  // number of used components
  int n_ob = ob.get_n_obs();
  int nd = d.rows();
  // cout<<"xx  "<<d.rows()<<endl;

  // cout<<nd<<endl<<n_ob<<endl<<n_com<<endl;

  assert(nd % n_ob == 0);
  assert(nd / n_ob == n_com);

  for (int j = 0; j < n_com; j++) {
    string file_name = out_name + "_" + strs[field_label[j]];
    ofstream out_s(file_name);
    for (int i = 0; i < n_ob; i++) {
      const Point& p = ob(i);
      out_s << scientific;
      out_s << setw(30) << setprecision(15) << left << p.x() << setw(30) << left
            << p.y();
      out_s << scientific;
      out_s << setw(30) << setprecision(15) << left << d(i + j * n_ob) << endl;
    }
  }
}

void InversionBase::result2vtk(string filename) {
  this->set_density_to_mesh();
  this->set_reference_model_to_mesh();
  vector<string> parameter_name = {"model", "m0", "m0s"};
  mesh.out_model_vtk(filename + string(".vtk"), 3, parameter_name);
}

#ifdef USE_NETCDF
void InversionBase::result2netcdf(string filename) {
  this->set_density_to_mesh();
  this->mesh.out_model_netcdf(filename + string(".nc"));
}
#endif

void InversionBase::create_crg_model_from_data(string filename,
                                               int x_size,
                                               int y_size,
                                               int z_size,
                                               string data_order,
                                               int fast_dimension) {
  this->use_cross_gradient_constraint = true;
  // cout<<use_cross_gradient_constraint<<endl;
  if (this->interpolator_m0s != NULL) {
    delete interpolator_m0s;
    interpolator_m0s = NULL;
  }

  ifstream input_file;
  input_file.open(filename);
  assert(input_file.good());
  cout << "Read cross-gradient constraint model from " << filename << endl;

  vector<double> grid_x;  // x
  vector<double> grid_y;  // y
  vector<double> grid_z;  // z

  grid_x.resize(x_size);
  grid_y.resize(y_size);
  grid_z.resize(z_size);

  // the size of the grid in each dimension
  array<int, 3> grid_sizes;
  grid_sizes[0] = grid_x.size();
  grid_sizes[1] = grid_y.size();
  grid_sizes[2] = grid_z.size();

  int num_elements = grid_sizes[0] * grid_sizes[1] * grid_sizes[2];
  std::vector<double> f_values(num_elements);

  double x, y, z, val;  // each column

  if (fast_dimension == 0) {
    for (int k = 0; k < grid_z.size(); k++) {
      for (int j = 0; j < grid_y.size(); j++) {
        for (int i = 0; i < grid_x.size(); i++) {
          if (data_order == "yxz") {
            input_file >> y >> x >> z >> val;
          } else if (data_order == "xyz") {
            input_file >> x >> y >> z >> val;
          } else if (data_order == "zxy") {
            input_file >> z >> x >> y >> val;
          } else if (data_order == "zyx") {
            input_file >> z >> y >> x >> val;
          } else if (data_order == "yzx") {
            input_file >> y >> z >> x >> val;
          } else if (data_order == "xzy") {
            input_file >> x >> z >> y >> val;
          } else {
            cout << "data_order should be one of the following:" << endl;
            cout << "xyz yxz zxy zyx xzy yzx" << endl;
            std::abort();
          }

          if (k == 0 && j == 0) {
            grid_x[i] = x;  // x
          }

          f_values[i * grid_sizes[1] * grid_sizes[2] + j * grid_sizes[2] + k] =
              val;
        }
        if (k == 0) {
          grid_y[j] = y;  // y
        }
      }
      grid_z[k] = z;  // z
    }
  } else if (fast_dimension == 1) {
    for (int k = 0; k < grid_z.size(); k++) {
      for (int i = 0; i < grid_x.size(); i++) {
        for (int j = 0; j < grid_y.size(); j++) {
          if (data_order == "yxz") {
            input_file >> y >> x >> z >> val;
          } else if (data_order == "xyz") {
            input_file >> x >> y >> z >> val;
          } else if (data_order == "zxy") {
            input_file >> z >> x >> y >> val;
          } else if (data_order == "zyx") {
            input_file >> z >> y >> x >> val;
          } else if (data_order == "yzx") {
            input_file >> y >> z >> x >> val;
          } else if (data_order == "xzy") {
            input_file >> x >> z >> y >> val;
          } else {
            cout << "data_order should be one of the following:" << endl;
            cout << "xyz yxz zxy zyx xzy yzx" << endl;
            std::abort();
          }

          if (k == 0 && i == 0) {
            grid_y[j] = y;  // y
          }

          f_values[i * grid_sizes[1] * grid_sizes[2] + j * grid_sizes[2] + k] =
              val;
        }
        if (k == 0) {
          grid_x[i] = x;  // x
        }
      }
      grid_z[k] = z;  // z
    }
  } else {
    cout << "Parameter fast_dimension"
         << " should be 0 or 1" << endl;
  }
  input_file.close();
  input_file.clear();

  // construct the grid in each dimension.
  // note that we will pass in a sequence of iterators pointing to the beginning
  // of each grid
  std::vector<std::vector<double>::iterator> grid_iter_list;
  grid_iter_list.push_back(grid_x.begin());
  grid_iter_list.push_back(grid_y.begin());
  grid_iter_list.push_back(grid_z.begin());

  // construct the interpolator. the last two arguments are pointers to the
  // underlying data
  this->interpolator_m0s = new InterpMultilinear<3, double>(
      grid_iter_list.begin(), grid_sizes.begin(), f_values.data(),
      f_values.data() + num_elements);
  // interpolator_m0s
  Mesh mesh_crg;
  double x_space =
      (grid_x[grid_x.size() - 1] - grid_x[0]) / (grid_x.size() - 1);
  double y_space =
      (grid_y[grid_y.size() - 1] - grid_y[0]) / (grid_y.size() - 1);
  double z_space =
      (grid_z[grid_z.size() - 1] - grid_z[0]) / (grid_z.size() - 1);  // km
  double x_model[2] = {grid_x[0] - 0.5 * x_space,
                       grid_x[grid_x.size() - 1] + 0.5 * x_space};
  double y_model[2] = {grid_y[0] - 0.5 * y_space,
                       grid_y[grid_y.size() - 1] + 0.5 * y_space};
  double z_model[2] = {grid_z[0] - 0.5 * z_space,
                       grid_z[grid_z.size() - 1] + 0.5 * z_space};
  mesh_crg.generate_regular_mesh(x_model, x_size, y_model, y_size, z_model,
                                 z_size);

  for (int i = 0; i < mesh_crg.n_elems(); i++) {
    Cell* c = mesh_crg.leaf_cells[i];
    double xc, yc, zc;
    c->get_center(xc, yc, zc);

    array<double, 3> args = {xc, yc, zc};
    double val = (*interpolator_m0s).interp(args.begin());
    c->set_parameter(val);
  }
  mesh_crg.out_model_vtk("crg_model.vtk", 1,
                         vector<string>(1,"crg"));

#ifdef USE_NETCDF
  mesh_crg.out_model_netcdf(string("crg_model.nc"), 0, "crg", "");
#endif  
}

void InversionBase::create_ref_model_from_data(string filename,
                                               int x_size,
                                               int y_size,
                                               int z_size,
                                               string data_order,
                                               int fast_dimension) {
  this->use_petrophysical_constraint = true;
  if (this->interpolator_m0 != NULL) {
    delete interpolator_m0;
    interpolator_m0 = NULL;
  }

  ifstream input_file;
  input_file.open(filename);
  assert(input_file.good());
  cout << "Read converted density model from " << filename << endl;

  vector<double> grid_x;  // x
  vector<double> grid_y;  // y
  vector<double> grid_z;  // z

  grid_x.resize(x_size);
  grid_y.resize(y_size);
  grid_z.resize(z_size);

  // the size of the grid in each dimension
  array<int, 3> grid_sizes;
  grid_sizes[0] = grid_x.size();
  grid_sizes[1] = grid_y.size();
  grid_sizes[2] = grid_z.size();

  int num_elements = grid_sizes[0] * grid_sizes[1] * grid_sizes[2];
  std::vector<double> f_values(num_elements);

  double x, y, z, val;  // each column

  if (fast_dimension == 0) {
    for (int k = 0; k < grid_z.size(); k++) {
      for (int j = 0; j < grid_y.size(); j++) {
        for (int i = 0; i < grid_x.size(); i++) {
          if (data_order == "yxz") {
            input_file >> y >> x >> z >> val;
          } else if (data_order == "xyz") {
            input_file >> x >> y >> z >> val;
          } else if (data_order == "zxy") {
            input_file >> z >> x >> y >> val;
          } else if (data_order == "zyx") {
            input_file >> z >> y >> x >> val;
          } else if (data_order == "yzx") {
            input_file >> y >> z >> x >> val;
          } else if (data_order == "xzy") {
            input_file >> x >> z >> y >> val;
          } else {
            cout << "data_order should be one of the following:" << endl;
            cout << "xyz yxz zxy zyx xzy yzx" << endl;
            std::abort();
          }

          if (k == 0 && j == 0) {
            grid_x[i] = x;  // x
          }

          f_values[i * grid_sizes[1] * grid_sizes[2] + j * grid_sizes[2] + k] =
              val;
        }
        if (k == 0) {
          grid_y[j] = y;  // y
        }
      }
      grid_z[k] = z;  // z
    }
  } else if (fast_dimension == 1) {
    for (int k = 0; k < grid_z.size(); k++) {
      for (int i = 0; i < grid_x.size(); i++) {
        for (int j = 0; j < grid_y.size(); j++) {
          if (data_order == "yxz") {
            input_file >> y >> x >> z >> val;
          } else if (data_order == "xyz") {
            input_file >> x >> y >> z >> val;
          } else if (data_order == "zxy") {
            input_file >> z >> x >> y >> val;
          } else if (data_order == "zyx") {
            input_file >> z >> y >> x >> val;
          } else if (data_order == "yzx") {
            input_file >> y >> z >> x >> val;
          } else if (data_order == "xzy") {
            input_file >> x >> z >> y >> val;
          } else {
            cout << "data_order should be one of the following:" << endl;
            cout << "xyz yxz zxy zyx xzy yzx" << endl;
            std::abort();
          }

          if (k == 0 && i == 0) {
            grid_y[j] = y;  // y
          }

          f_values[i * grid_sizes[1] * grid_sizes[2] + j * grid_sizes[2] + k] =
              val;
        }
        if (k == 0) {
          grid_x[i] = x;  // x
        }
      }
      grid_z[k] = z;  // z
    }
  } else {
    cout << "Parameter fast_dimension"
         << " should be 0 or 1" << endl;
  }
  input_file.close();
  input_file.clear();

  // construct the grid in each dimension.
  // note that we will pass in a sequence of iterators pointing to the beginning
  // of each grid
  std::vector<std::vector<double>::iterator> grid_iter_list;
  grid_iter_list.push_back(grid_x.begin());
  grid_iter_list.push_back(grid_y.begin());
  grid_iter_list.push_back(grid_z.begin());

  // construct the interpolator. the last two arguments are pointers to the
  // underlying data
  this->interpolator_m0 = new InterpMultilinear<3, double>(
      grid_iter_list.begin(), grid_sizes.begin(), f_values.data(),
      f_values.data() + num_elements);

  Mesh mesh_ref;
  double x_space =
      (grid_x[grid_x.size() - 1] - grid_x[0]) / (grid_x.size() - 1);
  double y_space =
      (grid_y[grid_y.size() - 1] - grid_y[0]) / (grid_y.size() - 1);
  double z_space =
      (grid_z[grid_z.size() - 1] - grid_z[0]) / (grid_z.size() - 1);  // km
  double x_model[2] = {grid_x[0] - 0.5 * x_space,
                       grid_x[grid_x.size() - 1] + 0.5 * x_space};
  double y_model[2] = {grid_y[0] - 0.5 * y_space,
                       grid_y[grid_y.size() - 1] + 0.5 * y_space};
  double z_model[2] = {(grid_z[0] - 0.5 * z_space),
                       (grid_z[grid_z.size() - 1] + 0.5 * z_space)};
  mesh_ref.generate_regular_mesh(x_model, x_size, y_model, y_size, z_model,
                                 z_size);

  for (int i = 0; i < mesh_ref.n_elems(); i++) {
    Cell* c = mesh_ref.leaf_cells[i];
    double xc, yc, zc;
    c->get_center(xc, yc, zc);

    array<double, 3> args = {xc, yc, zc};
    double val = (*interpolator_m0).interp(args.begin());
    c->set_parameter(val);
  }
  mesh_ref.out_model_vtk("ref_model.vtk", 1,
                         vector<string>(1, "reference_model"));

#ifdef USE_NETCDF
  mesh_ref.out_model_netcdf(string("ref_model.nc"), 0, "ref", "");
#endif
}