#include "GaussNewtonInversion.h"
void GaussNewtonInversion::invert() {
  // set reference model using interpolation
  if (use_petrophysical_constraint == true) {
    for (int i = 0; i < Nm; i++) {
      Cell* c = this->mesh.leaf_cells[i];
      double xc, yc, zc;
      c->get_center(xc, yc, zc);
      array<double, 3> args = {xc, yc, zc};
      double val = (*interpolator_m0).interp(args.begin());
      c->set_parameter(val, 1);
    }
    // cout << "m0" << endl;
    this->m0.resize(Nm);
    mesh.get_model_parameter_from_mesh(m0, 1);
  }

  if (use_cross_gradient_constraint == true) {
    cout << "Cross gradient constraint is used" << endl;
    for (int i = 0; i < Nm; i++) {
      Cell* c = this->mesh.leaf_cells[i];
      double xc, yc, zc;
      c->get_center(xc, yc, zc);
      array<double, 3> args = {xc, yc, zc};
      double val = (*interpolator_m0s).interp(args.begin());
      c->set_parameter(val, 2);
    }
    this->m0_s.resize(Nm);
    mesh.get_model_parameter_from_mesh(m0_s, 2);
    this->set_Tmatrix();
    // update_S_crg();
  }

  SMatrix WL1_s(Nm, Nm);
  SMatrix WL1_x(Nm, Nm);
  SMatrix WL1_y(Nm, Nm);
  SMatrix WL1_z(Nm, Nm);
  WL1_s.setIdentity();
  WL1_x.setIdentity();
  WL1_y.setIdentity();
  WL1_z.setIdentity();

  SMatrix Ws = a_s * S_s * V * D_s * Z;
  SMatrix Wx = a_x * S_x * V * D_x1 * Z;
  SMatrix Wy = a_y * S_y * V * D_y1 * Z;
  SMatrix Wz = a_z * S_z * V * D_z1 * Z;
  Ws.makeCompressed();
  Wx.makeCompressed();
  Wy.makeCompressed();
  Wz.makeCompressed();

  VectorXd Wsm_m0;
  VectorXd Wxm;
  VectorXd Wym;
  VectorXd Wzm;

  VectorXd Ws_m0 = Ws * m0;

  int nrow =
      (use_cross_gradient_constraint == true) ? (Nd + 7 * Nm) : (Nd + 4 * Nm);
  SMatrix A(nrow, Nm);
  A.reserve(Nd * Nm + Ws.nonZeros() + Wx.nonZeros() + Wy.nonZeros() +
            Wz.nonZeros() + 6 * Nm);
  VectorXd b(nrow);
  A.setZero();
  b.setZero();
  // A.topRows(Nd) = (Wd * G).sparseView();
  b.head(Nd) = Wd * dobs;
  // cout<<"x"<<Nm<<endl;
  SMatrix P(Nm, Nm);
  P.setZero();
  P.reserve(Nm);
  P.makeCompressed();

  double lambda = this->max_lambda;

  VectorXd delta_x(Nm), m_trial(Nm);  // m_trial_last(Nm), m_trial_last2(Nm);
  double misfit, misfit_last1, misfit_last2, misfit_last_iteration;

  
  if (use_petrophysical_constraint) {
    m = m0;
    for (int i = 0; i < Nm; i++) {
      if (m(i) > m_max(i) || abs(m(i) - m_max(i)) < 1e-7) {
        m(i) = m_max(i) - 0.01;
      }
      if (m(i) < m_min(i) || abs(m(i) - m_min(i)) < 1e-7) {
        m(i) = m_min(i) + 0.01;
      }
    }
  } else {
    m = m_ini;
    for (int i = 0; i < Nm; i++) {
      if (m(i) > m_max(i) || abs(m(i) - m_max(i)) < 1e-7) {
        m(i) = 0.5*(m_max(i)+m_min(i));
      }
      if (m(i) < m_min(i) || abs(m(i) - m_min(i)) < 1e-7) {
        m(i) = 0.5*(m_max(i)+m_min(i));
      }
    }
  }

  misfit = (Wd * (G * m - dobs)).squaredNorm() / Nd;

  misfit_last_iteration = 5 * misfit;

  for (int id = 0; id < Nm; id++) {
    P.coeffRef(id, id) =
        (m_max(id) - m(id)) * (m(id) - m_min(id)) / (m_max(id) - m_min(id));
  }
  LeastSquaresConjugateGradient<SMatrix> lscg;
  lscg.setTolerance(cg_tol);

  ofstream out_lambda_misfit("lambda_misfit_GN");
  out_lambda_misfit << setw(25) << "lambda" << setw(25) << "misfit" << endl;

  ofstream out_iteration_misfit("Iteration_misfit_GN");
  out_iteration_misfit << setw(25) << "Iteranation number" << setw(25)
                       << "number of tried regulartion parameters" << setw(25)
                       << "misfit" << endl;

  int i, j;
  double lambda_opt;
  double misfit_opt;
  VectorXd m_opt;
  misfit_opt = misfit;
  m_opt = m;

  int c_refinement = 0;
  int i_between_refinement = 0;

  int nnz =
      (use_cross_gradient_constraint == true)
          ? (Nd * Nm + Ws.nonZeros() + Wx.nonZeros() + Wy.nonZeros() +
             Wz.nonZeros() + T_x.nonZeros() + T_y.nonZeros() + T_z.nonZeros())
          : (Nd * Nm + Ws.nonZeros() + Wx.nonZeros() + Wy.nonZeros() +
             Wz.nonZeros());
  A.reserve(nnz);
  cout << "Set G to A" << endl;
  for (size_t r_id = 0; r_id < Nd; r_id++) {
    for (size_t c_id = 0; c_id < Nm; c_id++) {
      A.insert(r_id, c_id) = Wd.coeffRef(r_id, r_id) * G(r_id, c_id);
    }
  }
  b.head(Nd) = Wd * dobs;

  for (i = 0; i < GN_iter; i++) {
    cout << "The " << i + 1 << " th"
         << " Gauss-Newton Iteration:" << endl;
    cout << "Number of elements: " << Nm << endl;

    double min_cell_dx, min_cell_dy, min_cell_dz;
    int max_lev;

    mesh.get_minimum_size(min_cell_dx, min_cell_dy, min_cell_dz, max_lev);
    cout << "The smallest cell: (dx= " << min_cell_dx << " m, "
         << "dy= " << min_cell_dy << " m, "
         << "dz= " << min_cell_dz << " m)" << endl;

    m_trial = m;
    misfit_last1 = 10 * misfit;   // misfit for last lambda
    misfit_last2 = 100 * misfit;  // misfit for the lambda before last lambda
    misfit_last_iteration = misfit;
    lambda = this->max_lambda;

    A.middleRows(Nd, Nm) = sqrt(lambda) * WL1_s * Ws;
    A.middleRows(Nd + Nm, Nm) = sqrt(lambda) * WL1_z * Wz;
    A.middleRows(Nd + 2 * Nm, Nm) = sqrt(lambda) * WL1_x * Wx;
    A.middleRows(Nd + 3 * Nm, Nm) = sqrt(lambda) * WL1_y * Wy;
    // cout << S_crg.rows() << endl;
    // cout << S_crg.cols() << endl;
    if (use_cross_gradient_constraint == true) {
      A.middleRows(Nd + 4 * Nm, Nm) = this->a_crg * S_crg * T_z;
      A.middleRows(Nd + 5 * Nm, Nm) = this->a_crg * S_crg * T_x;
      A.middleRows(Nd + 6 * Nm, Nm) = this->a_crg * S_crg * T_y;
    }
    A.makeCompressed();

    b.segment(Nd, Nm) = sqrt(lambda) * WL1_s * Ws_m0;
    // b.segment(Nd + Nm, Nm) = sqrt(lambda) * Wr_m0;
    // b.segment(Nd + 2 * Nm, Nm) = sqrt(lambda) * Wtheta_m0;
    // b.segment(Nd + 3 * Nm, Nm) = sqrt(lambda) * Wphi_m0;

    if (cg_iteration_factor < 1) {
      lscg.setMaxIterations(cg_iteration_factor * Nm);
    } else {
      lscg.setMaxIterations(int(cg_iteration_factor));
    }
    for (int id = 0; id < Nm; id++) {
      int n = 1.0;
      P.coeffRef(id, id) = n * (m_max(id) - m(id)) * (m(id) - m_min(id)) /
                           (m_max(id) - m_min(id));
    }
    delta_x.setZero();
    for (j = 0; j < n_lambda; j++) {
      lscg.compute(A * P);
      // delta_x = lscg.solveWithGuess(b - A * m, delta_x);
      delta_x = lscg.solve(b - A * m);

      // m_trial_last2 = m_trial_last;
      // m_trial_last = m_trial;
      double temp, temp2, temp3;
      // update
      for (int id = 0; id < Nm; id++) {
        temp = exp(delta_x(id));
        temp2 = (m_min(id) * (m_max(id) - m(id)) +
                 m_max(id) * (m(id) - m_min(id)) * temp);
        temp3 = ((m_max(id) - m(id)) + (m(id) - m_min(id)) * temp);
        if (std::isinf(temp) || std::isinf(temp2)) {
          m_trial(id) = (m_max(id) - 1e-6 + m_trial(id)) / 2.0;
          if (id == 0) {
            m_trial(id) = 0.5 * (m_trial(id) + m_trial(id + 1));
          } else {
            m_trial(id) = 0.5 * (m_trial(id) + m_trial(id - 1));
          }
        } else if (fabs(temp3) < 1e-10) {
          m_trial(id) = (m_min(id) + 1e-6 + m_trial(id)) / 2.0;
          if (id == 0) {
            m_trial(id) = 0.5 * (m_trial(id) + m_trial(id + 1));
          } else {
            m_trial(id) = 0.5 * (m_trial(id) + m_trial(id - 1));
          }
        } else {
          m_trial(id) = temp2 / temp3;
        }
      }

      misfit_last2 = misfit_last1;
      misfit_last1 = misfit;
      misfit = (Wd * (G * m_trial - dobs)).squaredNorm() / Nd;
      // cout << "  Lambda=" << lambda << ", ";
      std::cout << "  #" << j << " lambda=" << lambda << ", misfit=" << misfit
                << " (CG iterations: " << lscg.iterations() << ","
                << " error: " << lscg.error() << ")" << std::endl;
      out_lambda_misfit << scientific << setw(25) << setprecision(14) << lambda
                        << scientific << setw(25) << setprecision(14)
                        << (misfit) << endl;

      if (j == 0 || misfit < misfit_opt) {
        misfit_opt = misfit;
        m_opt = m_trial;
        lambda_opt = lambda;
      }
      if ((std::abs(misfit - misfit_last1) / min(misfit, misfit_last1) <
           stag_tol) &&
          (std::abs(misfit_last1 - misfit_last2) /
               min(misfit_last1, misfit_last2) <
           stag_tol) &&
          j > 1) {
        cout << "  Misfit stagnates, go to next Gauss-Newton iteration."
             << endl;
        break;
      }
      if (misfit > misfit_last1 && misfit_last1 > misfit_last2 && j > 1) {
        cout
            << "  Misfit starts to increase, go to next Gauss-Newton iteration."
            << endl;
        break;
      } else {
        lambda = lambda * lambda_decreasing_rate;

        Wsm_m0 = Ws * (m_opt - m0);
        Wxm = Wx * m_opt;
        Wym = Wy * m_opt;
        Wzm = Wz * m_opt;

        for (int i = 0; i < Nm; i++) {
          WL1_s.coeffRef(i, i) =
              pow(Wsm_m0(i) * Wsm_m0(i) + epsilon2, -(2 - Lp) / 4.0);
          WL1_x.coeffRef(i, i) =
              pow(Wxm(i) * Wxm(i) + epsilon2, -(2 - Lp) / 4.0);
          WL1_y.coeffRef(i, i) =
              pow(Wym(i) * Wym(i) + epsilon2, -(2 - Lp) / 4.0);
          WL1_z.coeffRef(i, i) =
              pow(Wzm(i) * Wzm(i) + epsilon2, -(2 - Lp) / 4.0);
        }

        A.middleRows(Nd, Nm) = sqrt(lambda) * WL1_s * Ws;
        A.middleRows(Nd + Nm, Nm) = sqrt(lambda) * WL1_z * Wz;
        A.middleRows(Nd + 2 * Nm, Nm) = sqrt(lambda) * WL1_x * Wx;
        A.middleRows(Nd + 3 * Nm, Nm) = sqrt(lambda) * WL1_y * Wy;
        if (use_cross_gradient_constraint == true) {
          A.middleRows(Nd + 4 * Nm, Nm) = this->a_crg * S_crg * T_z;
          A.middleRows(Nd + 5 * Nm, Nm) = this->a_crg * S_crg * T_x;
          A.middleRows(Nd + 6 * Nm, Nm) = this->a_crg * S_crg * T_y;
        }

        b.segment(Nd, Nm) = sqrt(lambda) * WL1_s * Ws_m0;
      }
      if (misfit < target_misfit) {
        break;
      }
    }
    // lambda=lambda_opt;
    m = m_opt;
    misfit = misfit_opt;

    out_iteration_misfit << fixed << setw(25) << i + 1 << fixed << setw(25)
                         << ((j == n_lambda) ? n_lambda : (j + 1)) << scientific
                         << setw(25) << setprecision(14) << misfit << endl;

    if (record_process) {
      this->set_density_to_mesh();
      this->set_reference_model_to_mesh();

#ifdef USE_NETCDF
      if (use_cross_gradient_constraint == true) {
        this->mesh.out_model_netcdf(
            string("Structural_constraint_at_") + to_string(i) + string(".nc"),
            2, "crg", "");
      }
      if (use_petrophysical_constraint == true) {
        this->mesh.out_model_netcdf(string("Converted_density_model_at_") +
                                        to_string(i) + string(".nc"),
                                    1, "ref", "");
      }

      this->mesh.out_model_netcdf(string("result_at_") + to_string(i) +
                                  string(".nc"));
#endif
      this->mesh.out_model_vtk(string("result_at_") + to_string(i) +
                               string(".vtk"));
    }

    if (misfit < target_misfit || i == GN_iter - 1) {
      if (misfit < target_misfit) {
        cout << "Stop iteration, because the target misift has been achieved"
             << endl;
      } else if (i == GN_iter - 1) {
        cout << "Stop iteration, because the maximum iteration has been reached"
             << endl;
      }
      cout << "Gauss-Newton iteration number:" << i + 1 << endl;
      break;
    }

    bool flag1 = abs(misfit - misfit_last_iteration) /
                     min(misfit, misfit_last_iteration) <
                 stag_tol;
    if (c_refinement >= max_refinement_number && flag1) {
      cout << "Stop iteration, because the misfit stagnated" << endl;
      cout << "Gauss-Newton iteration number:" << i + 1 << endl;
      break;
    }

    if ((c_refinement < max_refinement_number) &&
        (++i_between_refinement % interval_between_refinements == 0 || flag1)) {
      cout << "\n";
      i_between_refinement = 0;
      ++c_refinement;
      if (flag1) {
        cout << "Misfit stagnates, refine mesh ..." << endl;
      } else {
        cout << interval_between_refinements
             << ((i_between_refinement == 1) ? string(" iteration has")
                                             : string(" iterations have"))
             << " passed since last refinement, refine mesh ..." << endl;
      }
      cout << c_refinement
           << (c_refinement == 1
                   ? ("st")
                   : (c_refinement == 2
                          ? ("nd")
                          : (c_refinement == 3 ? ("rd") : ("th"))))
           << " refinement." << endl
           << endl;

      if (use_cross_gradient_constraint && use_petrophysical_constraint) {
        refine_mesh(refinement_percentage, *interpolator_m0s, *interpolator_m0);
        // cout<<"stop here"<<endl;
        // abort();
      } else if (use_cross_gradient_constraint &&
                 (!use_petrophysical_constraint)) {
        refine_mesh(refinement_percentage, *interpolator_m0s, "crg");
        // cout<<"2"<<endl;
      } else if (use_petrophysical_constraint &&
                 (!use_cross_gradient_constraint)) {
        refine_mesh(refinement_percentage, *interpolator_m0, "pet");
        // cout << "2" << endl;
      } else {
        refine_mesh(refinement_percentage);
        // cout<<"2"<<endl;
      }

      cout << "Finished refinement" << endl;

      cout << "Initialize Ws" << endl;
      Ws.resize(Nm, Nm);
      Ws.reserve(2 * Nm);
      Ws = a_s * S_s * V * D_s * Z;
      Ws.makeCompressed();
      Ws.data().squeeze();

      cout << "Initialize Wx" << endl;
      Wx.resize(Nm, Nm);
      Wx.reserve(2 * Nm);
      Wx = a_x * S_x * V * D_x1 * Z;
      Wx.makeCompressed();
      Wx.data().squeeze();

      cout << "Initialize Wy" << endl;
      Wy.resize(Nm, Nm);
      Wy.reserve(2 * Nm);
      Wy = a_y * S_y * V * D_y1 * Z;
      Wy.makeCompressed();
      Wy.data().squeeze();

      cout << "Initialize Wz" << endl;
      Wz.resize(Nm, Nm);
      Wz.reserve(2 * Nm);
      Wz = a_z * S_z * V * D_z1 * Z;
      Wz.makeCompressed();
      Wz.data().squeeze();

      // cout<<"1"<<endl;
      WL1_s.resize(Nm, Nm);
      WL1_s.reserve(Nm);
      WL1_x.resize(Nm, Nm);
      WL1_x.reserve(Nm);
      WL1_y.resize(Nm, Nm);
      WL1_y.reserve(Nm);
      WL1_z.resize(Nm, Nm);
      WL1_z.reserve(Nm);
      WL1_s.setIdentity();
      WL1_x.setIdentity();
      WL1_y.setIdentity();
      WL1_z.setIdentity();

      Wsm_m0.resize(Nm);
      Wxm.resize(Nm);
      Wym.resize(Nm);
      Wzm.resize(Nm);
      Wsm_m0 = Ws * (m - m0);
      Wxm = Wx * m;
      Wym = Wy * m;
      Wzm = Wz * m;

      for (int i = 0; i < Nm; i++) {
        WL1_s.coeffRef(i, i) =
            pow(Wsm_m0(i) * Wsm_m0(i) + epsilon2, -(2 - Lp) / 4.0);
        WL1_x.coeffRef(i, i) = pow(Wxm(i) * Wxm(i) + epsilon2, -(2 - Lp) / 4.0);
        WL1_y.coeffRef(i, i) = pow(Wym(i) * Wym(i) + epsilon2, -(2 - Lp) / 4.0);
        WL1_z.coeffRef(i, i) = pow(Wzm(i) * Wzm(i) + epsilon2, -(2 - Lp) / 4.0);
      }

      // update_S_crg();
      // cout<<S_crg<<endl;

      Ws_m0 = Ws * m0;
      if (use_cross_gradient_constraint == true) {
        this->set_Tmatrix();
      }

      nrow = (use_cross_gradient_constraint == true) ? (Nd + 7 * Nm)
                                                     : (Nd + 4 * Nm);
      A.resize(nrow, Nm);
      A.data().squeeze();
      nnz = (use_cross_gradient_constraint == true)
                ? (Nd * Nm + Ws.nonZeros() + Wx.nonZeros() + Wy.nonZeros() +
                   Wz.nonZeros() + T_x.nonZeros() + T_y.nonZeros() +
                   T_z.nonZeros())
                : (Nd * Nm + Ws.nonZeros() + Wx.nonZeros() + Wy.nonZeros() +
                   Wz.nonZeros());
      A.reserve(nnz);
      cout << "Set G to A" << endl;
      for (size_t r_id = 0; r_id < Nd; r_id++) {
        for (size_t c_id = 0; c_id < Nm; c_id++) {
          A.insert(r_id, c_id) = Wd.coeffRef(r_id, r_id) * G(r_id, c_id);
        }
      }

      b.resize(nrow, 1);
      b.setZero();
      b.head(Nd) = Wd * dobs;

      P.resize(Nm, Nm);
      P.setZero();
      P.reserve(Nm);
      P.data().squeeze();

      delta_x.resize(Nm, 1);
      m_trial.resize(Nm, 1);
    }
    // cout << i << endl;
  }
  if (i == GN_iter) {
    cout << "Gauss-Newton iteration number:" << GN_iter << endl;
  }

  this->set_density_to_mesh();

  final_lambda = lambda;
  final_misfit = misfit;
  cout << "Lambda=" << lambda << ", ";
  cout << "Misfit=" << misfit << endl;
}

void GaussNewtonInversion::set_CG_parameter(double cg_tol_,
                                            double cg_iteration_factor_) {
  cg_tol = cg_tol_;
  cg_iteration_factor_ = abs(cg_iteration_factor_);
  this->cg_iteration_factor = cg_iteration_factor_;
}

void GaussNewtonInversion::display_inversion_parameters() const {
  display_info_fields();
  cout << endl;

  cout << "alpha_s=" << a_s << endl;
  cout << "alpha_z=" << a_z << endl;
  cout << "alpha_x=" << a_x << endl;
  cout << "alpha_y=" << a_y << endl;

  cout << "alpha_crg=" << a_crg << endl;

  cout << endl;
  cout << "Maximum regularizaton parameter is " << max_lambda << endl;
  cout << "Maximum number of regularization parameters (lambda) used for "
          "trials: "
       << n_lambda << endl
       << endl;

  cout << "Target misfit is " << target_misfit << endl;

  cout << "Stagnation factor is " << stag_tol
       << " (the inversion stagnate when the relative difference of misfits "
          "at "
          "2 consecutive iterations is smaller than this factor)"
       << endl;

  cout << endl;
  if (cg_iteration_factor < 1) {
    cout << "Maximum number of iterations of the conjugate gradient method is "
         << cg_iteration_factor * 100 << "% of the element number" << endl;
  } else {
    cout << "Maximum number of iterations of the conjugate gradient method is "
         << cg_iteration_factor << endl;
  }
  cout << "Tolerance value for the conjugate gradient method is " << cg_tol
       << endl;
  cout << endl;

  cout << "Maximum number of Gauss-Newton iterations: " << GN_iter << endl;

  cout << "Model density contrast limits: "
       << "[" << m_min.minCoeff() << ", " << m_max.maxCoeff() << "]"
       << "kg/m3" << endl;

  if (max_refinement_number > 0) {
    cout << endl;
    cout << "The inversion mesh will be adaptively refined at every "
         << interval_between_refinements << " iteration." << endl;
    cout << "The mesh may be refined for " << max_refinement_number
         << " times AT MOST!" << endl;
    cout << "For each time of refinement, " << refinement_percentage * 100
         << "% of cells are marked for refinement. BUT the exact number of "
            "cells that are refined might be different from this number "
            "because the "
            "refinement conditions imposed in the algorithm (i.e. cells "
            "smaller than the specified minimum "
            "size cannot be refined; the maximum difference in the levels of "
            "adjacent elements in "
            "the tree must not exceed 1)."
         << endl
         << endl;
    cout << "Cells smaller than (dx=" << min_dx << "m, dy=" << min_dy
         << "m, dz=" << min_dz << "m) cannot be refined" << endl;
  } else {
    cout << endl;
    cout << "The inversion mesh will not be refined." << endl;
  }
  cout << "Will the inversion model be recorded at each iteration?";
  if (record_process == 1) {
    cout << " Yes" << endl;
  } else {
    cout << " No" << endl;
  }
  cout << endl;
  cout << "Use cross-gradient constraint model? "
       << ((use_cross_gradient_constraint == true) ? ("Yes") : ("No")) << endl;
  cout << "Use petrophysical constraint model? "
       << ((use_petrophysical_constraint == true) ? ("Yes") : ("No")) << endl;
  cout << endl;
}