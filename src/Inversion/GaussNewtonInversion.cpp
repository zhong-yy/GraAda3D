#include "GaussNewtonInversion.h"
void GaussNewtonInversion::invert()
{
  // this->set_use_wavelet(false);
  // this->compute_G();
  // this->invert_with_Eigen_CG();
  if (this->use_wavelet)
  {
    method_id = 2;
  }
  if (method_id == 0)
  {
    this->set_use_wavelet(false);
    this->compute_G();
    this->invert_with_Eigen_CG(); // not recommended , because the code will get stuck when there are a large number of cells
  }
  else if (method_id == 1)
  {
    this->set_use_wavelet(false);
    this->compute_G();
    this->invert_with_own_CG();
  }
  else if (method_id == 2)
  {
    this->set_use_wavelet(true);
    this->compute_G_wavelet();
    this->invert_with_own_CG();
  }
}

void GaussNewtonInversion::invert_with_Eigen_CG()
{
  // set reference model using interpolation
  if (use_petrophysical_constraint == true)
  {
    for (int i = 0; i < Nm; i++)
    {
      Cell *c = this->mesh.leaf_cells[i];
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

  if (use_cross_gradient_constraint == true)
  {
    std::cout << "Cross gradient constraint is used" << endl;
    for (int i = 0; i < Nm; i++)
    {
      Cell *c = this->mesh.leaf_cells[i];
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
  SMatrix WL_s(Nm, Nm);
  SMatrix WL_x(Nm, Nm);
  SMatrix WL_y(Nm, Nm);
  SMatrix WL_z(Nm, Nm);
  SMatrix WL_tx(Nm, Nm);
  SMatrix WL_ty(Nm, Nm);
  SMatrix WL_tz(Nm, Nm);
  WL_s.setIdentity();
  WL_x.setIdentity();
  WL_y.setIdentity();
  WL_z.setIdentity();
  WL_tx.setIdentity();
  WL_ty.setIdentity();
  WL_tz.setIdentity();

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
  VectorXd Txm;
  VectorXd Tym;
  VectorXd Tzm;

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
  SMatrix P(Nm, Nm);
  P.setZero();
  P.reserve(Nm);
  P.makeCompressed();

  double lambda = this->max_lambda;

  VectorXd delta_x(Nm), m_trial(Nm); // m_trial_last(Nm), m_trial_last2(Nm);
  double misfit, misfit_last1, misfit_last2, misfit_last_iteration;

  if (use_petrophysical_constraint)
  {
    m = m0;
    for (int i = 0; i < Nm; i++)
    {
      if (m(i) > m_max(i) || abs(m(i) - m_max(i)) < 1e-7)
      {
        m(i) = m_max(i) - 0.01;
      }
      if (m(i) < m_min(i) || abs(m(i) - m_min(i)) < 1e-7)
      {
        m(i) = m_min(i) + 0.01;
      }
    }
  }
  else
  {
    m = m_ini;
    for (int i = 0; i < Nm; i++)
    {
      if (m(i) > m_max(i) || abs(m(i) - m_max(i)) < 1e-7)
      {
        m(i) = 0.5 * (m_max(i) + m_min(i));
      }
      if (m(i) < m_min(i) || abs(m(i) - m_min(i)) < 1e-7)
      {
        m(i) = 0.5 * (m_max(i) + m_min(i));
      }
    }
  }

  misfit = sqrt((Wd * (G * m - dobs)).squaredNorm() / Nd);

  misfit_last_iteration = 5 * misfit;

  for (int id = 0; id < Nm; id++)
  {
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
  std::cout << "Set G to A" << endl;
  for (size_t r_id = 0; r_id < Nd; r_id++)
  {
    for (size_t c_id = 0; c_id < Nm; c_id++)
    {
      // cout << r_id << ", " << c_id << endl;
      A.insert(r_id, c_id) = Wd.coeffRef(r_id, r_id) * G(r_id, c_id);
    }
  }
  b.head(Nd) = Wd * dobs;

  for (i = 0; i < GN_iter; i++)
  {
    std::cout << "The " << i + 1 << " th"
              << " Gauss-Newton Iteration:" << endl;
    std::cout << "Number of elements: " << Nm << endl;

    double min_cell_dx, min_cell_dy, min_cell_dz;
    int max_lev;

    mesh.get_minimum_size(min_cell_dx, min_cell_dy, min_cell_dz, max_lev);
    std::cout << "The smallest cell: (dx= " << min_cell_dx << " m, "
              << "dy= " << min_cell_dy << " m, "
              << "dz= " << min_cell_dz << " m)" << endl;

    m_trial = m;
    misfit_last1 = 10 * misfit;  // misfit for last lambda
    misfit_last2 = 100 * misfit; // misfit for the lambda before last lambda
    misfit_last_iteration = misfit;
    lambda = this->max_lambda;

    A.middleRows(Nd, Nm) = sqrt(lambda) * WL_s * Ws;
    A.middleRows(Nd + Nm, Nm) = sqrt(lambda) * WL_z * Wz;
    A.middleRows(Nd + 2 * Nm, Nm) = sqrt(lambda) * WL_x * Wx;
    A.middleRows(Nd + 3 * Nm, Nm) = sqrt(lambda) * WL_y * Wy;
    // cout << S_crg.rows() << endl;
    // cout << S_crg.cols() << endl;
    if (use_cross_gradient_constraint == true)
    {
      A.middleRows(Nd + 4 * Nm, Nm) = this->a_crg * WL_tz * S_crg * T_z;
      A.middleRows(Nd + 5 * Nm, Nm) = this->a_crg * WL_tx * S_crg * T_x;
      A.middleRows(Nd + 6 * Nm, Nm) = this->a_crg * WL_ty * S_crg * T_y;
    }
    A.makeCompressed();

    b.segment(Nd, Nm) = sqrt(lambda) * WL_s * Ws_m0;
    // b.segment(Nd + Nm, Nm) = sqrt(lambda) * Wr_m0;
    // b.segment(Nd + 2 * Nm, Nm) = sqrt(lambda) * Wtheta_m0;
    // b.segment(Nd + 3 * Nm, Nm) = sqrt(lambda) * Wphi_m0;

    if (cg_iteration_factor < 1)
    {
      lscg.setMaxIterations(cg_iteration_factor * Nm);
    }
    else
    {
      lscg.setMaxIterations(int(cg_iteration_factor));
    }
    for (int id = 0; id < Nm; id++)
    {
      int n = 1.0;
      P.coeffRef(id, id) = n * (m_max(id) - m(id)) * (m(id) - m_min(id)) /
                           (m_max(id) - m_min(id));
    }
    delta_x.setZero();
    for (j = 0; j < n_lambda; j++)
    {
      lscg.compute(A * P);
      // delta_x = lscg.solveWithGuess(b - A * m, delta_x);
      delta_x = lscg.solve(b - A * m);

      // m_trial_last2 = m_trial_last;
      // m_trial_last = m_trial;
      double temp, temp2, temp3;
      // update
      for (int id = 0; id < Nm; id++)
      {
        temp = exp(delta_x(id));
        temp2 = (m_min(id) * (m_max(id) - m(id)) +
                 m_max(id) * (m(id) - m_min(id)) * temp);
        temp3 = ((m_max(id) - m(id)) + (m(id) - m_min(id)) * temp);
        if (std::isinf(temp) || std::isinf(temp2))
        {
          m_trial(id) = (m_max(id) - 1e-6 + m_trial(id)) / 2.0;
          if (id == 0)
          {
            m_trial(id) = 0.5 * (m_trial(id) + m_trial(id + 1));
          }
          else
          {
            m_trial(id) = 0.5 * (m_trial(id) + m_trial(id - 1));
          }
        }
        else if (fabs(temp3) < 1e-10)
        {
          m_trial(id) = (m_min(id) + 1e-6 + m_trial(id)) / 2.0;
          if (id == 0)
          {
            m_trial(id) = 0.5 * (m_trial(id) + m_trial(id + 1));
          }
          else
          {
            m_trial(id) = 0.5 * (m_trial(id) + m_trial(id - 1));
          }
        }
        else
        {
          m_trial(id) = temp2 / temp3;
        }
      }

      misfit_last2 = misfit_last1;
      misfit_last1 = misfit;
      misfit = sqrt((Wd * (G * m_trial - dobs)).squaredNorm() / Nd);
      // cout << "  Lambda=" << lambda << ", ";
      std::cout << "  #" << j << " lambda=" << lambda << ", misfit=" << misfit
                << " (CG iterations: " << lscg.iterations() << ","
                << " error: " << lscg.error() << ")" << std::endl;
      out_lambda_misfit << scientific << setw(25) << setprecision(14) << lambda
                        << scientific << setw(25) << setprecision(14)
                        << (misfit) << endl;

      if (j == 0 || misfit < misfit_opt)
      {
        misfit_opt = misfit;
        m_opt = m_trial;
        lambda_opt = lambda;
      }
      if ((std::abs(misfit - misfit_last1) / min(misfit, misfit_last1) <
           stag_tol) &&
          (std::abs(misfit_last1 - misfit_last2) /
               min(misfit_last1, misfit_last2) <
           stag_tol) &&
          j > 1)
      {
        std::cout << "  Misfit stagnates, go to next Gauss-Newton iteration."
                  << endl;
        break;
      }
      if (misfit > misfit_last1 && misfit_last1 > misfit_last2 && j > 1)
      {
        std::cout
            << "  Misfit starts to increase, go to next Gauss-Newton iteration."
            << endl;
        break;
      }
      else
      {
        lambda = lambda * lambda_decreasing_rate;

        Wsm_m0 = Ws * (m_opt - m0);
        Wxm = Wx * m_opt;
        Wym = Wy * m_opt;
        Wzm = Wz * m_opt;

        for (int i = 0; i < Nm; i++)
        {

          WL_s.coeffRef(i, i) =
              pow(Wsm_m0(i) * Wsm_m0(i) + epsilon2, -(2 - Lp_s) / 4.0);
          WL_x.coeffRef(i, i) =
              pow(Wxm(i) * Wxm(i) + epsilon2, -(2 - Lp_x) / 4.0);
          WL_y.coeffRef(i, i) =
              pow(Wym(i) * Wym(i) + epsilon2, -(2 - Lp_y) / 4.0);
          WL_z.coeffRef(i, i) =
              pow(Wzm(i) * Wzm(i) + epsilon2, -(2 - Lp_z) / 4.0);
        }
        if ((Lp_t != 2) && (use_cross_gradient_constraint == true))
        {
          Txm = T_x * m_opt;
          Tym = T_y * m_opt;
          Tzm = T_z * m_opt;
          for (int i = 0; i < Nm; i++)
          {
            WL_tx.coeffRef(i, i) = pow(Txm(i) * Txm(i) + epsilon2, -(2 - Lp_t) / 4.0);
            WL_ty.coeffRef(i, i) = pow(Tym(i) * Tym(i) + epsilon2, -(2 - Lp_t) / 4.0);
            WL_tz.coeffRef(i, i) = pow(Tzm(i) * Tzm(i) + epsilon2, -(2 - Lp_t) / 4.0);
          }
        }

        A.middleRows(Nd, Nm) = sqrt(lambda) * WL_s * Ws;
        A.middleRows(Nd + Nm, Nm) = sqrt(lambda) * WL_z * Wz;
        A.middleRows(Nd + 2 * Nm, Nm) = sqrt(lambda) * WL_x * Wx;
        A.middleRows(Nd + 3 * Nm, Nm) = sqrt(lambda) * WL_y * Wy;
        if (use_cross_gradient_constraint == true)
        {
          A.middleRows(Nd + 4 * Nm, Nm) = this->a_crg * S_crg * WL_tz * T_z;
          A.middleRows(Nd + 5 * Nm, Nm) = this->a_crg * S_crg * WL_tx * T_x;
          A.middleRows(Nd + 6 * Nm, Nm) = this->a_crg * S_crg * WL_ty * T_y;
        }

        b.segment(Nd, Nm) = sqrt(lambda) * WL_s * Ws_m0;
      }
      if (misfit < target_misfit)
      {
        break;
      }
    }
    // lambda=lambda_opt;
    m = m_opt;
    misfit = misfit_opt;

    out_iteration_misfit << fixed << setw(25) << i + 1 << fixed << setw(25)
                         << ((j == n_lambda) ? n_lambda : (j + 1)) << scientific
                         << setw(25) << setprecision(14) << misfit << endl;

    if (record_process)
    {
      this->set_density_to_mesh();
      this->set_reference_model_to_mesh();

#ifdef USE_NETCDF
      if (use_cross_gradient_constraint == true)
      {
        this->mesh.out_model_netcdf(
            string("Structural_constraint_at_") + to_string(i) + string(".nc"),
            2, "crg", "");
      }
      if (use_petrophysical_constraint == true)
      {
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

    if (misfit < target_misfit || i == GN_iter - 1)
    {
      if (misfit < target_misfit)
      {
        std::cout << "Stop iteration, because the target misift has been achieved"
                  << endl;
      }
      else if (i == GN_iter - 1)
      {
        std::cout << "Stop iteration, because the maximum iteration has been reached"
                  << endl;
      }
      std::cout << "Gauss-Newton iteration number:" << i + 1 << endl;
      break;
    }

    bool flag1 = abs(misfit - misfit_last_iteration) /
                     min(misfit, misfit_last_iteration) <
                 stag_tol;
    if (c_refinement >= max_refinement_number && flag1)
    {
      std::cout << "Stop iteration, because the misfit stagnated" << endl;
      std::cout << "Gauss-Newton iteration number:" << i + 1 << endl;
      break;
    }

    if ((c_refinement < max_refinement_number) &&
        (++i_between_refinement % interval_between_refinements == 0 || flag1))
    {
      std::cout << "\n";
      i_between_refinement = 0;
      ++c_refinement;
      if (flag1)
      {
        std::cout << "Misfit stagnates, refine mesh ..." << endl;
      }
      else
      {
        std::cout << interval_between_refinements
                  << ((i_between_refinement == 1) ? string(" iteration has")
                                                  : string(" iterations have"))
                  << " passed since last refinement, refine mesh ..." << endl;
      }
      std::cout << c_refinement
                << (c_refinement == 1
                        ? ("st")
                        : (c_refinement == 2
                               ? ("nd")
                               : (c_refinement == 3 ? ("rd") : ("th"))))
                << " refinement." << endl
                << endl;

      if (use_cross_gradient_constraint && use_petrophysical_constraint)
      {
        refine_mesh(refinement_percentage, *interpolator_m0s, *interpolator_m0);
        // std::cout<<"stop here"<<endl;
        // abort();
      }
      else if (use_cross_gradient_constraint &&
               (!use_petrophysical_constraint))
      {
        refine_mesh(refinement_percentage, *interpolator_m0s, "crg");
        // std::cout<<"2"<<endl;
      }
      else if (use_petrophysical_constraint &&
               (!use_cross_gradient_constraint))
      {
        refine_mesh(refinement_percentage, *interpolator_m0, "pet");
        // std::cout << "2" << endl;
      }
      else
      {
        refine_mesh(refinement_percentage);
        // std::cout<<"2"<<endl;
      }

      std::cout << "Finished refinement" << endl;

      std::cout << "Initialize Ws" << endl;
      Ws.resize(Nm, Nm);
      Ws.reserve(2 * Nm);
      Ws = a_s * S_s * V * D_s * Z;
      Ws.makeCompressed();
      Ws.data().squeeze();

      std::cout << "Initialize Wx" << endl;
      Wx.resize(Nm, Nm);
      Wx.reserve(2 * Nm);
      Wx = a_x * S_x * V * D_x1 * Z;
      Wx.makeCompressed();
      Wx.data().squeeze();

      std::cout << "Initialize Wy" << endl;
      Wy.resize(Nm, Nm);
      Wy.reserve(2 * Nm);
      Wy = a_y * S_y * V * D_y1 * Z;
      Wy.makeCompressed();
      Wy.data().squeeze();

      std::cout << "Initialize Wz" << endl;
      Wz.resize(Nm, Nm);
      Wz.reserve(2 * Nm);
      Wz = a_z * S_z * V * D_z1 * Z;
      Wz.makeCompressed();
      Wz.data().squeeze();

      // std::cout<<"1"<<endl;
      // WL_d.reserve(Nm);
      WL_s.resize(Nm, Nm);
      WL_s.reserve(Nm);
      WL_x.resize(Nm, Nm);
      WL_x.reserve(Nm);
      WL_y.resize(Nm, Nm);
      WL_y.reserve(Nm);
      WL_z.resize(Nm, Nm);
      WL_z.reserve(Nm);
      WL_tx.resize(Nm, Nm);
      WL_tx.reserve(Nm);
      WL_ty.resize(Nm, Nm);
      WL_ty.reserve(Nm);
      WL_tz.resize(Nm, Nm);
      WL_tz.reserve(Nm);

      WL_s.setIdentity();
      WL_x.setIdentity();
      WL_y.setIdentity();
      WL_z.setIdentity();
      WL_tx.setIdentity();
      WL_ty.setIdentity();
      WL_tz.setIdentity();

      Wsm_m0.resize(Nm);
      Wxm.resize(Nm);
      Wym.resize(Nm);
      Wzm.resize(Nm);
      Txm.resize(Nm);
      Tym.resize(Nm);
      Tzm.resize(Nm);
      Wsm_m0 = Ws * (m - m0);
      Wxm = Wx * m;
      Wym = Wy * m;
      Wzm = Wz * m;

      for (int i = 0; i < Nm; i++)
      {
        WL_s.coeffRef(i, i) =
            pow(Wsm_m0(i) * Wsm_m0(i) + epsilon2, -(2 - Lp_s) / 4.0);
        WL_x.coeffRef(i, i) = pow(Wxm(i) * Wxm(i) + epsilon2, -(2 - Lp_x) / 4.0);
        WL_y.coeffRef(i, i) = pow(Wym(i) * Wym(i) + epsilon2, -(2 - Lp_y) / 4.0);
        WL_z.coeffRef(i, i) = pow(Wzm(i) * Wzm(i) + epsilon2, -(2 - Lp_z) / 4.0);
      }

      // update_S_crg();
      // std::cout<<S_crg<<endl;

      Ws_m0 = Ws * m0;
      if (use_cross_gradient_constraint == true)
      {
        this->set_Tmatrix();
        if (Lp_t != 2)
        {
          Txm = T_x * m;
          Tym = T_y * m;
          Tzm = T_z * m;
          for (int i = 0; i < Nm; i++)
          {
            WL_tx.coeffRef(i, i) = pow(Txm(i) * Txm(i) + epsilon2, -(2 - Lp_t) / 4.0);
            WL_ty.coeffRef(i, i) = pow(Tym(i) * Tym(i) + epsilon2, -(2 - Lp_t) / 4.0);
            WL_tz.coeffRef(i, i) = pow(Tzm(i) * Tzm(i) + epsilon2, -(2 - Lp_t) / 4.0);
          }
        }
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
      std::cout << "Set G to A" << endl;
      for (size_t r_id = 0; r_id < Nd; r_id++)
      {
        for (size_t c_id = 0; c_id < Nm; c_id++)
        {
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
    // std::cout << i << endl;
  }
  if (i == GN_iter)
  {
    std::cout << "Gauss-Newton iteration number:" << GN_iter << endl;
  }

  this->set_density_to_mesh();

  final_lambda = lambda;
  final_misfit = misfit;
  std::cout << "Lambda=" << lambda << ", ";
  std::cout << "Misfit=" << misfit << endl;
}

void GaussNewtonInversion::invert_with_own_CG()
{
  // set reference model using interpolation
  if (use_petrophysical_constraint == true)
  {
    for (int i = 0; i < Nm; i++)
    {
      Cell *c = this->mesh.leaf_cells[i];
      double xc, yc, zc;
      c->get_center(xc, yc, zc);
      array<double, 3> args = {xc, yc, zc};
      double val = (*interpolator_m0).interp(args.begin());
      c->set_parameter(val, 1);
    }
    // std::cout << "m0" << endl;
    this->m0.resize(Nm);
    mesh.get_model_parameter_from_mesh(m0, 1);
  }
  // cout << "1" << endl;

  if (use_cross_gradient_constraint == true)
  {
    std::cout << "Cross gradient constraint is used" << endl;
    for (int i = 0; i < Nm; i++)
    {
      Cell *c = this->mesh.leaf_cells[i];
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
  // cout << "x" << endl;

  SMatrix WL_s(Nm, Nm);
  SMatrix WL_x(Nm, Nm);
  SMatrix WL_y(Nm, Nm);
  SMatrix WL_z(Nm, Nm);
  SMatrix WL_tx(Nm, Nm);
  SMatrix WL_ty(Nm, Nm);
  SMatrix WL_tz(Nm, Nm);
  WL_s.setIdentity();
  WL_x.setIdentity();
  WL_y.setIdentity();
  WL_z.setIdentity();
  WL_tx.setIdentity();
  WL_ty.setIdentity();
  WL_tz.setIdentity();

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
  VectorXd Txm;
  VectorXd Tym;
  VectorXd Tzm;

  VectorXd Ws_m0 = Ws * m0;

  // int nrow =
  //     (use_cross_gradient_constraint == true) ? (Nd + 7 * Nm) : (Nd + 4 * Nm);
  // SMatrix A(nrow, Nm);
  // A.reserve(Nd * Nm + Ws.nonZeros() + Wx.nonZeros() + Wy.nonZeros() +
  //           Wz.nonZeros() + 6 * Nm);
  // VectorXd b(nrow);
  // A.setZero();
  // b.setZero();
  // // A.topRows(Nd) = (Wd * G).sparseView();
  // b.head(Nd) = Wd * dobs;
  // std::cout<<"x"<<Nm<<endl;

  SMatrix P(Nm, Nm);
  P.setZero();
  P.reserve(Nm);
  P.makeCompressed();

  double lambda = this->max_lambda;

  VectorXd delta_x(Nm), m_trial(Nm); // m_trial_last(Nm), m_trial_last2(Nm);
  double misfit, misfit_last1, misfit_last2, misfit_last_iteration;

  if (use_petrophysical_constraint)
  {
    m = m0;
    for (int i = 0; i < Nm; i++)
    {
      if (m(i) > m_max(i) || abs(m(i) - m_max(i)) < 1e-7)
      {
        m(i) = m_max(i) - 0.01;
      }
      if (m(i) < m_min(i) || abs(m(i) - m_min(i)) < 1e-7)
      {
        m(i) = m_min(i) + 0.01;
      }
    }
  }
  else
  {
    m = m_ini;
    for (int i = 0; i < Nm; i++)
    {
      if (m(i) > m_max(i) || abs(m(i) - m_max(i)) < 1e-7)
      {
        m(i) = 0.5 * (m_max(i) + m_min(i));
      }
      if (m(i) < m_min(i) || abs(m(i) - m_min(i)) < 1e-7)
      {
        m(i) = 0.5 * (m_max(i) + m_min(i));
      }
    }
  }
  VectorXd dpre(Nd);
  this->G_vec_mul(m, dpre);
  misfit = sqrt((Wd * (dpre - dobs)).squaredNorm() / Nd);

  misfit_last_iteration = 5 * misfit;

  for (int id = 0; id < Nm; id++)
  {
    P.coeffRef(id, id) =
        (m_max(id) - m(id)) * (m(id) - m_min(id)) / (m_max(id) - m_min(id));
  }
  // LeastSquaresConjugateGradient<SMatrix> lscg;
  // lscg.setTolerance(cg_tol);

  ofstream out_lambda_misfit("lambda_misfit_GN");
  out_lambda_misfit << setw(25) << "lambda" << setw(25) << "misfit" << endl;

  ofstream out_iteration_misfit("Iteration_misfit_GN");
  out_iteration_misfit << setw(25) << "Iteranation number" << setw(25)
                       << "number of tried regulartion parameters" << setw(25)
                       << "misfit" << endl;
  // cout << "B" << endl;

  int i, j;
  double lambda_opt;
  double misfit_opt;
  VectorXd m_opt;
  misfit_opt = misfit;
  m_opt = m;

  int c_refinement = 0;
  int i_between_refinement = 0;

  for (i = 0; i < GN_iter; i++)
  {
    std::cout << "The " << i + 1 << " th"
              << " Gauss-Newton Iteration:" << endl;
    std::cout << "Number of elements: " << Nm << endl;

    double min_cell_dx, min_cell_dy, min_cell_dz;
    int max_lev;

    mesh.get_minimum_size(min_cell_dx, min_cell_dy, min_cell_dz, max_lev);
    std::cout << "The smallest cell: (dx= " << min_cell_dx << " m, "
              << "dy= " << min_cell_dy << " m, "
              << "dz= " << min_cell_dz << " m)" << endl;

    m_trial = m;
    misfit_last1 = 10 * misfit;  // misfit for last lambda
    misfit_last2 = 100 * misfit; // misfit for the lambda before last lambda
    misfit_last_iteration = misfit;
    lambda = this->max_lambda;

    int cg_maxit = 0;
    if (cg_iteration_factor < 1)
    {
      cg_maxit = cg_iteration_factor * Nm;
    }
    else
    {
      cg_maxit = int(cg_iteration_factor);
    }
    for (int id = 0; id < Nm; id++)
    {
      int n = 1.0;
      P.coeffRef(id, id) = n * (m_max(id) - m(id)) * (m(id) - m_min(id)) /
                           (m_max(id) - m_min(id));
    }
    delta_x.setZero();
    for (j = 0; j < n_lambda; j++)
    {
      // lscg.compute(A * P);
      // delta_x = lscg.solve(b - A * m);
      double cg_error;
      int cg_iterations;
      delta_x = solve_cg(cg_tol, cg_maxit, cg_error, cg_iterations,
                         lambda, m, m0, Wd, dobs, P, Ws, Wx, Wy, Wz,
                         T_x, T_y, T_z, WL_s, WL_x, WL_y, WL_z, WL_tx, WL_ty, WL_tz);

      // m_trial_last2 = m_trial_last;
      // m_trial_last = m_trial;
      double temp, temp2, temp3;
      // update
      for (int id = 0; id < Nm; id++)
      {
        temp = exp(delta_x(id));
        temp2 = (m_min(id) * (m_max(id) - m(id)) +
                 m_max(id) * (m(id) - m_min(id)) * temp);
        temp3 = ((m_max(id) - m(id)) + (m(id) - m_min(id)) * temp);
        if (std::isinf(temp) || std::isinf(temp2))
        {
          m_trial(id) = (m_max(id) - 1e-6 + m_trial(id)) / 2.0;
          if (id == 0)
          {
            m_trial(id) = 0.5 * (m_trial(id) + m_trial(id + 1));
          }
          else
          {
            m_trial(id) = 0.5 * (m_trial(id) + m_trial(id - 1));
          }
        }
        else if (fabs(temp3) < 1e-10)
        {
          m_trial(id) = (m_min(id) + 1e-6 + m_trial(id)) / 2.0;
          if (id == 0)
          {
            m_trial(id) = 0.5 * (m_trial(id) + m_trial(id + 1));
          }
          else
          {
            m_trial(id) = 0.5 * (m_trial(id) + m_trial(id - 1));
          }
        }
        else
        {
          m_trial(id) = temp2 / temp3;
        }
      }

      misfit_last2 = misfit_last1;
      misfit_last1 = misfit;
      this->G_vec_mul(m_trial, dpre);
      misfit = sqrt((Wd * (dpre - dobs)).squaredNorm() / Nd);
      // std::cout << "  Lambda=" << lambda << ", ";
      std::cout << "  #" << j << " lambda=" << lambda << ", misfit=" << misfit
                << " (CG iterations: " << cg_iterations << ","
                << " error: " << cg_error << ")" << std::endl;
      out_lambda_misfit << scientific << setw(25) << setprecision(14) << lambda
                        << scientific << setw(25) << setprecision(14)
                        << (misfit) << endl;

      if (j == 0 || misfit < misfit_opt)
      {
        misfit_opt = misfit;
        m_opt = m_trial;
        lambda_opt = lambda;
      }
      if ((std::abs(misfit - misfit_last1) / min(misfit, misfit_last1) <
           stag_tol) &&
          (std::abs(misfit_last1 - misfit_last2) /
               min(misfit_last1, misfit_last2) <
           stag_tol) &&
          j > 1)
      {
        std::cout << "  Misfit stagnates, go to next Gauss-Newton iteration."
                  << endl;
        break;
      }
      if (misfit > misfit_last1 && misfit_last1 > misfit_last2 && j > 1)
      {
        std::cout
            << "  Misfit starts to increase, go to next Gauss-Newton iteration."
            << endl;
        break;
      }
      else
      {
        // cout << "Reduce lambda" << endl;
        lambda = lambda * lambda_decreasing_rate;

        if (Lp_s != 2)
        {
          Wsm_m0 = Ws * (m_opt - m0);
          for (int ir = 0; ir < Nm; ir++)
          {
            WL_s.coeffRef(ir, ir) =
                pow(Wsm_m0(ir) * Wsm_m0(ir) + epsilon2, -(2 - Lp_s) / 4.0);
          }
        }
        if (Lp_x != 2)
        {
          Wxm = Wx * m_opt;
          for (int ir = 0; ir < Nm; ir++)
          {
            WL_x.coeffRef(ir, ir) =
                pow(Wxm(ir) * Wxm(ir) + epsilon2, -(2 - Lp_x) / 4.0);
          }
        }
        if (Lp_y != 2)
        {
          Wym = Wy * m_opt;
          for (int ir = 0; ir < Nm; ir++)
          {
            WL_y.coeffRef(ir, ir) =
                pow(Wym(ir) * Wym(ir) + epsilon2, -(2 - Lp_y) / 4.0);
          }
        }
        if (Lp_z != 2)
        {
          Wzm = Wz * m_opt;
          for (int ir = 0; ir < Nm; ir++)
          {
            WL_z.coeffRef(ir, ir) =
                pow(Wzm(ir) * Wzm(ir) + epsilon2, -(2 - Lp_z) / 4.0);
          }
        }
        // for (int i = 0; i < Nm; i++)
        // {
        //   WL_s.coeffRef(i, i) =
        //       pow(Wsm_m0(i) * Wsm_m0(i) + epsilon2, -(2 - Lp_s) / 4.0);
        //   WL_x.coeffRef(i, i) =
        //       pow(Wxm(i) * Wxm(i) + epsilon2, -(2 - Lp_x) / 4.0);
        //   WL_y.coeffRef(i, i) =
        //       pow(Wym(i) * Wym(i) + epsilon2, -(2 - Lp_y) / 4.0);
        //   WL_z.coeffRef(i, i) =
        //       pow(Wzm(i) * Wzm(i) + epsilon2, -(2 - Lp_z) / 4.0);
        // }

        if ((Lp_t != 2) && (use_cross_gradient_constraint == true))
        {

          Txm = T_x * m_opt;
          Tym = T_y * m_opt;
          Tzm = T_z * m_opt;
          for (int i = 0; i < Nm; i++)
          {
            WL_tx.coeffRef(i, i) = pow(Txm(i) * Txm(i) + epsilon2, -(2 - Lp_t) / 4.0);
            WL_ty.coeffRef(i, i) = pow(Tym(i) * Tym(i) + epsilon2, -(2 - Lp_t) / 4.0);
            WL_tz.coeffRef(i, i) = pow(Tzm(i) * Tzm(i) + epsilon2, -(2 - Lp_t) / 4.0);
          }
        }

        // A.middleRows(Nd, Nm) = sqrt(lambda) * WL_s * Ws;
        // A.middleRows(Nd + Nm, Nm) = sqrt(lambda) * WL_z * Wz;
        // A.middleRows(Nd + 2 * Nm, Nm) = sqrt(lambda) * WL_x * Wx;
        // A.middleRows(Nd + 3 * Nm, Nm) = sqrt(lambda) * WL_y * Wy;
        // if (use_cross_gradient_constraint == true)
        // {
        //   A.middleRows(Nd + 4 * Nm, Nm) = this->a_crg * S_crg * T_z;
        //   A.middleRows(Nd + 5 * Nm, Nm) = this->a_crg * S_crg * T_x;
        //   A.middleRows(Nd + 6 * Nm, Nm) = this->a_crg * S_crg * T_y;
        // }

        // b.segment(Nd, Nm) = sqrt(lambda) * WL_s * Ws_m0;
      }
      if (misfit < target_misfit)
      {
        break;
      }
    }
    // lambda=lambda_opt;
    m = m_opt;
    misfit = misfit_opt;

    out_iteration_misfit << fixed << setw(25) << i + 1 << fixed << setw(25)
                         << ((j == n_lambda) ? n_lambda : (j + 1)) << scientific
                         << setw(25) << setprecision(14) << misfit << endl;

    if (record_process)
    {
      this->set_density_to_mesh();
      this->set_reference_model_to_mesh();

#ifdef USE_NETCDF
      if (use_cross_gradient_constraint == true)
      {
        this->mesh.out_model_netcdf(
            string("Structural_constraint_at_") + to_string(i) + string(".nc"),
            2, "crg", "");
      }
      if (use_petrophysical_constraint == true)
      {
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

    if (misfit < target_misfit || i == GN_iter - 1)
    {
      if (misfit < target_misfit)
      {
        std::cout << "Stop iteration, because the target misift has been achieved"
                  << endl;
      }
      else if (i == GN_iter - 1)
      {
        std::cout << "Stop iteration, because the maximum iteration has been reached"
                  << endl;
      }
      std::cout << "Gauss-Newton iteration number:" << i + 1 << endl;
      break;
    }

    bool flag1 = abs(misfit - misfit_last_iteration) /
                     min(misfit, misfit_last_iteration) <
                 stag_tol;
    if (c_refinement >= max_refinement_number && flag1)
    {
      std::cout << "Stop iteration, because the misfit stagnated" << endl;
      std::cout << "Gauss-Newton iteration number:" << i + 1 << endl;
      break;
    }

    if ((c_refinement < max_refinement_number) &&
        (++i_between_refinement % interval_between_refinements == 0 || flag1))
    {
      std::cout << "\n";
      i_between_refinement = 0;
      ++c_refinement;
      if (flag1)
      {
        std::cout << "Misfit stagnates, refine mesh ..." << endl;
      }
      else
      {
        std::cout << interval_between_refinements
                  << ((i_between_refinement == 1) ? string(" iteration has")
                                                  : string(" iterations have"))
                  << " passed since last refinement, refine mesh ..." << endl;
      }
      std::cout << c_refinement
                << (c_refinement == 1
                        ? ("st")
                        : (c_refinement == 2
                               ? ("nd")
                               : (c_refinement == 3 ? ("rd") : ("th"))))
                << " refinement." << endl
                << endl;

      if (use_cross_gradient_constraint && use_petrophysical_constraint)
      {
        refine_mesh(refinement_percentage, *interpolator_m0s, *interpolator_m0);
        // std::cout<<"stop here"<<endl;
        // abort();
      }
      else if (use_cross_gradient_constraint &&
               (!use_petrophysical_constraint))
      {
        refine_mesh(refinement_percentage, *interpolator_m0s, "crg");
        // std::cout<<"2"<<endl;
      }
      else if (use_petrophysical_constraint &&
               (!use_cross_gradient_constraint))
      {
        refine_mesh(refinement_percentage, *interpolator_m0, "pet");
        // std::cout << "2" << endl;
      }
      else
      {
        refine_mesh(refinement_percentage);
        // std::cout<<"2"<<endl;
      }

      std::cout << "Finished refinement" << endl;

      // std::cout << "Initialize Ws" << endl;
      Ws.resize(Nm, Nm);
      Ws.reserve(2 * Nm);
      Ws = a_s * S_s * V * D_s * Z;
      Ws.makeCompressed();
      Ws.data().squeeze();

      // std::cout << "Initialize Wx" << endl;
      Wx.resize(Nm, Nm);
      Wx.reserve(2 * Nm);
      Wx = a_x * S_x * V * D_x1 * Z;
      Wx.makeCompressed();
      Wx.data().squeeze();

      // std::cout << "Initialize Wy" << endl;
      Wy.resize(Nm, Nm);
      Wy.reserve(2 * Nm);
      Wy = a_y * S_y * V * D_y1 * Z;
      Wy.makeCompressed();
      Wy.data().squeeze();

      // std::cout << "Initialize Wz" << endl;
      Wz.resize(Nm, Nm);
      Wz.reserve(2 * Nm);
      Wz = a_z * S_z * V * D_z1 * Z;
      Wz.makeCompressed();
      Wz.data().squeeze();

      // std::cout<<"1"<<endl;
      WL_s.resize(Nm, Nm);
      WL_s.reserve(Nm);
      WL_x.resize(Nm, Nm);
      WL_x.reserve(Nm);
      WL_y.resize(Nm, Nm);
      WL_y.reserve(Nm);
      WL_z.resize(Nm, Nm);
      WL_z.reserve(Nm);
      WL_tx.resize(Nm, Nm);
      WL_tx.reserve(Nm);
      WL_ty.resize(Nm, Nm);
      WL_ty.reserve(Nm);
      WL_tz.resize(Nm, Nm);
      WL_tz.reserve(Nm);

      WL_s.setIdentity();
      WL_x.setIdentity();
      WL_y.setIdentity();
      WL_z.setIdentity();
      WL_tx.setIdentity();
      WL_ty.setIdentity();
      WL_tz.setIdentity();

      Wsm_m0.resize(Nm);
      Wxm.resize(Nm);
      Wym.resize(Nm);
      Wzm.resize(Nm);
      // Wsm_m0 = Ws * (m - m0);
      // Wxm = Wx * m;
      // Wym = Wy * m;
      // Wzm = Wz * m;
      if (Lp_s != 2)
      {
        Wsm_m0 = Ws * (m - m0);
        for (int ir = 0; ir < Nm; ir++)
        {
          WL_s.coeffRef(ir, ir) =
              pow(Wsm_m0(ir) * Wsm_m0(ir) + epsilon2, -(2 - Lp_s) / 4.0);
        }
      }
      if (Lp_x != 2)
      {
        Wxm = Wx * m;
        for (int ir = 0; ir < Nm; ir++)
        {
          WL_x.coeffRef(ir, ir) =
              pow(Wxm(ir) * Wxm(ir) + epsilon2, -(2 - Lp_x) / 4.0);
        }
      }
      if (Lp_y != 2)
      {
        Wym = Wy * m;
        for (int ir = 0; ir < Nm; ir++)
        {
          WL_y.coeffRef(ir, ir) =
              pow(Wym(ir) * Wym(ir) + epsilon2, -(2 - Lp_y) / 4.0);
        }
      }
      if (Lp_z != 2)
      {
        Wzm = Wz * m;
        for (int ir = 0; ir < Nm; ir++)
        {
          WL_z.coeffRef(ir, ir) =
              pow(Wzm(ir) * Wzm(ir) + epsilon2, -(2 - Lp_z) / 4.0);
        }
      }
      // for (int i = 0; i < Nm; i++)
      // {
      //   WL_s.coeffRef(i, i) =
      //       pow(Wsm_m0(i) * Wsm_m0(i) + epsilon2, -(2 - Lp_s) / 4.0);
      //   WL_x.coeffRef(i, i) = pow(Wxm(i) * Wxm(i) + epsilon2, -(2 - Lp_x) / 4.0);
      //   WL_y.coeffRef(i, i) = pow(Wym(i) * Wym(i) + epsilon2, -(2 - Lp_y) / 4.0);
      //   WL_z.coeffRef(i, i) = pow(Wzm(i) * Wzm(i) + epsilon2, -(2 - Lp_z) / 4.0);
      // }

      // update_S_crg();
      // std::cout<<S_crg<<endl;

      Ws_m0 = Ws * m0;
      if (use_cross_gradient_constraint == true)
      {
        this->set_Tmatrix();
        if (Lp_t != 2)
        {
          Txm = T_x * m;
          Tym = T_y * m;
          Tzm = T_z * m;
          for (int i = 0; i < Nm; i++)
          {
            WL_tx.coeffRef(i, i) = pow(Txm(i) * Txm(i) + epsilon2, -(2 - Lp_t) / 4.0);
            WL_ty.coeffRef(i, i) = pow(Tym(i) * Tym(i) + epsilon2, -(2 - Lp_t) / 4.0);
            WL_tz.coeffRef(i, i) = pow(Tzm(i) * Tzm(i) + epsilon2, -(2 - Lp_t) / 4.0);
          }
        }
      }

      // b.resize(nrow, 1);
      // b.setZero();
      // b.head(Nd) = Wd * dobs;

      P.resize(Nm, Nm);
      P.setZero();
      P.reserve(Nm);
      P.data().squeeze();

      delta_x.resize(Nm, 1);
      m_trial.resize(Nm, 1);
    }
    // std::cout << i << endl;
  }
  if (i == GN_iter)
  {
    std::cout << "Gauss-Newton iteration number:" << GN_iter << endl;
  }

  this->set_density_to_mesh();

  final_lambda = lambda;
  final_misfit = misfit;
  std::cout << "Lambda=" << lambda << ", ";
  std::cout << "Misfit=" << misfit << endl;
}

void GaussNewtonInversion::set_CG_parameter(double cg_tol_,
                                            double cg_iteration_factor_)
{
  cg_tol = cg_tol_;
  cg_iteration_factor_ = abs(cg_iteration_factor_);
  this->cg_iteration_factor = cg_iteration_factor_;
}

void GaussNewtonInversion::display_inversion_parameters() const
{
  display_info_fields();
  std::cout << endl;

  std::cout << "alpha_s=" << a_s << endl;
  std::cout << "alpha_z=" << a_z << endl;
  std::cout << "alpha_x=" << a_x << endl;
  std::cout << "alpha_y=" << a_y << endl;

  cout << "alpha_crg=" << a_crg << endl;
  // cout << endl;

  if (use_wavelet)
  {
    cout << endl
         << "Wavelet transform is used to accelerate the computation. The "
            "relative threshold is "
         << compression_threshold << endl;
  }

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
  if (cg_iteration_factor < 1)
  {
    cout << "Maximum number of iterations of the conjugate gradient method is "
         << cg_iteration_factor * 100 << "% of the element number" << endl;
  }
  else
  {
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

  if (max_refinement_number > 0)
  {
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
  }
  else
  {
    cout << endl;
    cout << "The inversion mesh will not be refined." << endl;
  }
  cout << "Will the inversion model be recorded at each iteration?";
  if (record_process == 1)
  {
    cout << " Yes" << endl;
  }
  else
  {
    cout << " No" << endl;
  }
  cout << endl;
  cout << "Use cross-gradient constraint model? "
       << ((use_cross_gradient_constraint == true) ? ("Yes") : ("No")) << endl;
  cout << "Use petrophysical constraint model? "
       << ((use_petrophysical_constraint == true) ? ("Yes") : ("No")) << endl;
  cout << endl;
}

VectorXd GaussNewtonInversion::solve_cg(const double tol,
                                        const int &maxit,
                                        double &err,
                                        int &iterations,
                                        const double &lambda,
                                        const VectorXd &mk,
                                        const VectorXd &mref,
                                        const SMatrix &Wd,
                                        const VectorXd &dobs,
                                        const SMatrix &P,
                                        const SMatrix &Ws,
                                        const SMatrix &Wx,
                                        const SMatrix &Wy,
                                        const SMatrix &Wz,
                                        const SMatrix &Tx0,
                                        const SMatrix &Ty0,
                                        const SMatrix &Tz0,
                                        const SMatrix &WL_s,
                                        const SMatrix &WL_x,
                                        const SMatrix &WL_y,
                                        const SMatrix &WL_z,
                                        const SMatrix &WL_tx,
                                        const SMatrix &WL_ty,
                                        const SMatrix &WL_tz)
{
  SMatrix Tz;
  SMatrix Tx;
  SMatrix Ty;
  if (use_cross_gradient_constraint)
  {
    Tz = this->a_crg * S_crg * Tz0;
    Tx = this->a_crg * S_crg * Tx0;
    Ty = this->a_crg * S_crg * Ty0;
  }
  int nrow =
      (use_cross_gradient_constraint == true) ? (Nd + 7 * Nm) : (Nd + 4 * Nm);
  VectorXd s(nrow);
  s.setZero();
  VectorXd dx(Nm);
  dx.setZero();
  VectorXd b(nrow);

  VectorXd Gmk;
  this->G_vec_mul(mk, Gmk);
  b.segment(0, Nd) = Wd * (dobs - Gmk);
  b.segment(Nd, Nm) = -sqrt(lambda) * WL_s * Ws * (mk - mref);
  b.segment(Nd + Nm, Nm) = -sqrt(lambda) * WL_z * Wz * mk;
  b.segment(Nd + 2 * Nm, Nm) = -sqrt(lambda) * WL_x * Wx * mk;
  b.segment(Nd + 3 * Nm, Nm) = -sqrt(lambda) * WL_y * Wy * mk;
  if (use_cross_gradient_constraint == true)
  {
    b.segment(Nd + 4 * Nm, Nm) = -WL_tz * Tz * mk;
    b.segment(Nd + 5 * Nm, Nm) = -WL_tx * Tx * mk;
    b.segment(Nd + 6 * Nm, Nm) = -WL_ty * Ty * mk;
  }
  s = -b;
  VectorXd r(Nm);
  VectorXd vec1(Nm);
  VectorXd vec2(Nm);
  VectorXd vec3(Nm);
  VectorXd vec4(Nm);
  VectorXd vec5(Nm);
  VectorXd vec6(Nm);
  VectorXd vec7(Nm);
  VectorXd vec8(Nm);
  r.setZero();

  VectorXd Wds = Wd * s.segment(0, Nd);
  VectorXd GTwds;
  this->GT_vec_mul(Wds, GTwds);
#pragma omp parallel sections
  {
#pragma omp section
    {
      vec1 = P * GTwds;
      vec2 = sqrt(lambda) * P * Ws.transpose() * WL_s.transpose() * s.segment(Nd, Nm);
    }
#pragma omp section
    {
      vec3 = sqrt(lambda) * P * Wz.transpose() * WL_z.transpose() * s.segment(Nd + Nm, Nm);
      vec4 = sqrt(lambda) * P * Wx.transpose() * WL_x.transpose() *
             s.segment(Nd + 2 * Nm, Nm);
    }
#pragma omp section
    {
      vec5 = sqrt(lambda) * P * Wy.transpose() * WL_y.transpose() *
             s.segment(Nd + 3 * Nm, Nm);
      if (use_cross_gradient_constraint == true)
      {
        vec6 = P * Tz.transpose() * WL_tz.transpose() * s.segment(Nd + 4 * Nm, Nm);
        vec7 = P * Tx.transpose() * WL_tx.transpose() * s.segment(Nd + 5 * Nm, Nm);
        vec8 = P * Ty.transpose() * WL_y.transpose() * s.segment(Nd + 6 * Nm, Nm);
      }
    }
  }
#pragma omp parallel for
  for (int k = 0; k < Nm; k++)
  {
    r(k) = vec1(k) + vec2(k) + vec3(k) + vec4(k) + vec5(k);
    if (use_cross_gradient_constraint == true)
    {
      r(k) = r(k) + vec6(k) + vec7(k) + vec8(k);
    }
  }
  double beta = 0.0;
  VectorXd p(Nm);
  p.setZero();
  VectorXd t(nrow);
  t.setZero();
  double alpha = 0.0;
  VectorXd r_last(Nm);
  VectorXd ATb(Nm);
  ATb.setZero();
  // VectorXd Wdb = Wd * s.segment(0, Nd);
  VectorXd Wdb = Wd * b.segment(0, Nd);
  VectorXd GTwdb;
  this->GT_vec_mul(Wdb, GTwdb);
  ATb += P * GTwdb;
  ATb += sqrt(lambda) * P * Ws.transpose() * WL_s.transpose() * b.segment(Nd, Nm);
  ATb += sqrt(lambda) * P * Wz.transpose() * WL_z.transpose() * b.segment(Nd + Nm, Nm);
  ATb += sqrt(lambda) * P * Wx.transpose() * WL_x.transpose() * b.segment(Nd + 2 * Nm, Nm);
  ATb += sqrt(lambda) * P * Wy.transpose() * WL_y.transpose() * b.segment(Nd + 3 * Nm, Nm);
  if (use_cross_gradient_constraint == true)
  {
    ATb += P * Tz.transpose() * WL_tz.transpose() * b.segment(Nd + 4 * Nm, Nm);
    ATb += P * Tx.transpose() * WL_tx.transpose() * b.segment(Nd + 5 * Nm, Nm);
    ATb += P * Ty.transpose() * WL_ty.transpose() * b.segment(Nd + 6 * Nm, Nm);
  }
  int i;
  double r2, t2, rlast2, ATb2;
  VectorXd Pp;
  VectorXd GPp;

  ATb2 = 0;
#pragma omp parallel for reduction(+ : ATb2)
  for (int k = 0; k < Nm; ++k)
  {
    ATb2 += ATb(k) * ATb(k);
  }
  for (i = 0; i < maxit; i++)
  {
    p = -r + beta * p;
    // t=A*p
    Pp = P * p;

#pragma omp parallel sections
    {
#pragma omp section
      {
        this->G_vec_mul(Pp, GPp);
        // t.segment(0, Nd) = Wd * G * P * p;
        t.segment(0, Nd) = Wd * GPp;
        // t.segment(Nd, Nm) = sqrt(lambda) * Ws * P * p;
      }
#pragma omp section
      {
        t.segment(Nd, Nm) = sqrt(lambda) * Ws * Pp;
        t.segment(Nd + Nm, Nm) = sqrt(lambda) * WL_z * Wz * Pp;
      }
#pragma omp section
      {
        t.segment(Nd + 2 * Nm, Nm) = sqrt(lambda) * WL_x * Wx * Pp;
        t.segment(Nd + 3 * Nm, Nm) = sqrt(lambda) * WL_y * Wy * Pp;
        if (use_cross_gradient_constraint == true)
        {
          t.segment(Nd + 4 * Nm, Nm) = WL_z * Tz * Pp;
          t.segment(Nd + 5 * Nm, Nm) = WL_x * Tx * Pp;
          t.segment(Nd + 6 * Nm, Nm) = WL_y * Ty * Pp;
        }
      }
    }
    r2 = 0;
#pragma omp parallel for reduction(+ : r2)
    for (int k = 0; k < Nm; ++k)
    {
      r2 += r(k) * r(k);
    }

    t2 = 0;
#pragma omp parallel for reduction(+ : t2)
    for (int k = 0; k < nrow; ++k)
    {
      t2 += t(k) * t(k);
    }
    // alpha = r.transpose() * r;
    // alpha /= t.transpose() * t;
    alpha = r2 / t2;
#pragma omp parallel for
    for (int k = 0; k < Nd; k++)
    {
      s(k) += alpha * t(k);
    }
#pragma omp parallel for
    for (int k = 0; k < Nm; k++)
    {
      dx(k) += alpha * p(k);
      // s(k) += alpha * t(k);
      s(k + Nd) += alpha * t(k + Nd);
      s(k + Nd + Nm) += alpha * t(k + Nd + Nm);
      s(k + Nd + 2 * Nm) += alpha * t(k + Nd + 2 * Nm);
      s(k + Nd + 3 * Nm) += alpha * t(k + Nd + 3 * Nm);
      // s(k + Nd + 4 * Nm) += alpha * t(k + Nd + 4 * Nm);
    }
    if (use_cross_gradient_constraint == true)
    {
#pragma omp parallel for
      for (int k = 0; k < Nm; k++)
      {
        s(k + Nd + 4 * Nm) += alpha * t(k + Nd + 4 * Nm);
        s(k + Nd + 5 * Nm) += alpha * t(k + Nd + 5 * Nm);
        s(k + Nd + 6 * Nm) += alpha * t(k + Nd + 6 * Nm);
      }
    }

    // dx += alpha * p;
    // s += alpha * t;

    r_last = r;
    r.setZero();
    Wds = Wd * s.segment(0, Nd);
    this->GT_vec_mul(Wds, GTwds);
    // r += P * G.transpose() * Wd * s.segment(0, Nd);
#pragma omp parallel sections
    {
#pragma omp section
      {
        vec1 = P * GTwds;
        vec2 = sqrt(lambda) * P * Ws.transpose() * WL_s.transpose() * s.segment(Nd, Nm);
      }
#pragma omp section
      {
        vec3 =
            sqrt(lambda) * P * Wz.transpose() * WL_z.transpose() * s.segment(Nd + Nm, Nm);
        vec4 = sqrt(lambda) * P * Wx.transpose() * WL_x.transpose() *
               s.segment(Nd + 2 * Nm, Nm);
      }
#pragma omp section
      {
        vec5 = sqrt(lambda) * P * Wy.transpose() * WL_y.transpose() *
               s.segment(Nd + 3 * Nm, Nm);
        if (use_cross_gradient_constraint == true)
        {
          vec6 = P * Tz.transpose() * WL_tz.transpose() * s.segment(Nd + 4 * Nm, Nm);
          vec7 = P * Tx.transpose() * WL_tx.transpose() * s.segment(Nd + 5 * Nm, Nm);
          vec8 = P * Ty.transpose() * WL_ty.transpose() * s.segment(Nd + 6 * Nm, Nm);
        }
      }
    }
#pragma omp parallel for
    for (int k = 0; k < Nm; k++)
    {
      r(k) = vec1(k) + vec2(k) + vec3(k) + vec4(k) + vec5(k);
      if (use_cross_gradient_constraint == true)
      {
        r(k) = r(k) + vec6(k) + vec7(k) + vec8(k);
      }
    }
    rlast2 = r2;
    r2 = 0.0;
#pragma omp parallel for reduction(+ : r2)
    for (int k = 0; k < Nm; ++k)
    {
      r2 += r(k) * r(k);
    }
    // beta = r.transpose() * r;
    // beta /= r_last.transpose() * r_last;
    // err = r.dot(r) / ATb.dot(ATb);
    beta = r2 / rlast2;
    err = r2 / ATb2;

    if (err < tol)
    {
      iterations = i + 1;
      break;
    }
  }
  if (i == maxit)
  {
    iterations = maxit;
  }
  return dx;
}