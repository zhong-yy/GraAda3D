#include "LinearInversion.h"
void LinearInversion::invert() {
  SMatrix Ws = a_s * S_s * V * D_s * Z;
  SMatrix Wz = a_z * S_z * V * D_z1 * Z;
  SMatrix Wx = a_x * S_x * V * D_x1 * Z;
  SMatrix Wy = a_y * S_y * V * D_y1 * Z;

  VectorXd Ws_m0 = Ws * m0;

  double lambda = this->max_lambda;

  SMatrix A(Nd + 4 * Nm, Nm);
  VectorXd b(Nd + 4 * Nm);

  A.setZero();
  b.setZero();
  A.topRows(Nd) = (Wd * G).sparseView();

  b.head(Nd) = Wd * dobs;
  A.middleRows(Nd, Nm) = sqrt(lambda) * Ws;
  A.middleRows(Nd + Nm, Nm) = sqrt(lambda) * Wz;
  A.middleRows(Nd + 2 * Nm, Nm) = sqrt(lambda) * Wx;
  A.middleRows(Nd + 3 * Nm, Nm) = sqrt(lambda) * Wy;
  b.segment(Nd, Nm) = sqrt(lambda) * Ws_m0;

  string output_lambda_vs_misfit("lambda_vs_misfit");
  ofstream out_lambda_misfit(output_lambda_vs_misfit.c_str());
  out_lambda_misfit << setw(25) << "lambda" << setw(25) << "misfit" << endl;

  // LeastSquaresConjugateGradient<SMatrix,Eigen::IdentityPreconditioner> lscg;
  LeastSquaresConjugateGradient<SMatrix> lscg;
  if (cg_iteration_factor < 5) {
    lscg.setMaxIterations(cg_iteration_factor * Nm);
  } else {
    lscg.setMaxIterations(int(cg_iteration_factor));
  }
  lscg.setTolerance(cg_tol);
  double misfit = 0;
  m = m_ini;

  // cout<<n_lambda<<endl;
  for (int i = 0; i < n_lambda; i++) {
    lscg.compute(A);
    m = lscg.solveWithGuess(b, m);

    double last_misfit = misfit;
    misfit = (Wd * (G * m - dobs)).squaredNorm() / Nd;

    cout << "  #" << i << " lambda=" << lambda << ", misfit=" << misfit
         << " (CG iterations: " << lscg.iterations() << ","
         << " error: " << lscg.error() << ")" << std::endl;

    out_lambda_misfit << scientific << setw(25) << setprecision(14) << lambda
                      << scientific << setw(25) << setprecision(14) << (misfit)
                      << endl;
    if (misfit < target_misfit ||
        (std::abs(last_misfit - misfit) / min(misfit, last_misfit) <
         stag_tol)) {
      cout << "Stop, because the target misift has been achieved or the misfit "
              "stagnated"
           << endl;
      break;
    } else {
      lambda = lambda * lambda_decreasing_rate;

      A.middleRows(Nd, Nm) = sqrt(lambda) * Ws;
      A.middleRows(Nd + Nm, Nm) = sqrt(lambda) * Wz;
      A.middleRows(Nd + 2 * Nm, Nm) = sqrt(lambda) * Wx;
      A.middleRows(Nd + 3 * Nm, Nm) = sqrt(lambda) * Wy;
      b.segment(Nd, Nm) = sqrt(lambda) * Ws_m0;
    }
  }
  cout << "Lambda=" << lambda << endl;
  cout << "Misfit=" << misfit << endl;
  final_lambda = lambda;
  final_misfit = misfit;
}