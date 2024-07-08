#ifndef GAUSSNEWTONINVERSION_H
#define GAUSSNEWTONINVERSION_H
#include "AdaptiveInversion.h"
class GaussNewtonInversion : public AdaptiveInversion
{
public:
  GaussNewtonInversion()
      : AdaptiveInversion(),
        record_process(false),
        cg_tol(1e-6),
        stag_tol(0.025),
        GN_iter(10),
        Lp_s(2),
        Lp_x(2),
        Lp_y(2),
        Lp_z(2),
        Lp_t(2),
        epsilon2(1e-9),
        method_id(1) {}
  GaussNewtonInversion(const Mesh &mesh_,
                       const Observation &ob_,
                       unsigned long long field_flag_)
      : AdaptiveInversion(mesh_, ob_, field_flag_),
        record_process(false),
        cg_tol(1e-6),
        stag_tol(0.025),
        GN_iter(10),
        Lp_s(2),
        Lp_x(2),
        Lp_y(2),
        Lp_z(2),
        Lp_t(2),
        epsilon2(1e-9),
        method_id(1) {}

  void invert() override;
  void invert_with_Eigen_CG();
  void invert_with_own_CG();

  void record_every_iteration() { this->record_process = true; }
  void set_method_id(int method_id0) { this->method_id = method_id0; }

  // when cg_iteration_factor is greater than 1, it's the iteration
  // number of conjugate gradient method; when it's lower than 1, CG
  // iteration number is  cg_iteration_factor times the cells number
  void set_CG_parameter(double tol, double cg_iteration_factor);
  void set_stagnation_tolerance(double x) { stag_tol = x; }
  void set_max_GN_iterations(int x) { GN_iter = x; }
  void set_Lp_inversion_parameter(int Lp_, double epsilon2_)
  {

    this->Lp_s = Lp_;
    this->Lp_x = Lp_;
    this->Lp_y = Lp_;
    this->Lp_z = Lp_;
    this->Lp_t = Lp_;
    this->epsilon2 = epsilon2_;
    cout << "Lp_s: " << Lp_s << ", Lp_x: " << Lp_x << ", Lp_y: " << Lp_y << ", Lp_z: " << Lp_z << ", Lp_t: " << Lp_t << ", epsilon**2: " << epsilon2 << endl;
  }
  void set_Lp_inversion_parameter(int Lp_s_, int Lp_x_, int Lp_y_, int Lp_z_, int Lp_t_, double epsilon2_)
  {
    this->Lp_s = Lp_s_;
    this->Lp_x = Lp_x_;
    this->Lp_y = Lp_y_;
    this->Lp_z = Lp_z_;
    this->Lp_t = Lp_t_;
    this->epsilon2 = epsilon2_;
    cout << "Lp_s: " << Lp_s << ", Lp_x: " << Lp_x << ", Lp_y: " << Lp_y << ", Lp_z: " << Lp_z << ", Lp_t" << Lp_t << ", epsilon**2: " << epsilon2 << endl;
  }
  VectorXd solve_cg(const double tol,
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
                    const SMatrix &WL_tz);

  void display_inversion_parameters() const;
  // void invert_Gauss_Newton_ada(int maxit1, int maxit2, int maxit3, VectorXd m_ini,
  //                              double tol1, double tol2, double stagnate_tol, double lambda0, int ada_level = 5, double percentage = 0.1, double damp = 0.7, double n = 1);

protected:
  int Lp_s; // Lp norm for distance from reference model
  int Lp_x; // Lp norm for smoothness constraint
  int Lp_y; // Lp norm for smoothness constraint
  int Lp_z; // Lp norm for smoothness constraint
  int Lp_t; // Lp norm for cross-gradient constraint
  double epsilon2;
  bool record_process;
  double cg_tol;              // tolerance for conjugate gradient
  double cg_iteration_factor; // CG iteration number
  double stag_tol;            // stagnation
  int GN_iter;                // Gauss-Newton iteration number
  int method_id;
};

#endif