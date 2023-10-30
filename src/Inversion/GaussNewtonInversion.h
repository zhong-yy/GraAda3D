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
        Lp(2),
        epsilon2(1e-9) {}
  GaussNewtonInversion(const Mesh &mesh_,
                       const Observation &ob_,
                       unsigned long long field_flag_)
      : AdaptiveInversion(mesh_, ob_, field_flag_),
        record_process(false),
        cg_tol(1e-6),
        stag_tol(0.025),
        GN_iter(10),
        Lp(2),
        epsilon2(1e-9) {}

  void invert() override;

  void record_every_iteration() { this->record_process = true; }

  // when cg_iteration_factor is greater than 1, it's the iteration
  // number of conjugate gradient method; when it's lower than 1, CG
  // iteration number is  cg_iteration_factor times the cells number
  void set_CG_parameter(double tol, double cg_iteration_factor);
  void set_stagnation_tolerance(double x) { stag_tol = x; }
  void set_max_GN_iterations(int x) { GN_iter = x; }
  void set_Lp_inversion_parameter(int Lp_, double epsilon2_)
  {
    this->Lp = Lp_;
    this->epsilon2 = epsilon2_;
  }

  void display_inversion_parameters() const;
  // void invert_Gauss_Newton_ada(int maxit1, int maxit2, int maxit3, VectorXd m_ini,
  //                              double tol1, double tol2, double stagnate_tol, double lambda0, int ada_level = 5, double percentage = 0.1, double damp = 0.7, double n = 1);

protected:
  int Lp;
  double epsilon2;
  bool record_process;
  double cg_tol;              // tolerance for conjugate gradient
  double cg_iteration_factor; // CG iteration number
  double stag_tol;            // stagnation
  int GN_iter;                // Gauss-Newton iteration number
};

#endif