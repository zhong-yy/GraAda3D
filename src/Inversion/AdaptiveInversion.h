#ifndef ADAPTIVE_INVERSION_H
#define ADAPTIVE_INVERSION_H
#include <algorithm>
#include <fstream>

#include "InversionBase.h"

// AdaptiveInversion--Adaptive inversion interface
class AdaptiveInversion : public InversionBase
{
protected:
  double refinement_percentage;     // percentage of cells to be refined every time
  int max_refinement_number;        // maximal times of refinement
  int interval_between_refinements; // interval of iterations beween successive
                                    // refinements
  double min_dx;
  double min_dy;
  double min_dz;

public:
  AdaptiveInversion();
  AdaptiveInversion(const Mesh &mesh_,
                    const Observation &ob_,
                    unsigned long long field_flag_);
  virtual ~AdaptiveInversion();

  void set_max_refinement_number(int n) { max_refinement_number = n; }

  void set_refinement_percentage(double x) { refinement_percentage = x; }
  void set_min_cell_size_in_adaptive_mesh(double min_dx0, double min_dy0, double min_dz0)
  {
    this->min_dx = min_dx0;
    this->min_dy = min_dy0;
    this->min_dz = min_dz0;
  }
  void set_interval_between_refinements(int n)
  {
    interval_between_refinements = n;
  }

  // update sensitivity matrix of gravity forward modelling after mesh is
  // refined
  void expand_G(const map<unsigned int, Cell *> &split_cells);

  // calculate refinement index for each grid cell
  void indicator_calculator(VectorXd &indicator);

  void refine_mesh(double a);

  /**
   * @brief refine inversion mesh and constraint model mesh. The constraint
   * model values in the updated mesh are interpolated from the initial
   * pre-defined interpolator which has been built from gridded data from file.
   *
   * @param  a                  refinement percentage
   * @param  interp             interpolator
   * @param  reference_surface  radius of the reference shpere
   * @param  flag               It must be "crg" or "pet", specifying
   * cross-gradient model or reference density model is interpolated for a
   * predefined "interpolator"
   *
   * see: InversionBase class, InversionBase::create_crg_model_from_data(...)
   *      InversionBase::create_ref_model_from_data(...)
   */
  void refine_mesh(double a,
                   InterpMultilinear<3, double> &interp,
                   string flag = "crg");

  void refine_mesh(double a,
                   InterpMultilinear<3, double> &interp_m0s,
                   InterpMultilinear<3, double> &interp_m0);

  void sort_vec(const VectorXd &vec, VectorXd &sorted_vec, VectorXi &ind)
  {
    //[0 1 2 3 ... N-1]
    ind = VectorXi::LinSpaced(vec.size(), 0, vec.size() - 1);

    auto rule = [vec](int i, int j) -> bool
    {
      return vec(i) > vec(j);
    };
    std::sort(ind.data(), ind.data() + ind.size(), rule);
    // data() return the pointer to the first element of VectorXd, which is similar to begin()
    sorted_vec.resize(vec.size());
    for (int i = 0; i < vec.size(); i++)
    {
      sorted_vec(i) = vec(ind(i));
    }
  }
};

#endif