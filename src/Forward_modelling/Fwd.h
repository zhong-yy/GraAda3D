#ifndef _FWD
#define _FWD
#include <bitset>
#include"GravFormula.h"
#include "gs.h"
#include "Observation.h"
#include "Mesh.h"
using namespace std;

class Fwd
{
public:
  Fwd();
  Fwd(const Mesh& mesh_, const Observation& ob_, unsigned long long field_flag_ = Compute_g_z);
  Fwd(const Mesh& mesh_, const Observation& ob_, bitset<10> field_flag);
  virtual ~Fwd();

  VectorXd compute_gobs(const VectorXd& rho){
    VectorXd d=this->G*rho;
    return d;
  }

  void set_field_flag(unsigned long long field_flag1){
    field_flag=field_flag1;
  }

  void compute_G();

  //bind mesh of tesseroids
  void set_mesh(const Mesh& mesh0);

  //bind computation points
  void set_observation(const Observation& ob0);

  const Eigen::MatrixXd &get_G() { return G; }

  unsigned int get_n_fields()
  {
    return field_flag.count();
  }

protected:
  Eigen::MatrixXd G;
  Mesh mesh;
  Observation ob;
  int Nm;    //number of parameters
  int N_obs; //number of data points
  int Nd;    //number of observations

  bitset<10> field_flag;

  void display_info_fields()const;
  /** 
   * 0 V
   * 1 g_z
   * 2 g_x
   * 3 g_y
   * 4 T_zz
   * 5 T_xz
   * 6 T_yz
   * 7 T_xx
   * 8 T_xy
   * 9 T_yy
   */
};
#endif