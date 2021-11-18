#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include "GaussNewtonInversion.h"
#include "timer.h"
class GraAdaInv {
 public:
  GraAdaInv();
  ~GraAdaInv();

  void read_data_parameters(string data_para);
  void read_model_parameters(string model_para);
  void read_inversion_parameters(string inversion_para);
  void set_config_file(string file_name) { this->config_file = file_name; }
  void start_inversion();
  void write_result();

  int read_data_from_file(string file_name,
                          const vector<unsigned int>& data_order);
  istream& next_valid_line(istream& is, string& str);

  /**
   * @brief  Deal with comments, spaces and empty lines
   *
   * @param line a line of text to be processed
   * @param comment_str a character or a string used to begin a comment
   */
  void line_process(std::string& line, const std::string comment_str = "#");

 protected:
  double height;
  int n_obs;
  Observation ob;
  VectorXd dobs;
  Mesh inv_mesh;
  GaussNewtonInversion* inv;

  int n_fields;
  unsigned long long field_flag;

  string data_file;
  string config_file;

  double Lp_inversion_p;
  double Lp_inversion_eps;
  double min_size_dx, min_size_dy, min_size_dz;
  double x_model[2];
  double y_model[2];
  double z_model[2];
  double n_z, n_x, n_y;

  vector<double> noise_percentage;
  vector<double> equipment_noise;

  double as, az, ax, ay, acrg;
  double depth_weighting_exponent;

  double target_misfit;  // data misfit

  double cg_tol;  // conjugate gradient method tolerance
  double cg_iteration_factor;

  double stagnate_tol;
  int gauss_newton_iterations;

  double start_lambda;
  int n_lambda;
  double decreasing_rate;

  double min_value, max_value;

  double refinement_tol;
  int interval_between_two_refinements;
  int max_refinement;
  int show_process_or_not;

  string crg_model_file;
  string ref_model_file;

  int n_x_crg_model, n_y_crg_model, n_z_crg_model;
  int n_x_pet_model, n_y_pet_model, n_z_pet_model;

  bool use_crg;
  bool use_petr;

  string format_of_coordinates_crg;
  string format_of_coordinates_petr;
  int fast_dimension_crg;
  int fast_dimension_petr;

  string output_model_name;
};

int main(int argc, char** argv) {
  if (argc < 2) {
    printf("Usage: %s input_paramters_filename\n", argv[0]);
    return 1;
  }

  assert(argc == 2);
  ifstream config(argv[1]);
  if (!config.good()) {
    cout << "Cannot read " << argv[1] << ", please check whether " << argv[1]
         << " exists" << endl;
    return 1;
  }
  
  cout<<R"(  ____                     _          _           _____ ____  )"<<endl
      <<R"( / ___|  _ __    __ _     / \      __| |   __ _  |___ /|  _ \ )"<<endl
      <<R"(| |  _  | '__|  / _` |   / _ \    / _` |  / _` |   |_ \| | | |)"<<endl
      <<R"(| |_| | | |    | (_| |  / ___ \  | (_| | | (_| |  ___) | |_| |)"<<endl
      <<R"( \____| |_|     \__,_| /_/   \_\  \__,_|  \__,_| |____/|____/ )"<<endl;

  cout << "Using inversion parameters from configuration file: " << argv[1]
       << endl;
       
  

  GraAdaInv inv_case;

  string data_para;
  string model_para;
  string inversion_para;

  string line;

  inv_case.next_valid_line(config, line);
  istringstream iss1(line);
  iss1 >> data_para;

  inv_case.next_valid_line(config, line);
  istringstream iss2(line);
  iss2 >> model_para;

  inv_case.next_valid_line(config, line);
  istringstream iss3(line);
  iss3 >> inversion_para;

  // this must called before reading model parameters
  inv_case.read_data_parameters(data_para);
  inv_case.read_model_parameters(model_para);
  inv_case.read_inversion_parameters(inversion_para);

  inv_case.start_inversion();
  inv_case.write_result();

  return 0;
}

GraAdaInv::GraAdaInv() {
  inv = NULL;
  use_crg = false;
  use_petr = false;
  format_of_coordinates_crg = "yxz";
  format_of_coordinates_petr = "yxz";
  fast_dimension_crg = 0;
  fast_dimension_petr = 0;
  as = 1;
  az = 1;
  ax = 1;
  ay = 1;
  acrg = 1;
}
GraAdaInv::~GraAdaInv() {
  if (inv != NULL)
    delete inv;
}

void GraAdaInv::line_process(std::string& line, const std::string comment_str) {
  for (char& c : line)  // C++11以上版本的语法
  {
    //制表符 tab，逗号，分号都当作有效的分隔符，统一转成空格
    //为了避免错误，回车符和换行符也转为空格（否则无法处理空行）
    if (c == '\t' || c == ',' || c == ';' || c == '\r' || c == '\n')
      c = ' ';
  }

  //查找注释符所在位置，如果不存在，则得到string::npos
  int n_comment_start = line.find_first_of(comment_str);
  if (n_comment_start != std::string::npos)  //这一句必须的
    line.erase(n_comment_start);             //删除注释

  line.erase(0, line.find_first_not_of(" "));  //删除行首空格
  line.erase(line.find_last_not_of(" ") + 1);  //删除行末空格
  //调用的string& erase (size_t pos = 0, size_t len = npos);
  // len为默认参数,
  // size_t可以当作无符号整数，npos是string内置的静态常量，为size_t的最大值

  /****************************************************************
   *处理完毕。如果这一行只有空格，制表符 tab，注释，那么处理后line为空；
   *如果行首有多个空格(或者空格和tab交错)，行尾为注释，如
   *“   a b c#坐标”
   *那么处理后字符串line的行首多个空格(和tab)和行尾注释被删掉，得到
   *“a b c”
   ****************************************************************/
}

istream& GraAdaInv::next_valid_line(istream& input_stream, string& line) {
  while (std::getline(input_stream, line)) {
    line_process(line);
    if (line.empty()) {
      continue;
    } else {
      break;
    }
  }
  return input_stream;
}

void GraAdaInv::read_data_parameters(string data_para) {
  ifstream input_stream(data_para.c_str());

  string line;

  next_valid_line(input_stream, line);
  istringstream iss(line);
  iss >> data_file;

  next_valid_line(input_stream, line);
  iss.clear();
  iss.str("");
  iss.str(line);
  iss >> height;
  cout << "Height of data" << height << endl;

  next_valid_line(input_stream, line);
  iss.clear();
  iss.str("");
  iss.str(line);

  iss >> n_fields;
  cout << endl;
  cout << n_fields << " gravity component"
       << ((n_fields == 1) ? (" is") : ("s are")) << " used" << endl;
  // dobs.resize(n_obs * n_fileds);

  field_flag = 0;
  vector<unsigned int> field_label(n_fields);
  next_valid_line(input_stream, line);
  iss.clear();
  iss.str("");
  iss.str(line);
  for (int i = 0; i < n_fields; i++) {
    iss >> field_label[i];
    switch (field_label[i]) {
      case 0:
        field_flag = field_flag | Compute_V;
        break;
      case 1:
        field_flag = field_flag | Compute_g_z;
        break;
      case 2:
        field_flag = field_flag | Compute_g_x;
        break;
      case 3:
        field_flag = field_flag | Compute_g_y;
        break;
      case 4:
        field_flag = field_flag | Compute_T_zz;
        break;
      case 5:
        field_flag = field_flag | Compute_T_zx;
        break;
      case 6:
        field_flag = field_flag | Compute_T_zy;
        break;
      case 7:
        field_flag = field_flag | Compute_T_xx;
        break;
      case 8:
        field_flag = field_flag | Compute_T_xy;
        break;
      case 9:
        field_flag = field_flag | Compute_T_yy;
        break;
    }
  }

  vector<unsigned int> data_order;
  data_order.resize(n_fields);
  for (int i = 0; i < n_fields; i++) {
    data_order[i] = i;
  }

  auto sort_rule = [field_label](unsigned int i, unsigned int j) -> bool {
    return field_label[i] < field_label[j];
  };
  sort(data_order.begin(), data_order.end(), sort_rule);

  // double noise_percentage, equipment_noise;
  // config >> noise_percentage >> equipment_noise;

  noise_percentage.resize(n_fields);
  equipment_noise.resize(n_fields);

  vector<double> temp1 = noise_percentage;
  vector<double> temp2 = equipment_noise;

  for (int i = 0; i < n_fields; i++) {
    next_valid_line(input_stream, line);
    iss.clear();
    iss.str("");
    iss.str(line);
    iss >> temp1[i] >> temp2[i];
    // cout << temp1[i] << "," << temp2[i] << endl;
  }

  for (int i = 0; i < n_fields; i++) {
    noise_percentage[i] = temp1[data_order[i]];
    equipment_noise[i] = temp2[data_order[i]];
    // cout << noise_percentage[i] << "," << equipment_noise[i] << endl;
  }

  read_data_from_file(data_file, data_order);
}

int GraAdaInv::read_data_from_file(string data_file,
                                   const vector<unsigned int>& data_order) {
  ifstream input_stream(data_file.c_str());
  assert(input_stream.good());
  string line;
  int n_obs = 0;  // number of observation points
  vector<vector<double> > g_data;
  g_data.resize(n_fields);
  assert(data_order.size() == n_fields);
  while (std::getline(input_stream, line)) {
    line_process(line);
    if (line.empty()) {
      continue;
    } else {
      n_obs++;
      double x0, y0, z0;
      std::istringstream iss(line);
      iss >> x0 >> y0;
      ob.add_point(x0, y0, height);
      // double gx, gy, gz;
      for (int j = 0; j < n_fields; j++) {
        double tmp;
        iss >> tmp;
        g_data[j].push_back(tmp);
        // input_stream >> dobs(n_line + j * n_obs);
      }
    }
  }
  dobs.resize(n_obs * n_fields);
  for (int j = 0; j < n_fields; j++) {
    assert(g_data[j].size() == n_obs);
    for (int i = 0; i < n_obs; i++) {
      dobs(i + j * n_obs) = g_data[data_order[j]][i];
    }
  }
  cout << "Number of observation points: " << n_obs << endl;
  return n_obs;
}

void GraAdaInv::read_model_parameters(string model_para) {
  ifstream input_stream(model_para.c_str());
  string line;
  next_valid_line(input_stream, line);
  istringstream iss(line);
  iss >> x_model[0] >> x_model[1] >> n_x;

  next_valid_line(input_stream, line);
  istringstream iss2(line);
  iss2 >> y_model[0] >> y_model[1] >> n_y;

  next_valid_line(input_stream, line);
  istringstream iss3(line);
  iss3 >> z_model[0] >> z_model[1] >> n_z;
  cout << "Model extension in x direction (m): ";
  cout << x_model[0] << ", " << x_model[1] << endl;
  cout << "Model extension in y direction (m): ";
  cout << y_model[0] << ", " << y_model[1] << endl;
  cout << "Model extension in z direction (m):: ";
  cout << z_model[0] << ", " << z_model[1] << endl;

  int n_pad_x;
  int n_pad_y;
  int n_pad_z;
  double factor_x;
  double factor_y;
  double factor_z;

  next_valid_line(input_stream, line);
  iss.clear();
  iss.str("");
  iss.str(line);
  iss >> n_pad_x;

  next_valid_line(input_stream, line);
  iss.clear();
  iss.str("");
  iss.str(line);
  iss >> factor_x;

  next_valid_line(input_stream, line);
  iss.clear();
  iss.str("");
  iss.str(line);
  iss >> n_pad_y;

  next_valid_line(input_stream, line);
  iss.clear();
  iss.str("");
  iss.str(line);
  iss >> factor_y;

  next_valid_line(input_stream, line);
  iss.clear();
  iss.str("");
  iss.str(line);
  iss >> n_pad_z;

  next_valid_line(input_stream, line);
  iss.clear();
  iss.str("");
  iss.str(line);
  iss >> factor_z;

  if (n_pad_x == 0 && n_pad_y == 0 && n_pad_z == 0) {
    cout << "No padding cells" << endl;
  } else {
    cout<<"Padding is introduced around the region of interest."<<endl;
    cout << "Number of padding cells in x direction: " << n_pad_x << endl;
    cout << "Number of padding cells in y direction: " << n_pad_y << endl;
    cout << "Number of padding cells in z direction: " << n_pad_z << endl;

    cout << "Increasing factor of padding cells (x): " << factor_x << endl;
    cout << "Increasing factor of padding cells (y): " << factor_y << endl;
    cout << "Increasing factor of padding cells (z): " << factor_z << endl;
  }

  inv_mesh.generate_regular_mesh(x_model, n_x, y_model, n_y, z_model, n_z);
  inv_mesh.generate_regular_mesh_with_padding(
      x_model, n_x, y_model, n_y, z_model, n_z, n_pad_x, factor_x, n_pad_y,
      factor_y, n_pad_z, factor_z);

  inv_mesh.out_model_vtk("initial_mesh.vtk");
}

void GraAdaInv::read_inversion_parameters(string inversion_para) {
  ifstream input_stream(inversion_para.c_str());

  string line;
  next_valid_line(input_stream, line);
  istringstream iss(line);
  iss >> as >> az >> ax >> ay >> acrg;

  next_valid_line(input_stream, line);
  iss.clear();
  iss.str("");
  iss.str(line);
  iss >> Lp_inversion_p >> Lp_inversion_eps;

  next_valid_line(input_stream, line);
  iss.clear();
  iss.str("");
  iss.str(line);
  iss >> depth_weighting_exponent;

  next_valid_line(input_stream, line);
  iss.clear();
  iss.str("");
  iss.str(line);
  iss >> start_lambda >> n_lambda >> decreasing_rate;

  //   double final_lambda, final_misfit;
  next_valid_line(input_stream, line);
  iss.clear();
  iss.str("");
  iss.str(line);
  iss >> target_misfit;

  next_valid_line(input_stream, line);
  iss.clear();
  iss.str("");
  iss.str(line);
  iss >> stagnate_tol;

  next_valid_line(input_stream, line);
  iss.clear();
  iss.str("");
  iss.str(line);
  iss >> cg_tol >> cg_iteration_factor;

  next_valid_line(input_stream, line);
  iss.clear();
  iss.str("");
  iss.str(line);
  iss >> gauss_newton_iterations;

  next_valid_line(input_stream, line);
  iss.clear();
  iss.str("");
  iss.str(line);
  iss >> min_value >> max_value;

  next_valid_line(input_stream, line);
  iss.clear();
  iss.str("");
  iss.str(line);
  iss >> max_refinement;

  next_valid_line(input_stream, line);
  iss.clear();
  iss.str("");
  iss.str(line);
  iss >> min_size_dx >> min_size_dy >> min_size_dz;

  next_valid_line(input_stream, line);
  iss.clear();
  iss.str("");
  iss.str(line);
  iss >> interval_between_two_refinements;

  next_valid_line(input_stream, line);
  iss.clear();
  iss.str("");
  iss.str(line);
  iss >> refinement_tol;

  next_valid_line(input_stream, line);
  iss.clear();
  iss.str("");
  iss.str(line);
  iss >> show_process_or_not;

  next_valid_line(input_stream, line);
  iss.clear();
  iss.str("");
  iss.str(line);
  int flag = 0;
  iss >> flag;
  if (flag == 0) {
    use_crg = false;
  } else {
    use_crg = true;
  }

  next_valid_line(input_stream, line);
  iss.clear();
  iss.str("");
  iss.str(line);
  if (use_crg) {
    iss >> crg_model_file >> format_of_coordinates_crg;
    if (format_of_coordinates_crg == "xyz") {
      iss >> n_x_crg_model >> n_y_crg_model >> n_z_crg_model;
    } else if (format_of_coordinates_crg == "yxz") {
      iss >> n_y_crg_model >> n_x_crg_model >> n_z_crg_model;
    } else if (format_of_coordinates_crg == "zxy") {
      iss >> n_z_crg_model >> n_x_crg_model >> n_y_crg_model;
    } else if (format_of_coordinates_crg == "zyx") {
      iss >> n_z_crg_model >> n_y_crg_model >> n_x_crg_model;
    } else if (format_of_coordinates_crg == "xzy") {
      iss >> n_x_crg_model >> n_z_crg_model >> n_y_crg_model;
    } else if (format_of_coordinates_crg == "yzx") {
      iss >> n_y_crg_model >> n_z_crg_model >> n_x_crg_model;
    }
    iss >> fast_dimension_crg;
  }

  next_valid_line(input_stream, line);
  iss.clear();
  iss.str("");
  iss.str(line);
  flag = 0;
  iss >> flag;
  if (flag == 0) {
    use_petr = false;
  } else {
    use_petr = true;
  }

  next_valid_line(input_stream, line);
  iss.clear();
  iss.str("");
  iss.str(line);
  if (use_petr) {
    iss >> ref_model_file >> format_of_coordinates_petr;
    if (format_of_coordinates_petr == "xyz") {
      iss >> n_x_pet_model >> n_y_pet_model >> n_z_pet_model;
    } else if (format_of_coordinates_petr == "yxz") {
      iss >> n_y_pet_model >> n_x_pet_model >> n_z_pet_model;
    } else if (format_of_coordinates_petr == "zxy") {
      iss >> n_z_pet_model >> n_x_pet_model >> n_y_pet_model;
    } else if (format_of_coordinates_petr == "zyx") {
      iss >> n_z_pet_model >> n_y_pet_model >> n_x_pet_model;
    } else if (format_of_coordinates_petr == "xzy") {
      iss >> n_x_pet_model >> n_z_pet_model >> n_y_pet_model;
    } else if (format_of_coordinates_petr == "yzx") {
      iss >> n_y_pet_model >> n_z_pet_model >> n_x_pet_model;
    }
    iss >> fast_dimension_petr;
  }

  next_valid_line(input_stream, line);
  iss.clear();
  iss.str("");
  iss.str(line);
  iss >> output_model_name;
}
void GraAdaInv::start_inversion() {
  string line;

  Timer timer;
  timer.start();
  inv = new GaussNewtonInversion(inv_mesh, ob, field_flag);
  inv->set_dobs(dobs, noise_percentage, equipment_noise);
  inv->set_depth_weighting(depth_weighting_exponent);
  inv->set_weights_of_objectives(as, az, ax, ay, acrg);
  inv->set_max_lambda(start_lambda);
  inv->set_max_GN_iterations(gauss_newton_iterations);
  inv->set_lambda_decreasing_rate(decreasing_rate);
  inv->set_n_lambda(n_lambda);
  inv->set_CG_parameter(cg_tol, cg_iteration_factor);
  inv->set_interval_between_refinements(interval_between_two_refinements);
  inv->set_refinement_percentage(refinement_tol);
  inv->set_max_refinement_number(max_refinement);
  inv->set_min_cell_size_in_adaptive_mesh(min_size_dx, min_size_dy,
                                          min_size_dz);
  inv->set_stagnation_tolerance(stagnate_tol);
  inv->set_target_misfit(target_misfit);

  inv->set_Lp_inversion_parameter(Lp_inversion_p, Lp_inversion_eps);

  if (use_crg) {
    inv->create_crg_model_from_data(
        crg_model_file, n_x_crg_model, n_y_crg_model, n_z_crg_model,
        format_of_coordinates_crg, fast_dimension_crg);
  }
  if (use_petr) {
    inv->create_ref_model_from_data(
        ref_model_file, n_x_pet_model, n_y_pet_model, n_z_pet_model,
        format_of_coordinates_petr, fast_dimension_petr);
  }

  // if()

  inv->compute_G();

  int Nm = inv_mesh.n_elems();
  VectorXd m_min = VectorXd::Constant(Nm, min_value);
  VectorXd m_max = VectorXd::Constant(Nm, max_value);
  VectorXd m_ini = VectorXd::Constant(Nm, 0);

  inv->set_min_max(m_min, m_max);
  inv->set_m_ini(m_ini);

  cout << endl;

  cout << "Start inversion" << endl;

  double computation_time;

  VectorXd m_result, d_pre;

  if (show_process_or_not != 0) {
    // cout << "The inversion process will be recorded" << endl;
    inv->record_every_iteration();
  } else {
    // cout << "The inversion process will not be recorded" << endl;
  }

  inv->display_inversion_parameters();

  inv->invert();
  timer.stop();
  computation_time = timer.getElapsedTimeInSec();
  ofstream out_info("info");
  out_info << "1. Initial inverion mesh" << endl;
  out_info << "X range  ";
  out_info << x_model[0] << ", " << x_model[1] << endl;
  out_info << "Y range: ";
  out_info << y_model[0] << ", " << y_model[1] << endl;
  out_info << "Z range ";
  out_info << z_model[0] << ", " << z_model[1] << endl;
  out_info << "Number of elements along z axis in the initial mesh " << n_z
           << endl;
  out_info << "Number of elements along x axis in the initial mesh " << n_x
           << endl;
  out_info << "Number of elements along y axis in the initial mesh " << n_y
           << endl;
  out_info << "Total number of elements in the initial mesh " << n_z * n_x * n_y
           << endl
           << endl;
  out_info << "2. Final mesh" << endl;
  out_info << "Element number " << inv_mesh.cells.size() << endl << endl;

  out_info << "3. Evaluation of the inversion" << endl;
  out_info << "Final misfit " << inv->get_final_misfit() << endl;
  out_info << "Running time " << computation_time << " seconds" << endl;

  cout << "Running time of inversion:" << computation_time << "s" << endl;

  cout << "Inversion completed" << endl;
}

void GraAdaInv::write_result() {
  cout << "Writing the inversion model into file ..." << endl;
  inv->result2vtk(output_model_name);
  inv->result2text(output_model_name);
  

#ifdef USE_NETCDF
  inv->result2netcdf(output_model_name);
#endif

  cout << "Writing predicted data ..." << endl;
  inv->output_obs_data("dobs");
  inv->output_predicted_data("dpredicted");
  inv->output_obs_data_vtk("dobs.vtk");
  inv->output_predicted_data_vtk("dpredicted.vtk");
}
