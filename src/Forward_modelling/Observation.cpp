#include "Observation.h"

Observation::Observation() : obs(0)
{
  n_obs = 0;
}

Observation::~Observation() {}

void Observation::read_site(const string &name)
{
  ifstream is;
  is.open(name.c_str());
  assert(is.good());
  is >> n_obs;
  cout << "\nReading observation sites ...\nThe total number of sites is "
       << n_obs << ".\n";
  unsigned lab;
  obs.resize(n_obs);
  double x, y, z;
  for (unsigned i = 0; i < n_obs; i++)
  {
    assert(is.good());
    is >> lab;
    // assert(lab == i + 1);
    is >> x >> y >> z;
    obs[i].set_xyz(x, y, z);
  }
  cout << "Reading sites done!\n\n";
}

void Observation::add_point(const Point &p)
{
  obs.push_back(p);
  n_obs++;
}
void Observation::add_point(double x, double y, double z)
{
  obs.push_back(Point(x, y, z));
  n_obs++;
}
void Observation::add_point_pre(double x, double y, double z)
{
  obs.insert(obs.begin(), Point(x, y, z));
  n_obs++;
}

void Observation::generate(int nx,
                           double start_x,
                           double end_x,
                           int ny,
                           double start_y,
                           double end_y,
                           double z0)
{
  obs.clear();
  obs.resize(nx * ny);
  double xspace = 0;
  double yspace = 0;
  if (nx > 1)
  {
    xspace = (end_x - start_x) / (nx - 1);
  }
  if (ny > 1)
  {
    yspace = (end_y - start_y) / (ny - 1);
  }
  for (int i_y = 0; i_y < ny; i_y++)
  {
    for (int i_x = 0; i_x < nx; i_x++)
    {
      obs[i_y * nx + i_x].set_xyz(start_x + i_x * xspace,
                                  start_y + i_y * yspace, z0);
    }
  }
  this->n_obs = nx * ny;
}
void Observation::generate(int n_x,
                           double start_x,
                           double end_x,
                           int n_y,
                           double start_y,
                           double end_y,
                           int n_z,
                           double start_z,
                           double end_z)
{
  obs.clear();
  obs.resize(n_z * n_x * n_y);
  double z_space;
  if (n_z == 1)
  {
    z_space = 0;
  }
  else
  {
    z_space = (end_z - start_z) / (n_z - 1.0);
  }
  double x_space = (end_x - start_x) / (n_x - 1.0);
  double y_space = (end_y - start_y) / (n_y - 1.0);
  for (int i_z = 0; i_z < n_z; i_z++)
  {
    for (int i_y = 0; i_y < n_y; i_y++)
    {
      for (int i_x = 0; i_x < n_x; i_x++)
      {
        obs[i_z * (n_x * n_y) + i_y * n_x + i_x].set_xyz(
            start_x + i_x * x_space,
            start_y + i_y * y_space, start_z + i_z * z_space);
      }
    }
  }
  this->n_obs = n_z * n_x * n_y;
}

ostream &operator<<(ostream &output, Observation &observation)
{
  output << "#Number of observation points: " << observation.obs.size() << '\n';
  output << "#r\ttheta\tphi\n";
  for (int i = 0; i < observation.n_obs; i++)
  {
    output << '\n';
    double *xyz = observation.obs[i].get_xyz();
    output << left << setw(15) << xyz[0] << left << setw(15)
           << xyz[1] << left << setw(15)
           << xyz[2];
  }
  return output;
}
