#include "Cell.h"
Cell::Cell(double x0,
           double y0,
           double z0,
           double x1,
           double y1,
           double z1,
           int level_,
           int num_para,
           bool isleaf_)
    : RectPrism(x0, y0, z0, x1, y1, z1)
{
  this->id = -1;
  this->parameters.resize(num_para);
  for (int i = 0; i < num_para; i++)
  {
    this->parameters[i] = 0;
  }

  this->level = level_;
  this->isleaf = isleaf_;

  for (int i = 0; i < 8; i++)
  {
    child_cells[i] = NULL;
  }
  for (int i = 0; i < 2; i++)
  {
    external_faces_x[i] = NULL;
    external_faces_y[i] = NULL;
    external_faces_z[i] = NULL;
  }
  for (int i = 0; i < 4; i++)
  {
    internal_faces_x[i] = NULL;
    internal_faces_y[i] = NULL;
    internal_faces_z[i] = NULL;
  }
}

void Cell::set_internal_faces(Face *f[4], unsigned int normal_dirction)
{
  if (normal_dirction == NORTH_SOUTH)
  {
    // the face is perpendicular to the x axis
    assert(abs(f[0]->xc - f[1]->xc) < TOL && abs(f[0]->xc - f[2]->xc) < TOL &&
           abs(f[0]->xc - f[3]->xc) < TOL);
    for (int i = 0; i < 3; i++)
    {
      for (int j = 1; j < 4 - i; j++)
      { // bubble sort
        if (f[j - 1]->yc > f[j]->yc)
        {
          Face *temp = f[j - 1];
          f[j - 1] = f[j];
          f[j] = temp;
        }
      }
    }
    if (f[0]->zc > f[1]->zc)
    {
      Face *temp = f[0];
      f[0] = f[1];
      f[1] = temp;
    }
    if (f[2]->zc > f[3]->zc)
    {
      Face *temp = f[2];
      f[2] = f[3];
      f[3] = temp;
    }
    for (int i = 0; i < 4; i++)
    {
      this->internal_faces_x[i] = f[i];
    }
  }
  else if (normal_dirction == WEST_EAST)
  {
    // the face is perpendicular to the y axis
    assert(abs(f[0]->yc - f[1]->yc) < TOL && abs(f[0]->yc - f[2]->yc) < TOL &&
           abs(f[0]->yc - f[3]->yc) < TOL);
    for (int i = 0; i < 3; i++)
    {
      for (int j = 1; j < 4 - i; j++)
      {
        if (f[j - 1]->xc > f[j]->xc)
        {
          Face *temp = f[j - 1];
          f[j - 1] = f[j];
          f[j] = temp;
        }
      }
    }
    if (f[0]->zc > f[1]->zc)
    {
      Face *temp = f[0];
      f[0] = f[1];
      f[1] = temp;
    }
    if (f[2]->zc > f[3]->zc)
    {
      Face *temp = f[2];
      f[2] = f[3];
      f[3] = temp;
    }
    for (int i = 0; i < 4; i++)
    {
      this->internal_faces_y[i] = f[i];
    }
  }
  else if (normal_dirction == UP_DOWN)
  {
    // the face is perpendicular to the z axis
    assert(abs(f[0]->zc - f[1]->zc) < TOL && abs(f[0]->zc - f[2]->zc) < TOL &&
           abs(f[0]->zc - f[3]->zc) < TOL);
    for (int i = 0; i < 3; i++)
    {
      for (int j = 1; j < 4 - i; j++)
      {
        if (f[j - 1]->xc > f[j]->xc)
        {
          Face *temp = f[j - 1];
          f[j - 1] = f[j];
          f[j] = temp;
        }
      }
    }
    if (f[0]->yc > f[1]->yc)
    {
      Face *temp = f[0];
      f[0] = f[1];
      f[1] = temp;
    }
    if (f[2]->yc > f[3]->yc)
    {
      Face *temp = f[2];
      f[2] = f[3];
      f[3] = temp;
    }
    for (int i = 0; i < 4; i++)
    {
      this->internal_faces_z[i] = f[i];
    }
  }
  else
  {
    cerr << "Not NORTH_SOUTH/WEST_EAST/UP_DOWN" << endl;
  }
}
void Cell::set_external_faces(Face *f1,
                              Face *f2,
                              unsigned int normal_dirction)
{
  switch (normal_dirction)
  {
  case NORTH_SOUTH: // x
    external_faces_x[0] = f1;
    external_faces_x[1] = f2;
    if (f1->xc > f2->xc)
    {
      // cout << "x1>x2, f1: (" << f1->thetac << ", " << f1->phic << ", " <<
      // f1->rc << ")  f2: (" << f2->thetac << ", " << f2->phic << ", " <<
      // f2->rc << ")" << endl;
      external_faces_x[0] = f2;
      external_faces_x[1] = f1;
    }
    break;
  case WEST_EAST: // y
    external_faces_y[0] = f1;
    external_faces_y[1] = f2;
    if (f1->yc > f2->yc)
    {
      external_faces_y[0] = f2;
      external_faces_y[1] = f1;
    }
    break;
  case UP_DOWN: // z
    external_faces_z[0] = f1;
    external_faces_z[1] = f2;
    if (f1->zc > f2->zc)
    {
      external_faces_z[0] = f2;
      external_faces_z[1] = f1;
    }
    break;
  default:
    cerr << "Not NORTH_SOUTH/WEST_EAST/UP_DOWN" << endl;
  }
}
