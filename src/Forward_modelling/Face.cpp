#include "Face.h"
Face::Face() {
  xc = 0;
  yc = 0;
  zc = 0;
  for (int i = 0; i < 2; i++) {
    neigh_cells[i] = NULL;
  }
  for (int i = 0; i < 4; i++) {
    child_faces[i] = NULL;
  }
  level = 0;
  isleaf = true;
}
Face::Face(Cell* c1, Cell* c2, unsigned int direction, int level_, bool leaf) {
  this->direction = direction;
  xc = 0;
  yc = 0;
  zc = 0;
  this->set_neigh_cells(c1, c2, direction);
  this->level = level_;
  isleaf = leaf;
  for (int i = 0; i < 4; i++) {
    child_faces[i] = NULL;
  }
}
Face::Face(Cell* c1,
           Cell* c2,
           double x,
           double y,
           double z,
           unsigned int direction,
           int level_,
           bool leaf) {
  this->direction = direction;
  xc = x;
  yc = y;
  zc = z;
  this->set_neigh_cells(c1, c2, direction);
  this->level = level_;
  isleaf = leaf;
  if (direction == UP_DOWN) {
    if (c1 != NULL && c2 != NULL) {
      if (!(abs(neigh_cells[0]->_z[1] - zc) < 1e-10)) {
        display();
      }
      assert(abs(neigh_cells[0]->_z[1] - neigh_cells[1]->_z[0]) < 1e-10);
      assert(abs(neigh_cells[0]->_z[1] - zc) < 1e-10);
    }
  }
  for (int i = 0; i < 4; i++) {
    child_faces[i] = NULL;
  }
}

void Face::set_neigh_cells(Cell* c1, Cell* c2, unsigned int normal_direction) {
  this->neigh_cells[0] = c1;
  this->neigh_cells[1] = c2;

  //if c1's coordinates>c2's coordinates,swap
  double x1, y1, z1;
  double x2, y2, z2;
  if (c1 != NULL && c2 != NULL) {
    c1->get_center(x1, y1, z1);
    c2->get_center(x2, y2, z2);
    if (normal_direction == NORTH_SOUTH) {               //x
      if (x1 > x2) {
        neigh_cells[0] = c2;
        neigh_cells[1] = c1;
        neigh_cells[0]->get_center(x1, y1, z1);
        neigh_cells[1]->get_center(x2, y2, z2);
        assert(x1 < x2);
      }
    } else if (normal_direction == WEST_EAST) {          //y
      if (y1 > y2) {
        neigh_cells[0] = c2;
        neigh_cells[1] = c1;

        neigh_cells[0]->get_center(x1, y1, z1);
        neigh_cells[1]->get_center(x2, y2, z2);
        assert(y1 < y2);
      }
    } else if (normal_direction == UP_DOWN) {
      if (z1 > z2) {
        neigh_cells[0] = c2;
        neigh_cells[1] = c1;

        neigh_cells[0]->get_center(x1, y1, z1);
        neigh_cells[1]->get_center(x2, y2, z2);
        assert(z1 < z2);
      }
    } else {
      cout << normal_direction << endl;
      this->display();
      cerr << "not NORTH_SOUTH/WEST_EAST/UP_DOWN" << endl;
      abort();
    }
  }
}

void Face::set_child_faces(Face* f1, Face* f2, Face* f3, Face* f4) {
  this->child_faces[0] = f1;
  this->child_faces[1] = f2;
  this->child_faces[2] = f3;
  this->child_faces[3] = f4;
}
void Face::display() {
  double x, y, z;
  double dx, dy, dz;
  for (int i = 0; i < 2; i++) {
    cout << "Side " << i << ":\t";
    if (neigh_cells[i] == NULL) {
      cout << "NULL" << endl;
    } else {
      cout << "[(" << setprecision(10) << neigh_cells[i]->_x[0] << ", "
           << setprecision(10) << neigh_cells[i]->_x[1] << ")  "
           << "(" << neigh_cells[i]->_y[0] << ", "
           << neigh_cells[i]->_y[1]<< ")  "
           << "(" << neigh_cells[i]->_z[0]<< ", "
           << neigh_cells[i]->_z[1]<< ")]" << endl;
    }
  }
  cout << "center: (" << xc << ", " << yc<< ", "
       << zc<< ")" << endl
       << endl;
}

void Face::set_center(double x, double y, double z) {
  this->xc = x;
  this->yc = y;
  this->zc = z;
}