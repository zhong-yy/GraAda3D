#include "Mesh.h"
Mesh::Mesh()
    : leaf_cells(0),
      leaf_faces(0),
      cells(0),
      faces(0),
      num_cell(0),
      num_face(0) {
  nz = 0;
  nx = 0;
  ny = 0;
  num_leaf_cells = 0;
  num_leaf_faces = 0;
  n_parameters = 1;
  x_lim[0] = 0;
  x_lim[1] = 0;
  y_lim[0] = 0;
  y_lim[1] = 0;
  z_lim[0] = 0;
  z_lim[1] = 0;
}
Mesh::~Mesh() {
  this->clear_all();
}
Mesh::Mesh(const Mesh& source_mesh) {
  this->clear_all();  // clear

  // copy non-pointer type data
  this->nz = source_mesh.nz;
  this->nx = source_mesh.nx;
  this->ny = source_mesh.ny;
  for (int i = 0; i < 2; i++) {
    this->x_lim[i] = source_mesh.x_lim[i];
    this->y_lim[i] = source_mesh.y_lim[i];
    this->z_lim[i] = source_mesh.z_lim[i];
  }
  this->num_face = source_mesh.num_face;
  this->num_cell = source_mesh.num_cell;
  this->z_points = source_mesh.z_points;
  this->num_leaf_cells = source_mesh.num_leaf_cells;
  this->num_leaf_faces = source_mesh.num_leaf_faces;
  this->n_parameters = source_mesh.n_parameters;

  // cout<<"A"<<endl;
  // deep copy pointers
  unordered_map<Cell*, Cell*> cc_mp;  // key: old mesh, value: new mesh
  unordered_map<Face*, Face*> ff_mp;

  this->cells.resize(source_mesh.cells.size());
  for (int i = 0; i < source_mesh.cells.size(); i++) {
    this->cells[i].resize(source_mesh.cells[i].size());
    for (int j = 0; j < source_mesh.cells[i].size(); j++) {
      Cell* c = source_mesh.cells[i][j];
      this->cells[i][j] =
          new Cell(c->_x[0], c->_y[0], c->_z[0], c->_x[1], c->_y[1], c->_z[1],
                   c->level, c->parameters.size(), c->isleaf);
      this->cells[i][j]->set_id(c->id);
      // this->cells[i][j]->_density = c->_density;
      cc_mp[c] = this->cells[i][j];
    }
  }

  // cout<<"A"<<endl;

  this->faces.resize(source_mesh.faces.size());
  for (int i = 0; i < source_mesh.faces.size(); i++) {
    this->faces[i].resize(source_mesh.faces[i].size());
    for (int j = 0; j < source_mesh.faces[i].size(); j++) {
      Face* f = source_mesh.faces[i][j];
      this->faces[i][j] = new Face(NULL, NULL, f->xc, f->yc, f->zc,
                                   f->direction, f->level, f->isleaf);
      ff_mp[f] = this->faces[i][j];
    }
  }

  // set associated pointers
  for (int i = 0; i < cells.size(); i++) {
    for (int j = 0; j < cells[i].size(); j++) {
      Cell* c = source_mesh.cells[i][j];
      for (int k = 0; k < 8; k++) {
        if (c->child_cells[k] != NULL) {
          cc_mp[c]->child_cells[k] = cc_mp[c->child_cells[k]];
        } else {
          cc_mp[c]->child_cells[k] = NULL;
        }
      }
      // cout<<"X"<<endl;
      for (int k = 0; k < 2; k++) {
        if (c->external_faces_z[k] != NULL) {
          cc_mp[c]->external_faces_z[k] = ff_mp[c->external_faces_z[k]];
        } else {
          cc_mp[c]->external_faces_z[k] = NULL;
        }

        if (c->external_faces_x[k] != NULL) {
          cc_mp[c]->external_faces_x[k] = ff_mp[c->external_faces_x[k]];
        } else {
          cc_mp[c]->external_faces_x[k] = NULL;
        }

        if (c->external_faces_y[k] != NULL) {
          cc_mp[c]->external_faces_y[k] = ff_mp[c->external_faces_y[k]];
        } else {
          cc_mp[c]->external_faces_y[k] = NULL;
        }
      }
      for (int k = 0; k < 4; k++) {
        if (c->internal_faces_z[k] != NULL) {
          cc_mp[c]->internal_faces_z[k] = ff_mp[c->internal_faces_z[k]];
        } else {
          cc_mp[c]->internal_faces_z[k] = NULL;
        }

        if (c->internal_faces_x[k] != NULL) {
          cc_mp[c]->internal_faces_x[k] = ff_mp[c->internal_faces_x[k]];
        } else {
          cc_mp[c]->internal_faces_x[k] = NULL;
        }

        if (c->internal_faces_y[k] != NULL) {
          cc_mp[c]->internal_faces_y[k] = ff_mp[c->internal_faces_y[k]];
        } else {
          cc_mp[c]->internal_faces_y[k] = NULL;
        }
      }
    }
  }

  for (int i = 0; i < faces.size(); i++) {
    for (int j = 0; j < faces[i].size(); j++) {
      Face* f = source_mesh.faces[i][j];
      for (int k = 0; k < 4; k++) {
        if (f->child_faces[k] != NULL) {
          ff_mp[f]->child_faces[k] = ff_mp[f->child_faces[k]];
        } else {
          ff_mp[f]->child_faces[k] = NULL;
        }
      }
      for (int k = 0; k < 2; k++) {
        if (f->neigh_cells[k] != NULL) {
          ff_mp[f]->neigh_cells[k] = cc_mp[f->neigh_cells[k]];
        } else {
          ff_mp[f]->neigh_cells[k] = NULL;
        }
      }
    }
  }

  leaf_cells.resize(source_mesh.leaf_cells.size());
  leaf_faces.resize(source_mesh.leaf_faces.size());
  for (int i = 0; i < leaf_cells.size(); i++) {
    Cell* c = source_mesh.leaf_cells[i];
    this->leaf_cells[i] = cc_mp[c];
  }
  for (int i = 0; i < leaf_faces.size(); i++) {
    Face* f = source_mesh.leaf_faces[i];
    this->leaf_faces[i] = ff_mp[f];
  }
}

Mesh& Mesh::operator=(const Mesh& source_mesh) {
  this->clear_all();  // clear

  // copy non-pointer type data
  this->nz = source_mesh.nz;
  this->nx = source_mesh.nx;
  this->ny = source_mesh.ny;
  for (int i = 0; i < 2; i++) {
    this->x_lim[i] = source_mesh.x_lim[i];
    this->y_lim[i] = source_mesh.y_lim[i];
    this->z_lim[i] = source_mesh.z_lim[i];
  }
  this->num_face = source_mesh.num_face;
  this->num_cell = source_mesh.num_cell;
  this->z_points = source_mesh.z_points;
  this->num_leaf_cells = source_mesh.num_leaf_cells;
  this->num_leaf_faces = source_mesh.num_leaf_faces;
  this->n_parameters = source_mesh.n_parameters;

  // cout<<"A"<<endl;
  // deep copy pointers
  unordered_map<Cell*, Cell*> cc_mp;  // key: old mesh, value: new mesh
  unordered_map<Face*, Face*> ff_mp;

  this->cells.resize(source_mesh.cells.size());
  for (int i = 0; i < source_mesh.cells.size(); i++) {
    this->cells[i].resize(source_mesh.cells[i].size());
    for (int j = 0; j < source_mesh.cells[i].size(); j++) {
      Cell* c = source_mesh.cells[i][j];
      this->cells[i][j] =
          new Cell(c->_x[0], c->_y[0], c->_z[0], c->_x[1], c->_y[1], c->_z[1],
                   c->level, c->parameters.size(), c->isleaf);
      this->cells[i][j]->set_id(c->id);
      // this->cells[i][j]->_density = c->_density;
      cc_mp[c] = this->cells[i][j];
    }
  }

  // cout<<"A"<<endl;

  this->faces.resize(source_mesh.faces.size());
  for (int i = 0; i < source_mesh.faces.size(); i++) {
    this->faces[i].resize(source_mesh.faces[i].size());
    for (int j = 0; j < source_mesh.faces[i].size(); j++) {
      Face* f = source_mesh.faces[i][j];
      this->faces[i][j] = new Face(NULL, NULL, f->xc, f->yc, f->zc,
                                   f->direction, f->level, f->isleaf);
      ff_mp[f] = this->faces[i][j];
    }
  }

  // set associated pointers
  for (int i = 0; i < cells.size(); i++) {
    for (int j = 0; j < cells[i].size(); j++) {
      Cell* c = source_mesh.cells[i][j];
      for (int k = 0; k < 8; k++) {
        if (c->child_cells[k] != NULL) {
          cc_mp[c]->child_cells[k] = cc_mp[c->child_cells[k]];
        } else {
          cc_mp[c]->child_cells[k] = NULL;
        }
      }
      // cout<<"X"<<endl;
      for (int k = 0; k < 2; k++) {
        if (c->external_faces_z[k] != NULL) {
          cc_mp[c]->external_faces_z[k] = ff_mp[c->external_faces_z[k]];
        } else {
          cc_mp[c]->external_faces_z[k] = NULL;
        }

        if (c->external_faces_x[k] != NULL) {
          cc_mp[c]->external_faces_x[k] = ff_mp[c->external_faces_x[k]];
        } else {
          cc_mp[c]->external_faces_x[k] = NULL;
        }

        if (c->external_faces_y[k] != NULL) {
          cc_mp[c]->external_faces_y[k] = ff_mp[c->external_faces_y[k]];
        } else {
          cc_mp[c]->external_faces_y[k] = NULL;
        }
      }
      for (int k = 0; k < 4; k++) {
        if (c->internal_faces_z[k] != NULL) {
          cc_mp[c]->internal_faces_z[k] = ff_mp[c->internal_faces_z[k]];
        } else {
          cc_mp[c]->internal_faces_z[k] = NULL;
        }

        if (c->internal_faces_x[k] != NULL) {
          cc_mp[c]->internal_faces_x[k] = ff_mp[c->internal_faces_x[k]];
        } else {
          cc_mp[c]->internal_faces_x[k] = NULL;
        }

        if (c->internal_faces_y[k] != NULL) {
          cc_mp[c]->internal_faces_y[k] = ff_mp[c->internal_faces_y[k]];
        } else {
          cc_mp[c]->internal_faces_y[k] = NULL;
        }
      }
    }
  }

  for (int i = 0; i < faces.size(); i++) {
    for (int j = 0; j < faces[i].size(); j++) {
      Face* f = source_mesh.faces[i][j];
      for (int k = 0; k < 4; k++) {
        if (f->child_faces[k] != NULL) {
          ff_mp[f]->child_faces[k] = ff_mp[f->child_faces[k]];
        } else {
          ff_mp[f]->child_faces[k] = NULL;
        }
      }
      for (int k = 0; k < 2; k++) {
        if (f->neigh_cells[k] != NULL) {
          ff_mp[f]->neigh_cells[k] = cc_mp[f->neigh_cells[k]];
        } else {
          ff_mp[f]->neigh_cells[k] = NULL;
        }
      }
    }
  }

  leaf_cells.resize(source_mesh.leaf_cells.size());
  leaf_faces.resize(source_mesh.leaf_faces.size());
  for (int i = 0; i < leaf_cells.size(); i++) {
    Cell* c = source_mesh.leaf_cells[i];
    this->leaf_cells[i] = cc_mp[c];
  }
  for (int i = 0; i < leaf_faces.size(); i++) {
    Face* f = source_mesh.leaf_faces[i];
    this->leaf_faces[i] = ff_mp[f];
  }
  return *this;
}

void Mesh::clear_all() {
  for (int i = 0; i < num_cell.size(); i++) {
    for (int j = 0; j < num_cell[i]; j++) {
      delete cells[i][j];
      cells[i][j] = NULL;
    }
  }
  for (int i = 0; i < num_face.size(); i++) {
    for (int j = 0; j < num_face[i]; j++) {
      delete faces[i][j];
      faces[i][j] = NULL;
    }
  }

  // Pointer leaf_cells[i] points to the same Cell object as
  // a certain element of vector<vector<Cell*>>cells do, so it has been
  // deleted through the above loops, and should not be deleted again.
  // So it is with leaf_faces
  for (int i = 0; i < leaf_cells.size(); i++) {
    leaf_cells[i] = NULL;
  }
  for (int j = 0; j < leaf_faces.size(); j++) {
    leaf_faces[j] = NULL;
  }

  cells.clear();
  faces.clear();
  leaf_cells.clear();
  leaf_faces.clear();
  num_leaf_cells = 0;
  num_leaf_faces = 0;
  nz = 0;
  nx = 0;
  ny = 0;
  num_cell.clear();
  num_face.clear();
}
RectPrism& Mesh::get_elem(const unsigned int i) {
  assert(i < this->n_elems());
  assert(leaf_cells[i] != NULL);
  return *(this->leaf_cells[i]);
}
Cell* Mesh::get_element_level_0(unsigned int i_x,
                                unsigned int j_y,
                                unsigned int k_z) {
  assert(i_x < nx);
  assert(j_y < ny);
  assert(k_z < nz);
  return this->cells[0][i_x * ny * nz + j_y * nz + k_z];
}

void Mesh::generate_regular_mesh(double x_model[2],
                                 int n_x,
                                 double y_model[2],
                                 int n_y,
                                 double z_model[2],
                                 int n_z,
                                 int num_para) {
  this->clear_all();

  for (int i = 0; i < 2; i++) {
    this->x_lim[i] = x_model[i];
    this->y_lim[i] = y_model[i];
    this->z_lim[i] = z_model[i];
  }

  this->n_parameters = num_para;
  this->nz = n_z;
  this->nx = n_x;
  this->ny = n_y;

  z_points = VectorXd::LinSpaced(n_z + 1, z_model[0], z_model[1]);
  VectorXd x_points = VectorXd::LinSpaced(n_x + 1, x_model[0], x_model[1]);
  VectorXd y_points = VectorXd::LinSpaced(n_y + 1, y_model[0], y_model[1]);

  cells.push_back(vector<Cell*>(0));
  cells[0].resize(n_z * n_x * n_y);

  faces.push_back(vector<Face*>(0));
  faces[0].resize(3 * n_z * n_x * n_y + n_z * n_x + n_x * n_y + n_z * n_y);

  num_cell.push_back(cells[0].size());
  num_face.push_back(faces[0].size());

  // construct cells
  int id = 0;
  for (int i = 0; i < n_x; i++) {
    for (int j = 0; j < n_y; j++) {
      for (int k = 0; k < n_z; k++) {
        cells[0][id] =
            new Cell(x_points(i), y_points(j), z_points(k), x_points(i + 1),
                     y_points(j + 1), z_points(k + 1), 0, num_para, true);
        for (int ip = 0; ip < num_para; ip++) {
          cells[0][id]->set_parameter(0, ip);
        }
        cells[0][id]->set_id(id);
        id++;
      }
    }
  }
  num_leaf_cells = cells[0].size();
  leaf_cells.resize(num_leaf_cells);
  for (int i = 0; i < num_leaf_cells; i++) {
    leaf_cells[i] = cells[0][i];
  }
  // construct faces

  id = 0;
  double xc, yc, zc;
  double dx, dy, dz;
  for (int i = 0; i < n_x; i++) {
    for (int j = 0; j < n_y; j++) {
      for (int k = 0; k < n_z + 1; k++) {
        if (k == 0) {
          faces[0][id] =
              new Face(NULL, cells[0][i * n_y * n_z + j * n_z], UP_DOWN);
          cells[0][i * n_y * n_z + j * n_z]->get_size(dx, dy, dz);
          cells[0][i * n_y * n_z + j * n_z]->get_center(xc, yc, zc);
          faces[0][id]->set_center(xc, yc, zc - 0.5 * dz);
          assert(faces[0][id]->neigh_cells[0] == NULL &&
                 faces[0][id]->neigh_cells[1] != NULL);
        } else if (k == n_z) {
          faces[0][id] = new Face(cells[0][i * n_y * n_z + j * n_z + k - 1],
                                  NULL, UP_DOWN);
          cells[0][i * n_y * n_z + j * n_z + k - 1]->get_size(dx, dy, dz);
          cells[0][i * n_y * n_z + j * n_z + k - 1]->get_center(xc, yc, zc);
          faces[0][id]->set_center(xc, yc, zc + 0.5 * dz);
          assert(faces[0][id]->neigh_cells[0] != NULL &&
                 faces[0][id]->neigh_cells[1] == NULL);
        } else {
          double xc1, yc1, zc1;
          cells[0][i * n_y * n_z + j * n_z + k - 1]->get_center(xc, yc, zc);
          cells[0][i * n_y * n_z + j * n_z + k]->get_center(xc1, yc1, zc1);
          faces[0][id] =
              new Face(cells[0][i * n_y * n_z + j * n_z + k - 1],
                       cells[0][i * n_y * n_z + j * n_z + k], UP_DOWN);
          faces[0][id]->set_center(0.5 * (xc + xc1), 0.5 * (yc + yc1),
                                   0.5 * (zc + zc1));
        }

        if (k > 0) {
          int index = i * n_y * n_z + j * n_z + k - 1;
          cells[0][index]->set_external_faces(faces[0][id - 1], faces[0][id],
                                              UP_DOWN);
        }

        id++;
      }
    }
  }

  // id = 0;
  for (int i = 0; i < n_x + 1; i++) {
    for (int j = 0; j < n_y; j++) {
      for (int k = 0; k < n_z; k++) {
        if (i == 0) {
          faces[0][id] =
              new Face(NULL, cells[0][i * n_y * n_z + j * n_z], NORTH_SOUTH);
          cells[0][i * n_y * n_z + j * n_z]->get_size(dx, dy, dz);
          cells[0][i * n_y * n_z + j * n_z]->get_center(xc, yc, zc);
          faces[0][id]->set_center(xc - 0.5 * dx, yc, zc);
          assert(faces[0][id]->neigh_cells[0] == NULL &&
                 faces[0][id]->neigh_cells[1] != NULL);
        } else if (i == n_x) {
          faces[0][id] = new Face(cells[0][(i - 1) * n_y * n_z + j * n_z + k],
                                  NULL, NORTH_SOUTH);
          cells[0][(i - 1) * n_y * n_z + j * n_z + k]->get_size(dx, dy, dz);
          cells[0][(i - 1) * n_y * n_z + j * n_z + k]->get_center(xc, yc, zc);
          faces[0][id]->set_center(xc + 0.5 * dx, yc, zc);
          assert(faces[0][id]->neigh_cells[0] != NULL &&
                 faces[0][id]->neigh_cells[1] == NULL);
        } else {
          double xc1, yc1, zc1;
          faces[0][id] =
              new Face(cells[0][(i - 1) * n_y * n_z + j * n_z + k],
                       cells[0][i * n_y * n_z + j * n_z + k], NORTH_SOUTH);
          cells[0][(i - 1) * n_y * n_z + j * n_z + k]->get_center(xc, yc, zc);
          cells[0][i * n_y * n_z + j * n_z + k]->get_center(xc1, yc1, zc1);
          faces[0][id]->set_center(0.5 * (xc + xc1), 0.5 * (yc + yc1),
                                   0.5 * (zc + zc1));
        }
        if (i > 0) {
          int index = (i - 1) * n_y * n_z + j * n_z + k;
          cells[0][index]->set_external_faces(faces[0][id - n_y * n_z],
                                              faces[0][id], NORTH_SOUTH);
        }
        id++;
      }
    }
  }

  // id = 0;
  for (int i = 0; i < n_x; i++) {
    for (int j = 0; j < n_y + 1; j++) {
      for (int k = 0; k < n_z; k++) {
        if (j == 0) {
          faces[0][id] =
              new Face(NULL, cells[0][i * n_y * n_z + j * n_z + k], WEST_EAST);
          cells[0][i * n_y * n_z + j * n_z + k]->get_center(xc, yc, zc);
          cells[0][i * n_y * n_z + j * n_z + k]->get_size(dx, dy, dz);
          faces[0][id]->set_center(xc, yc - 0.5 * dy, zc);
          assert(faces[0][id]->neigh_cells[0] == NULL &&
                 faces[0][id]->neigh_cells[1] != NULL);
        } else if (j == n_y) {
          faces[0][id] = new Face(cells[0][i * n_y * n_z + (j - 1) * n_z + k],
                                  NULL, WEST_EAST);
          cells[0][i * n_y * n_z + (j - 1) * n_z + k]->get_center(xc, yc, zc);
          cells[0][i * n_y * n_z + (j - 1) * n_z + k]->get_size(dx, dy, dz);
          faces[0][id]->set_center(xc, yc + 0.5 * dy, zc);
          assert(faces[0][id]->neigh_cells[0] != NULL &&
                 faces[0][id]->neigh_cells[1] == NULL);
        } else {
          double xc1, yc1, zc1;
          faces[0][id] =
              new Face(cells[0][i * n_y * n_z + (j - 1) * n_z + k],
                       cells[0][i * n_y * n_z + j * n_z + k], WEST_EAST);
          cells[0][i * n_y * n_z + (j - 1) * n_z + k]->get_center(xc, yc, zc);
          cells[0][i * n_y * n_z + j * n_z + k]->get_center(xc1, yc1, zc1);
          faces[0][id]->set_center(0.5 * (xc + xc1), 0.5 * (yc + yc1),
                                   0.5 * (zc + zc1));
        }
        if (j > 0) {
          int index = i * n_y * n_z + (j - 1) * n_z + k;
          cells[0][index]->set_external_faces(faces[0][id - 1 * nz],
                                              faces[0][id], WEST_EAST);
        }
        id++;
      }
    }
  }

  // for (int i = 0; i < n_x; i++)
  // {
  //     for (int j = 0; j < n_y; j++)
  //     {
  //         for (int k = 0; k < n_z; k++)
  //         {

  //         }
  //     }
  // }
  num_leaf_faces = faces[0].size();
  leaf_faces.resize(num_leaf_faces);
  for (int i = 0; i < num_leaf_faces; i++) {
    leaf_faces[i] = faces[0][i];
  }

  this->set_n_parameter(num_para);
  // this->sort(0);
}
void Mesh::set_parameter_in_a_region(double xlim[2],double ylim[2],double zlim[2],double para_value,int i_th){
    for (int i = 0; i < cells[0].size(); i++) {
    Cell* c = cells[0][i];
    double cz0 = c->_z[0];
    double cz1 = c->_z[1];
    double cx0 = c->_x[0];
    double cx1 = c->_x[1];
    double cy0 = c->_y[0];
    double cy1 = c->_y[1];

    if ((cz0 > zlim[0] || abs(cz0 - zlim[0]) < 1e-10) &&
        (cz1 < zlim[1] || abs(cz1 - zlim[1]) < 1e-10)) {
      if ((cx0 > xlim[0] || abs(cx0 - xlim[0]) < 1e-10) &&
          (cx1 < xlim[1] || abs(cx1 - xlim[1]) < 1e-10)) {
        if ((cy0 > ylim[0] || abs(cy0 - ylim[0]) < 1e-10) &&
            (cy1 < ylim[1] || abs(cy1 - ylim[1]) < 1e-10)) {
          c->set_parameter(para_value, i_th);
        }
      }
    }
  }
}
void Mesh::set_parameter_in_a_region(double x0,
                                     double x1,
                                     double y0,
                                     double y1,
                                     double z0,
                                     double z1,
                                     double para_value,
                                     int i_th) {
  for (int i = 0; i < cells[0].size(); i++) {
    Cell* c = cells[0][i];
    double cz0 = c->_z[0];
    double cz1 = c->_z[1];
    double cx0 = c->_x[0];
    double cx1 = c->_x[1];
    double cy0 = c->_y[0];
    double cy1 = c->_y[1];

    if ((cz0 > z0 || abs(cz0 - z0) < 1e-10) &&
        (cz1 < z1 || abs(cz1 - z1) < 1e-10)) {
      if ((cx0 > x0 || abs(cx0 - x0) < 1e-10) &&
          (cx1 < x1 || abs(cx1 - x1) < 1e-10)) {
        if ((cy0 > y0 || abs(cy0 - y0) < 1e-10) &&
            (cy1 < y1 || abs(cy1 - y1) < 1e-10)) {
          c->set_parameter(para_value, i_th);
        }
      }
    }
  }
}

void Mesh::set_block_parameter(unsigned int i_min,
                               unsigned int i_max,
                               unsigned int j_min,
                               unsigned int j_max,
                               unsigned int k_min,
                               unsigned int k_max,
                               double para_value,
                               int i_th) {
  assert(i_min >= 0 && i_max < nx);
  assert(j_min >= 0 && j_max < ny);
  assert(k_min >= 0 && k_max < nz);
  for (int i = i_min; i <= i_max; i++) {
    for (int j = j_min; j <= j_max; j++) {
      for (int k = k_min; k <= k_max; k++) {
        int id = i * ny * nz + j * nz + k;
        cells[0][id]->set_parameter(para_value, i_th);
      }
    }
  }
}

void Mesh::out_model_vtk(string filename,
                         int n,
                         vector<string> parameter_name) {
  // prepare the node
  std::set<Point> v_set;
  std::vector<std::vector<Point>> prism(this->n_elems());
  // for (int i = 0; i < n_elems() && (this->get_elem(i)._phi[0] > 0.); i++)
  for (int i = 0; i < n_elems(); i++) {
    // build the point (8 points)
    Point v[8];
    double x1, x2, y1, y2, z1, z2;
    x1 = this->get_elem(i)._x[0];
    x2 = this->get_elem(i)._x[1];
    y1 = this->get_elem(i)._y[0];
    y2 = this->get_elem(i)._y[1];
    z1 = this->get_elem(i)._z[0];
    z2 = this->get_elem(i)._z[1];

    v[0] = Point(x1, y1, z1);
    v[1] = Point(x2, y1, z1);
    v[2] = Point(x2, y2, z1);
    v[3] = Point(x1, y2, z1);
    v[4] = Point(x1, y1, z2);
    v[5] = Point(x2, y1, z2);
    v[6] = Point(x2, y2, z2);
    v[7] = Point(x1, y2, z2);

    for (int j = 0; j < 8; j++)
      v_set.insert(v[j]);
    for (int j = 0; j < 8; j++)
      prism[i].push_back(v[j]);
  }
  std::map<Point, unsigned int> v_id_map;
  unsigned int counter = 0;
  for (std::set<Point>::iterator it = v_set.begin(); it != v_set.end(); it++) {
    v_id_map[(*it)] = counter;
    counter++;
  }
  const unsigned int total_points = v_set.size();
  const unsigned int total_cells = this->n_elems();
  //----------------------------------------------------------------------
  // Open the of stream
  std::ofstream vtk_mesh(filename.c_str());
  if (!vtk_mesh.good()) {
    std::cerr << "Can not open file:\t" << filename + ".vtk" << std::endl;
  } else {
    // Parts 1-2-3, mandatory
    vtk_mesh
        << "# vtk DataFile Version 3.0\n"  // File version and identifier
                                           // Header info, doublely cool data
        << "Gravity inversion using rectangular prisms\n"
        << "ASCII\n";  // ASCII data (not BINARY)

    // Part 4, Geometry/topology, unstructured mesh
    vtk_mesh << "DATASET UNSTRUCTURED_GRID\n";  // topography and geometry

    // POINTS info (0-->n-1)
    vtk_mesh << "\nPOINTS\t" << total_points << "\tdouble\n";
    // Loop POINTS to write out coordinates
    for (std::set<Point>::iterator it = v_set.begin(); it != v_set.end();
         it++) {
      double x = (*it)(0);
      double y = (*it)(1);
      double z = (*it)(2);
      vtk_mesh << x << "\t"   // x-coordinate
               << y << "\t"   // y-coordinate
               << z << "\n";  // z-coordinate
    }

    // CELL info (0-->m-1)
    typedef std::map<Point, unsigned int>::iterator IT;

    vtk_mesh << "\nCELLS\t" << total_cells << "\t" << total_cells * (8 + 1)
             << "\n";
    for (unsigned int i = 0; i < total_cells; i++) {
      std::vector<Point>& T = prism[i];  // 20 vertex
      assert(T.size() == 8);
      unsigned int T_ID[8];
      for (int j = 0; j < 8; j++) {
        IT it = v_id_map.find(T[j]);
        assert(it != v_id_map.end());
        T_ID[j] = (*it).second;
      }
      vtk_mesh << (unsigned int)8 << "\t";
      for (int j = 0; j < 8; j++)
        vtk_mesh << T_ID[j] << "\t";
      vtk_mesh << "\n";
    }

    // CELL types (m)
    vtk_mesh << "\nCELL_TYPES\t" << total_cells << "\n";
    for (unsigned int i = 0; i < total_cells; i++) {
      vtk_mesh << (unsigned int)12 << "\n";  // 12-hexahedron
                                             // figure 3 in vtk format file.
    }

    // Part 5, attributes
    vtk_mesh << "\nCELL_DATA\t" << total_cells << "\n";
    for (int j = 0; j < n; j++) {
      vtk_mesh << "SCALARS " << parameter_name[j] << " double 1\n"
               << "LOOKUP_TABLE "
               << "table" << j << endl;
      for (unsigned int i = 0; i < total_cells; i++) {
        double value = leaf_cells[i]->get_parameter(j);
        vtk_mesh << value << "\n";
      }
    }
    vtk_mesh << "\n";

  }  // file opened successfully

  vtk_mesh.close();
  cout << "The model has been written to vtk file: " << filename << endl;
}

map<unsigned int, Cell*> Mesh::refinement(Cell* c) {
  // cout<<"count="<<std::count(leaf_cells.begin(), leaf_cells.end(), c)<<endl;
  // cout<<"isleaf="<<c->isleaf<<endl;
  assert(std::count(leaf_cells.begin(), leaf_cells.end(), c) == 1);
  bool flag = false;
  bool flag1[2] = {true, true};
  bool flag2[2] = {true, true};
  bool flag3[2] = {true, true};
  int nei_index[2] = {0, 1};
  map<unsigned int, Cell*> split_cells;
  // cout << "AB" << endl;
  do {
    for (int i = 0; i < 2; i++) {
      Face* fx = c->external_faces_x[i];
      if ((fx->neigh_cells[0] != NULL) && (fx->neigh_cells[1] != NULL)) {
        // flag = flag && ((fx->isleaf == false) ||
        // (fx->neigh_cells[nei_index[i]]->level >= c->level));
        flag1[i] = (fx->isleaf == false) ||
                   (fx->neigh_cells[nei_index[i]]->level >= c->level);
      }

      Face* fy = c->external_faces_y[i];
      if ((fy->neigh_cells[0] != NULL) && (fy->neigh_cells[1] != NULL)) {
        // flag = flag && ((fy->isleaf == false) ||
        // fy->neigh_cells[nei_index[i]]->level >= c->level);
        flag2[i] = (fy->isleaf == false) ||
                   fy->neigh_cells[nei_index[i]]->level >= c->level;
      }

      Face* fz = c->external_faces_z[i];
      if ((fz->neigh_cells[0] != NULL) && (fz->neigh_cells[1] != NULL)) {
        flag3[i] = (fz->isleaf == false) ||
                   fz->neigh_cells[nei_index[i]]->level >= c->level;
      }
    }
    flag = flag1[0] && flag1[1] && flag2[0] && flag2[1] && flag3[0] && flag3[1];

    if (flag == false) {
      // refine neibouring cells with lower level
      for (int i = 0; i < 2; i++) {
        // cout<<i<<endl;
        Face* fx = c->external_faces_x[i];
        if (flag1[i] == false) {
          if (fx->neigh_cells[i] != c && fx->neigh_cells[i] != NULL) {
            map<unsigned int, Cell*> temp;
            temp = refinement(fx->neigh_cells[i]);
            split_cells.insert(temp.begin(), temp.end());
          }
          assert(fx->neigh_cells[i]->level >= c->level);
        }

        Face* fy = c->external_faces_y[i];
        if (flag2[i] == false) {
          if (fy->neigh_cells[i] != c && fy->neigh_cells[i] != NULL) {
            map<unsigned int, Cell*> temp;
            temp = refinement(fy->neigh_cells[i]);
            split_cells.insert(temp.begin(), temp.end());
          }
          assert(fy->neigh_cells[i]->level >= c->level);
        }

        Face* fz = c->external_faces_z[i];
        if (flag3[i] == false) {
          if (fz->neigh_cells[i] != c && fz->neigh_cells[i] != NULL) {
            map<unsigned int, Cell*> temp;
            temp = refinement(fz->neigh_cells[i]);
            split_cells.insert(temp.begin(), temp.end());
          }
          assert(fz->neigh_cells[i]->level >= c->level);
        }
      }
    }
  } while (flag == false);

  int child_level = c->level + 1;
  int max_level_previous = cells.size() - 1;
  if (child_level > max_level_previous) {
    cells.push_back(vector<Cell*>(0));
    faces.push_back(vector<Face*>(0));
    num_cell.push_back(0);
    num_face.push_back(0);
  }
  // cout << "CD" << endl;
  if (flag) {
    c->isleaf = false;

    double x[3], y[3], z[3];
    double dx, dy, dz;
    assert(c->_x[0] < c->_x[1]);
    assert(c->_y[0] < c->_y[1]);
    assert(c->_z[0] < c->_z[1]);
    z[0] = c->_z[0];
    x[0] = c->_x[0];
    y[0] = c->_y[0];
    c->get_center(x[1], y[1], z[1]);
    c->get_size(dx, dy, dz);
    z[2] = c->_z[1];
    x[2] = c->_x[1];
    y[2] = c->_y[1];

    Cell* c_child[8];
    int id = 0;
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        for (int k = 0; k < 2; k++) {
          c_child[id] = new Cell(x[i], y[j], z[k], x[i + 1], y[j + 1], z[k + 1],
                                 child_level, n_parameters, true);
          assert(c->parameters.size() == n_parameters);
          assert(c->parameters.size() == c_child[id]->parameters.size());
          for (int ip = 0; ip < n_parameters; ip++) {
            c_child[id]->set_parameter(c->get_parameter(ip), ip);
          }
          id++;
        }
      }
    }

    for (int i = 0; i < 8; i++) {
      c->child_cells[i] = c_child[i];
      cells[child_level].push_back(c_child[i]);
    }

    Face* internal_child_x[4];
    Face* internal_child_y[4];
    Face* internal_child_z[4];

    internal_child_x[0] =
        new Face(c_child[0], c_child[4], x[1], y[1] - 0.25 * dy,
                 z[1] - 0.25 * dz, NORTH_SOUTH, child_level, true);
    internal_child_x[1] =
        new Face(c_child[1], c_child[5], x[1], y[1] - 0.25 * dy,
                 z[1] + 0.25 * dz, NORTH_SOUTH, child_level, true);
    internal_child_x[2] =
        new Face(c_child[2], c_child[6], x[1], y[1] + 0.25 * dy,
                 z[1] - 0.25 * dz, NORTH_SOUTH, child_level, true);
    internal_child_x[3] =
        new Face(c_child[3], c_child[7], x[1], y[1] + 0.25 * dy,
                 z[1] + 0.25 * dz, NORTH_SOUTH, child_level, true);

    c->set_internal_faces(internal_child_x, NORTH_SOUTH);

    for (int i = 0; i < 4; i++) {
      faces[child_level].push_back(internal_child_x[i]);
      leaf_faces.push_back(internal_child_x[i]);
    }

    internal_child_y[0] =
        new Face(c_child[0], c_child[2], x[1] - 0.25 * dx, y[1],
                 z[1] - 0.25 * dz, WEST_EAST, child_level, true);
    internal_child_y[1] =
        new Face(c_child[1], c_child[3], x[1] - 0.25 * dx, y[1],
                 z[1] + 0.25 * dz, WEST_EAST, child_level, true);
    internal_child_y[2] =
        new Face(c_child[4], c_child[6], x[1] + 0.25 * dx, y[1],
                 z[1] - 0.25 * dz, WEST_EAST, child_level, true);
    internal_child_y[3] =
        new Face(c_child[5], c_child[7], x[1] + 0.25 * dx, y[1],
                 z[1] + 0.25 * dz, WEST_EAST, child_level, true);

    c->set_internal_faces(internal_child_y, WEST_EAST);

    for (int i = 0; i < 4; i++) {
      faces[child_level].push_back(internal_child_y[i]);
      leaf_faces.push_back(internal_child_y[i]);
    }

    internal_child_z[0] =
        new Face(c_child[0], c_child[1], x[1] - 0.25 * dx, y[1] - 0.25 * dy,
                 z[1], UP_DOWN, child_level, true);
    internal_child_z[1] =
        new Face(c_child[2], c_child[3], x[1] - 0.25 * dx, y[1] + 0.25 * dy,
                 z[1], UP_DOWN, child_level, true);
    internal_child_z[2] =
        new Face(c_child[4], c_child[5], x[1] + 0.25 * dx, y[1] - 0.25 * dy,
                 z[1], UP_DOWN, child_level, true);
    internal_child_z[3] =
        new Face(c_child[6], c_child[7], x[1] + 0.25 * dx, y[1] + 0.25 * dy,
                 z[1], UP_DOWN, child_level, true);

    c->set_internal_faces(internal_child_z, UP_DOWN);

    for (int i = 0; i < 4; i++) {
      faces[child_level].push_back(internal_child_z[i]);
      leaf_faces.push_back(internal_child_z[i]);
    }
    num_face[child_level] = num_face[child_level] + 12;

    Face* f;
    for (int i = 0; i < 2; i++) {
      f = c->external_faces_x[i];
      Face* f_child[4];
      if (f->isleaf == true) {
        f->isleaf = false;
        if (i == 1 && f->neigh_cells[1] == NULL) {
          f_child[0] =
              new Face(c_child[0 + i * 4], NULL, f->xc, f->yc - 0.25 * dy,
                       f->zc - 0.25 * dz, NORTH_SOUTH, child_level, true);
          f_child[1] =
              new Face(c_child[1 + i * 4], NULL, f->xc, f->yc - 0.25 * dy,
                       f->zc + 0.25 * dz, NORTH_SOUTH, child_level, true);
          f_child[2] =
              new Face(c_child[2 + i * 4], NULL, f->xc, f->yc + 0.25 * dy,
                       f->zc - 0.25 * dz, NORTH_SOUTH, child_level, true);
          f_child[3] =
              new Face(c_child[3 + i * 4], NULL, f->xc, f->yc + 0.25 * dy,
                       f->zc + 0.25 * dz, NORTH_SOUTH, child_level, true);
        } else {
          f_child[0] = new Face(f->neigh_cells[i], c_child[0 + i * 4], f->xc,
                                f->yc - 0.25 * dy, f->zc - 0.25 * dz,
                                NORTH_SOUTH, child_level, true);
          f_child[1] = new Face(f->neigh_cells[i], c_child[1 + i * 4], f->xc,
                                f->yc - 0.25 * dy, f->zc + 0.25 * dz,
                                NORTH_SOUTH, child_level, true);
          f_child[2] = new Face(f->neigh_cells[i], c_child[2 + i * 4], f->xc,
                                f->yc + 0.25 * dy, f->zc - 0.25 * dz,
                                NORTH_SOUTH, child_level, true);
          f_child[3] = new Face(f->neigh_cells[i], c_child[3 + i * 4], f->xc,
                                f->yc + 0.25 * dy, f->zc + 0.25 * dz,
                                NORTH_SOUTH, child_level, true);
        }

        f->set_child_faces(f_child[0], f_child[1], f_child[2], f_child[3]);

        for (int j = 0; j < 4; j++) {
          faces[child_level].push_back(f_child[j]);
        }
        num_face[child_level] = num_face[child_level] + 4;

        vector<Face*>::iterator it_to_be_deleted =
            find(leaf_faces.begin(), leaf_faces.end(), f);
        vector<Face*>::iterator it_insert_point =
            leaf_faces.erase(it_to_be_deleted);
        leaf_faces.insert(it_insert_point,
                          {f_child[0], f_child[1], f_child[2], f_child[3]});
      } else {
        f->child_faces[0]->set_neigh_cells(f->child_faces[0]->neigh_cells[i],
                                           c_child[0 + i * 4], NORTH_SOUTH);
        f->child_faces[1]->set_neigh_cells(f->child_faces[1]->neigh_cells[i],
                                           c_child[1 + i * 4], NORTH_SOUTH);
        f->child_faces[2]->set_neigh_cells(f->child_faces[2]->neigh_cells[i],
                                           c_child[2 + i * 4], NORTH_SOUTH);
        f->child_faces[3]->set_neigh_cells(f->child_faces[3]->neigh_cells[i],
                                           c_child[3 + i * 4], NORTH_SOUTH);
      }
      c_child[0 + i * 4]->set_external_faces(
          f->child_faces[0], c->internal_faces_x[0], NORTH_SOUTH);
      c_child[1 + i * 4]->set_external_faces(
          f->child_faces[1], c->internal_faces_x[1], NORTH_SOUTH);
      c_child[2 + i * 4]->set_external_faces(
          f->child_faces[2], c->internal_faces_x[2], NORTH_SOUTH);
      c_child[3 + i * 4]->set_external_faces(
          f->child_faces[3], c->internal_faces_x[3], NORTH_SOUTH);
    }

    for (int i = 0; i < 2; i++) {
      f = c->external_faces_y[i];
      Face* f_child[4];
      if (f->isleaf) {
        f->isleaf = false;
        if (i == 1 && f->neigh_cells[i] == NULL) {
          f_child[0] =
              new Face(c_child[0 + i * 2], NULL, f->xc - 0.25 * dx, f->yc,
                       f->zc - 0.25 * dz, WEST_EAST, child_level, true);
          f_child[1] =
              new Face(c_child[1 + i * 2], NULL, f->xc - 0.25 * dx, f->yc,
                       f->zc + 0.25 * dz, WEST_EAST, child_level, true);
          f_child[2] =
              new Face(c_child[4 + i * 2], NULL, f->xc + 0.25 * dx, f->yc,
                       f->zc - 0.25 * dz, WEST_EAST, child_level, true);
          f_child[3] =
              new Face(c_child[5 + i * 2], NULL, f->xc + 0.25 * dx, f->yc,
                       f->zc + 0.25 * dz, WEST_EAST, child_level, true);
        } else {
          f_child[0] =
              new Face(f->neigh_cells[i], c_child[0 + i * 2], f->xc - 0.25 * dx,
                       f->yc, f->zc - 0.25 * dz, WEST_EAST, child_level, true);
          f_child[1] =
              new Face(f->neigh_cells[i], c_child[1 + i * 2], f->xc - 0.25 * dx,
                       f->yc, f->zc + 0.25 * dz, WEST_EAST, child_level, true);
          f_child[2] =
              new Face(f->neigh_cells[i], c_child[4 + i * 2], f->xc + 0.25 * dx,
                       f->yc, f->zc - 0.25 * dz, WEST_EAST, child_level, true);
          f_child[3] =
              new Face(f->neigh_cells[i], c_child[5 + i * 2], f->xc + 0.25 * dx,
                       f->yc, f->zc + 0.25 * dz, WEST_EAST, child_level, true);
        }

        f->set_child_faces(f_child[0], f_child[1], f_child[2], f_child[3]);
        for (int j = 0; j < 4; j++) {
          faces[child_level].push_back(f_child[j]);
        }
        num_face[child_level] = num_face[child_level] + 4;

        vector<Face*>::iterator it_to_be_deleted =
            find(leaf_faces.begin(), leaf_faces.end(), f);
        vector<Face*>::iterator it_insert_point =
            leaf_faces.erase(it_to_be_deleted);
        leaf_faces.insert(it_insert_point,
                          {f_child[0], f_child[1], f_child[2], f_child[3]});
      } else {
        f->child_faces[0]->set_neigh_cells(f->child_faces[0]->neigh_cells[i],
                                           c_child[0 + i * 2], WEST_EAST);
        f->child_faces[1]->set_neigh_cells(f->child_faces[1]->neigh_cells[i],
                                           c_child[1 + i * 2], WEST_EAST);
        f->child_faces[2]->set_neigh_cells(f->child_faces[2]->neigh_cells[i],
                                           c_child[4 + i * 2], WEST_EAST);
        f->child_faces[3]->set_neigh_cells(f->child_faces[3]->neigh_cells[i],
                                           c_child[5 + i * 2], WEST_EAST);
      }
      c_child[0 + i * 2]->set_external_faces(f->child_faces[0],
                                             c->internal_faces_y[0], WEST_EAST);
      c_child[1 + i * 2]->set_external_faces(f->child_faces[1],
                                             c->internal_faces_y[1], WEST_EAST);
      c_child[4 + i * 2]->set_external_faces(f->child_faces[2],
                                             c->internal_faces_y[2], WEST_EAST);
      c_child[5 + i * 2]->set_external_faces(f->child_faces[3],
                                             c->internal_faces_y[3], WEST_EAST);
    }

    for (int i = 0; i < 2; i++) {
      f = c->external_faces_z[i];
      Face* f_child[4];
      if (f->isleaf) {
        f->isleaf = false;
        if (i == 1 && f->neigh_cells[i] == NULL) {
          f_child[0] =
              new Face(c_child[0 + i * 1], NULL, f->xc - 0.25 * dx,
                       f->yc - 0.25 * dy, f->zc, UP_DOWN, child_level, true);
          f_child[1] =
              new Face(c_child[2 + i * 1], NULL, f->xc - 0.25 * dx,
                       f->yc + 0.25 * dy, f->zc, UP_DOWN, child_level, true);
          f_child[2] =
              new Face(c_child[4 + i * 1], NULL, f->xc + 0.25 * dx,
                       f->yc - 0.25 * dy, f->zc, UP_DOWN, child_level, true);
          f_child[3] =
              new Face(c_child[6 + i * 1], NULL, f->xc + 0.25 * dx,
                       f->yc + 0.25 * dy, f->zc, UP_DOWN, child_level, true);
        } else {
          if (f->neigh_cells[0] != NULL && f->neigh_cells[1] != NULL) {
            if (!(abs(f->neigh_cells[0]->_z[1] - f->zc) < 1e-10)) {
              f->display();
            }

            assert(abs(f->neigh_cells[0]->_z[1] - f->neigh_cells[1]->_z[0]) <
                   1e-10);
            assert(abs(f->neigh_cells[0]->_z[1] - f->zc) < 1e-10);
          }

          f_child[0] =
              new Face(f->neigh_cells[i], c_child[0 + i * 1], f->xc - 0.25 * dx,
                       f->yc - 0.25 * dy, f->zc, UP_DOWN, child_level, true);
          f_child[1] =
              new Face(f->neigh_cells[i], c_child[2 + i * 1], f->xc - 0.25 * dx,
                       f->yc + 0.25 * dy, f->zc, UP_DOWN, child_level, true);
          f_child[2] =
              new Face(f->neigh_cells[i], c_child[4 + i * 1], f->xc + 0.25 * dx,
                       f->yc - 0.25 * dy, f->zc, UP_DOWN, child_level, true);
          f_child[3] =
              new Face(f->neigh_cells[i], c_child[6 + i * 1], f->xc + 0.25 * dx,
                       f->yc + 0.25 * dy, f->zc, UP_DOWN, child_level, true);
        }

        f->set_child_faces(f_child[0], f_child[1], f_child[2], f_child[3]);
        for (int j = 0; j < 4; j++) {
          faces[child_level].push_back(f_child[j]);
        }
        num_face[child_level] = num_face[child_level] + 4;

        vector<Face*>::iterator it_to_be_deleted =
            find(leaf_faces.begin(), leaf_faces.end(), f);
        vector<Face*>::iterator it_insert_point =
            leaf_faces.erase(it_to_be_deleted);
        leaf_faces.insert(it_insert_point,
                          {f_child[0], f_child[1], f_child[2], f_child[3]});
      } else {
        f->child_faces[0]->set_neigh_cells(f->child_faces[0]->neigh_cells[i],
                                           c_child[0 + i * 1], UP_DOWN);
        f->child_faces[1]->set_neigh_cells(f->child_faces[1]->neigh_cells[i],
                                           c_child[2 + i * 1], UP_DOWN);
        f->child_faces[2]->set_neigh_cells(f->child_faces[2]->neigh_cells[i],
                                           c_child[4 + i * 1], UP_DOWN);
        f->child_faces[3]->set_neigh_cells(f->child_faces[3]->neigh_cells[i],
                                           c_child[6 + i * 1], UP_DOWN);
      }
      c_child[0 + i * 1]->set_external_faces(f->child_faces[0],
                                             c->internal_faces_z[0], UP_DOWN);
      c_child[2 + i * 1]->set_external_faces(f->child_faces[1],
                                             c->internal_faces_z[1], UP_DOWN);
      c_child[4 + i * 1]->set_external_faces(f->child_faces[2],
                                             c->internal_faces_z[2], UP_DOWN);
      c_child[6 + i * 1]->set_external_faces(f->child_faces[3],
                                             c->internal_faces_z[3], UP_DOWN);
    }

    vector<Cell*>::iterator c_iter =
        std::find(leaf_cells.begin(), leaf_cells.end(), c);
    // assert(c_iter != leaf_cells.end());
    vector<Cell*>::iterator iter = leaf_cells.erase(c_iter);
    //vector<Cell*>::iterator iter2 = 
    leaf_cells.insert(
        iter, {c_child[0], c_child[1], c_child[2], c_child[3], c_child[4],
               c_child[5], c_child[6], c_child[7]});
    num_leaf_cells = leaf_cells.size();
    num_leaf_faces = leaf_faces.size();

    num_cell[child_level] = num_cell[child_level] + 8;

    // assert(iter2 != leaf_cells.end());
    // split_cells[c->id]=iter2;
    split_cells.insert(pair<unsigned int, Cell*>(c->id, c));
    // cout << "id=" << c->id << ", " << (*split_cells[c->id])->isleaf << endl;
    // return iter2;
  }

  return split_cells;
  // else
  // {
  //     cout << "leaf_cells.end()" << endl;
  //     return leaf_cells.end();
  // }
}

map<unsigned int, Cell*> Mesh::refinement(int index_to_be_refined) {
  Cell* c = leaf_cells[index_to_be_refined];
  map<unsigned int, Cell*> map_cells = this->refinement(c);
  return map_cells;
}

void Mesh::get_model_parameter_from_mesh(VectorXd& m, int ith) {
  m.resize(leaf_cells.size());
  for (int i = 0; i < leaf_cells.size(); i++) {
    m(i) = leaf_cells[i]->get_parameter(ith);
  }
}

void Mesh::rearrange_id() {
  for (int i = 0; i < leaf_cells.size(); i++) {
    leaf_cells[i]->set_id(i);
  }
}

bool Mesh::great_equal(long double left, long double right) {
  bool temp = false;
  // left>=right?
  long double a = left - right;
  if (a > 0.)
    temp = true;  // a>0, OK.
  if (std::abs(a) < TOL)
    temp = true;  // a=0, OK.
  return temp;
}

void Mesh::sort(int level) {
  for (int i = 0; i < cells[level].size(); i++) {
    Cell* c = cells[level][i];
    // sort faces
    if (c->external_faces_z[0]->zc > c->external_faces_z[1]->zc) {
      Face* t = c->external_faces_z[0];
      c->external_faces_z[0] = c->external_faces_z[1];
      c->external_faces_z[1] = t;
    }

    if (c->external_faces_x[0]->xc > c->external_faces_x[1]->xc) {
      Face* t = c->external_faces_x[0];
      c->external_faces_x[0] = c->external_faces_x[1];
      c->external_faces_x[1] = t;
    }

    if (c->external_faces_y[0]->yc > c->external_faces_y[1]->yc) {
      Face* t = c->external_faces_y[0];
      c->external_faces_y[0] = c->external_faces_y[1];
      c->external_faces_y[1] = t;
    }

    // sort neighbouring cells
    for (int j = 0; j < 2; j++) {
      Cell *n1, *n2;
      n1 = c->external_faces_z[j]->neigh_cells[0];
      n2 = c->external_faces_z[j]->neigh_cells[1];

      double zc1, xc1, yc1;
      double zc2, xc2, yc2;

      if (n1 != NULL && n2 != NULL) {
        n1->get_center(xc1, yc1, zc1);
        n2->get_center(zc2, xc2, yc2);
        if (zc1 > zc2) {
          c->external_faces_z[j]->neigh_cells[0] = n2;
          c->external_faces_z[j]->neigh_cells[1] = n1;
        }
      } else {
        double rc0, thetac0, phic0;
        c->get_center(rc0, thetac0, phic0);
        if (rc0 > (c->external_faces_z[j]->zc)) {
          if (n2 == NULL && n1 == c) {
            c->external_faces_z[j]->neigh_cells[0] = NULL;
            c->external_faces_z[j]->neigh_cells[1] = c;
          }
        } else {
          if (n1 == NULL && n2 == c) {
            c->external_faces_z[j]->neigh_cells[0] = c;
            c->external_faces_z[j]->neigh_cells[1] = NULL;
          }
        }
      }
    }

    for (int j = 0; j < 2; j++) {
      Cell *n1, *n2;
      n1 = c->external_faces_x[j]->neigh_cells[0];
      n2 = c->external_faces_x[j]->neigh_cells[1];

      double zc1, xc1, yc1;
      double zc2, xc2, yc2;

      if (n1 != NULL && n2 != NULL) {
        n1->get_center(xc1, yc1, zc1);
        n2->get_center(zc2, xc2, yc2);
        if (xc1 > xc2) {
          c->external_faces_x[j]->neigh_cells[0] = n2;
          c->external_faces_x[j]->neigh_cells[1] = n1;
        }
      } else {
        double rc0, thetac0, phic0;
        c->get_center(rc0, thetac0, phic0);
        if (thetac0 > (c->external_faces_x[j]->xc)) {
          if (n2 == NULL && n1 == c) {
            c->external_faces_x[j]->neigh_cells[0] = NULL;
            c->external_faces_x[j]->neigh_cells[1] = c;
          }
        } else {
          if (n1 == NULL && n2 == c) {
            c->external_faces_x[j]->neigh_cells[0] = c;
            c->external_faces_x[j]->neigh_cells[1] = NULL;
          }
        }
      }
    }

    for (int j = 0; j < 2; j++) {
      Cell *n1, *n2;
      n1 = c->external_faces_y[j]->neigh_cells[0];
      n2 = c->external_faces_y[j]->neigh_cells[1];

      double zc1, xc1, yc1;
      double zc2, xc2, yc2;

      if (n1 != NULL && n2 != NULL) {
        n1->get_center(xc1, yc1, zc1);
        n2->get_center(zc2, xc2, yc2);
        if (yc1 > yc2) {
          c->external_faces_y[j]->neigh_cells[0] = n2;
          c->external_faces_y[j]->neigh_cells[1] = n1;
        }
      } else {
        double rc0, thetac0, phic0;
        c->get_center(rc0, thetac0, phic0);
        if (phic0 > (c->external_faces_y[j]->yc)) {
          if (n2 == NULL && n1 == c) {
            c->external_faces_y[j]->neigh_cells[0] = NULL;
            c->external_faces_y[j]->neigh_cells[1] = c;
          }
        } else {
          if (n1 == NULL && n2 == c) {
            c->external_faces_y[j]->neigh_cells[0] = c;
            c->external_faces_y[j]->neigh_cells[1] = NULL;
          }
        }
      }
    }
  }
}

void Mesh::fill_data(int offset_i,
                     int offset_j,
                     int offset_k,
                     double*** data,
                     Cell* c,
                     int max_level,
                     int ith_para) {
  if (c->isleaf) {
    int level_difference = max_level - c->level;
    int N = pow(2, level_difference);
    double value = c->get_parameter(ith_para);
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        for (int k = 0; k < N; k++) {
          data[offset_i + i][offset_j + j][offset_k + k] = value;
        }
      }
    }
  } else {
    int level_difference = max_level - c->child_cells[0]->level;
    int N = pow(2, level_difference);
    //递归
    fill_data(offset_i, offset_j, offset_k, data, c->child_cells[0], max_level,
              ith_para);
    fill_data(offset_i, offset_j, offset_k + N, data, c->child_cells[1],
              max_level, ith_para);
    fill_data(offset_i, offset_j + N, offset_k, data, c->child_cells[2],
              max_level, ith_para);
    fill_data(offset_i, offset_j + N, offset_k + N, data, c->child_cells[3],
              max_level, ith_para);
    fill_data(offset_i + N, offset_j, offset_k, data, c->child_cells[4],
              max_level, ith_para);
    fill_data(offset_i + N, offset_j, offset_k + N, data, c->child_cells[5],
              max_level, ith_para);
    fill_data(offset_i + N, offset_j + N, offset_k, data, c->child_cells[6],
              max_level, ith_para);
    fill_data(offset_i + N, offset_j + N, offset_k + N, data, c->child_cells[7],
              max_level, ith_para);
  }
}
#ifdef USE_NETCDF
int Mesh::out_model_netcdf(string filename,
                           int ith_para,
                           string VAL_NAME,
                           string VAL_UNITS) {
  int max_level = cells.size() - 1;
  int N = std::pow(2, max_level);
  int NZ = nz * N;
  int NX = nx * N;
  int NY = ny * N;

  assert(cells[0].size() == nz * nx * ny);

  int n_total = std::pow(8, max_level) * nz * nx * ny;
  double*** DENSITY_DATA;
  DENSITY_DATA = new double**[NX];
  for (unsigned int i = 0; i < NX; i++) {
    DENSITY_DATA[i] = new double*[NY];
    for (unsigned int j = 0; j < NY; j++) {
      DENSITY_DATA[i][j] = new double[NZ];
    }
  }

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        fill_data(i * N, j * N, k * N, DENSITY_DATA,
                  cells[0][i * ny * nz + j * nz + k], max_level, ith_para);
      }
    }
  }

  // writh data to netcdf file
  double* xs = new double[NX];
  double* ys = new double[NY];
  double* zs = new double[NZ];

  double** xs_bnd = new double*[NX];
  double** ys_bnd = new double*[NY];
  double** z_bnd = new double*[NZ];

  double x_space = (x_lim[1] - x_lim[0]) / NX;
  for (unsigned int i = 0; i < NX; i++) {
    xs[i] = x_lim[0] + 0.5 * x_space + i * x_space;
    xs_bnd[i] = new double[2];
    xs_bnd[i][0] = x_lim[0] + i * x_space;
    xs_bnd[i][1] = x_lim[0] + (i + 1) * x_space;
  }

  double y_space = (y_lim[1] - y_lim[0]) / NY;
  for (unsigned int j = 0; j < NY; j++) {
    ys[j] = y_lim[0] + 0.5 * y_space + j * y_space;
    ys_bnd[j] = new double[2];

    ys_bnd[j][0] = y_lim[0] + j * y_space;
    ys_bnd[j][1] = y_lim[0] + (j + 1) * y_space;
  }
  for (unsigned int k = 0; k < nz; k++) {
    double z_space = (z_points(k + 1) - z_points(k)) / (1.0 * N);
    for (unsigned int k2 = 0; k2 < N; k2++) {
      int index = k * N + k2;
      zs[index] = z_points(k) + 0.5 * z_space + k2 * z_space;
      z_bnd[index] = new double[2];
      z_bnd[index][0] = z_points(k) + k2 * z_space;
      z_bnd[index][1] = z_points(k) + (k2 + 1) * z_space;
    }
  }
  // Names of things.
  const char* X_NAME = "x";
  const char* Y_NAME = "y";
  const char* VAL_NAME_C = VAL_NAME.c_str();
  const char* Z_NAME = "z";

  string UNITS = "units";
  string AXIS = "axis";
  // string POSITIVE = "positive";
  string DEGREES_EAST = "degrees_east";
  string DEGREES_NORTH = "degrees_north";
  string METER = "m";
  // string UP = "up";
  string DOWN = "down";
  // For the units attributes.

  string BOUNDS = "bounds";
  string XBND_NAME = "xbnd";
  string YBND_NAME = "ybnd";
  string ZBND_NAME = "zbnd";

  try {
    // Create the file. The Replace parameter tells netCDF to overwrite
    // this file, if it already exists.
    NcFile test(filename, NcFile::replace);
    test.putAtt("Conventions", "CF-1.7");
    test.putAtt("node_offset", ncInt, 1);
    // test.putAtt(NODE_OFFSET, "1");
    // test.putAtt(NODE_OFFSET, ncInt, 1, pixel_registration);
    // Define the dimensions. NetCDF will hand back an ncDim object for
    // each.
    NcDim zDim = test.addDim(Z_NAME, NZ);
    NcDim yDim = test.addDim(Y_NAME, NY);
    NcDim xDim = test.addDim(X_NAME, NX);

    NcVar zVar = test.addVar(Z_NAME, ncDouble, zDim);
    NcVar yVar = test.addVar(Y_NAME, ncDouble, yDim);
    NcVar xVar = test.addVar(X_NAME, ncDouble, xDim);
    // Define units attributes for coordinate vars. This attaches a
    // text attribute to each of the coordinate variables, containing
    // the units.
    // xVar.putAtt("long_name", "x-coordinate in Cartesian system");

    // xVar.putAtt("standard_name", "projection_x_coordinate");
    // xVar.putAtt(UNITS, METER);
    // xVar.putAtt(UNITS, METER);
    // xVar.putAtt(AXIS, "X");

    // yVar.putAtt(UNITS, METER);
    // yVar.putAtt(AXIS, "Y");
    double zrange[2] = {z_lim[0], z_lim[1]};
    zVar.putAtt("long_name", "z");
    zVar.putAtt("actual_range", ncDouble, 2, zrange);

    double yrange[2] = {y_lim[0], y_lim[1]};
    yVar.putAtt("long_name", "y");
    yVar.putAtt("actual_range", ncDouble, 2, yrange);

    double xrange[2] = {x_lim[0], x_lim[1]};
    xVar.putAtt("long_name", "x");
    xVar.putAtt("actual_range", ncDouble, 2, xrange);
    // yVar.putAtt("long_name", "y-coordinate in Cartesian system");
    // yVar.putAtt("standard_name", "projection_y_coordinate");

    // zVar.putAtt(UNITS, METER);
    // zVar.putAtt(POSITIVE, DOWN);

    // Write the coordinate variable data to the file.
    zVar.putVar(zs);
    yVar.putVar(ys);
    xVar.putVar(xs);

    // bounds for cells
    NcDim bndDim = test.addDim("bnd", 2);

    vector<NcDim> dimVector_zbnd;
    dimVector_zbnd.push_back(zDim);
    dimVector_zbnd.push_back(bndDim);
    NcVar zbnd = test.addVar(ZBND_NAME, ncDouble, dimVector_zbnd);
    zVar.putAtt(BOUNDS, ZBND_NAME);

    vector<NcDim> dimVector_ybnd;
    dimVector_ybnd.push_back(yDim);
    dimVector_ybnd.push_back(bndDim);
    NcVar ybnd = test.addVar(YBND_NAME, ncDouble, dimVector_ybnd);
    yVar.putAtt(BOUNDS, YBND_NAME);

    vector<NcDim> dimVector_xbnd;
    dimVector_xbnd.push_back(xDim);
    dimVector_xbnd.push_back(bndDim);
    NcVar xbnd = test.addVar(XBND_NAME, ncDouble, dimVector_xbnd);
    xVar.putAtt(BOUNDS, XBND_NAME);

    vector<size_t> startp_bnd, countp_bnd;
    startp_bnd.push_back(0);
    startp_bnd.push_back(0);
    countp_bnd.push_back(1);
    countp_bnd.push_back(1);
    for (int i = 0; i < NZ; i++) {
      startp_bnd[0] = i;
      for (int j = 0; j < 2; j++) {
        startp_bnd[1] = j;
        double a = z_bnd[i][j];
        zbnd.putVar(startp_bnd, countp_bnd, &a);
      }
    }
    for (int i = 0; i < NY; i++) {
      startp_bnd[0] = i;
      for (int j = 0; j < 2; j++) {
        startp_bnd[1] = j;
        double a = ys_bnd[i][j];
        ybnd.putVar(startp_bnd, countp_bnd, &a);
      }
    }
    for (int i = 0; i < NX; i++) {
      startp_bnd[0] = i;
      for (int j = 0; j < 2; j++) {
        startp_bnd[1] = j;
        double a = xs_bnd[i][j];
        xbnd.putVar(startp_bnd, countp_bnd, &a);
      }
    }

    // Define the netCDF variables for the density data
    vector<NcDim> dimVector;
    // dimVector.push_back(xDim);
    // dimVector.push_back(yDim);
    // dimVector.push_back(zDim);

    dimVector.push_back(zDim);
    dimVector.push_back(yDim);
    dimVector.push_back(xDim);
    NcVar denVar = test.addVar(VAL_NAME_C, ncDouble, dimVector);

    denVar.putAtt(UNITS, VAL_UNITS);
    // denVar.putAtt(NODE_OFFSET, "1");

    // Write the density data;
    vector<size_t> startp, countp;
    startp.push_back(0);
    startp.push_back(0);
    startp.push_back(0);

    countp.push_back(1);
    countp.push_back(1);
    countp.push_back(1);

    // for (size_t i = 0; i < NX; i++) {
    //   for (size_t j = 0; j < NY; j++) {
    //     for (size_t k = 0; k < NZ; k++) {
    // cout << "here" << endl;
    for (size_t i = 0; i < NZ; i++) {
      for (size_t j = 0; j < NY; j++) {
        for (size_t k = 0; k < NX; k++) {
          startp[0] = i;
          startp[1] = j;
          startp[2] = k;
          // double a = DENSITY_DATA[i][j][k];
          double a = DENSITY_DATA[k][j][i];

          // startp[0] = k;
          // startp[1] = j;
          // startp[2] = i;
          // double a = DENSITY_DATA[k][j][i];
          denVar.putVar(startp, countp, &a);
        }
      }
    }

    // free resources
    for (int i = 0; i < NX; i++) {
      for (int j = 0; j < NY; j++) {
        delete[] DENSITY_DATA[i][j];
        DENSITY_DATA[i][j] = NULL;
      }
      delete[] DENSITY_DATA[i];
      DENSITY_DATA[i] = NULL;
    }
    delete[] DENSITY_DATA;
    DENSITY_DATA = NULL;

    delete[] xs;
    xs = NULL;

    delete[] ys;
    ys = NULL;

    delete[] zs;
    zs = NULL;

    for (int i = 0; i < NX; i++) {
      delete[] xs_bnd[i];
      xs_bnd[i] = NULL;
    }
    delete[] xs_bnd;
    xs_bnd = NULL;
    for (int j = 0; j < NY; j++) {
      delete[] ys_bnd[j];
      ys_bnd[j] = NULL;
    }
    delete[] ys_bnd;
    ys_bnd = NULL;
    for (int k = 0; k < NZ; k++) {
      delete[] z_bnd[k];
      z_bnd[k] = NULL;
    }
    delete[] z_bnd;
    z_bnd = NULL;

    // denVar.putVar(DENSITY_DATA);
    cout << "The model has been written to NetCDF file: " << filename << endl;
    return 0;
  } catch (NcException& e) {
    e.what();
    return NC_ERR;
  }
}
#endif
