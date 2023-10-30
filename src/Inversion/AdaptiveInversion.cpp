#include "AdaptiveInversion.h"
AdaptiveInversion::AdaptiveInversion(const Mesh &mesh_, const Observation &ob_,
                                     unsigned long long field_flag_)
    : InversionBase(mesh_, ob_, field_flag_),
      refinement_percentage(0.1),
      max_refinement_number(0),
      interval_between_refinements(1),
      min_dx(0),
      min_dy(0),
      min_dz(0) {}

AdaptiveInversion::AdaptiveInversion()
    : InversionBase(),
      refinement_percentage(0.1),
      max_refinement_number(0),
      interval_between_refinements(1),
      min_dx(0),
      min_dy(0),
      min_dz(0) {}

AdaptiveInversion::~AdaptiveInversion() {}

void AdaptiveInversion::expand_G(const map<unsigned int, Cell *> &split_cells)
{
    int num_newcells = split_cells.size();
    // cout << "num_newcells=" << num_newcells << endl;
    MatrixXd G2 = this->G;
    this->G.resize(G2.rows(), G2.cols() + 7 * num_newcells);
    MatrixXd temp(G2.rows(), 8);

    // bitset<10> basic_field_flag;
    // bitset<10> computed_field_flag;
    // //目前只定义了10个基本分量，basic_field_flag就是computed_field_flag
    // //以后扩展组合分量(比如T_xx-T_yy)的时候，computed_field_flag可以包含basic_field_flag
    // //（未实现）
    // for (int i = 0; i < 10; i++) {
    //   basic_field_flag[i] = field_flag[i];
    //   computed_field_flag[i] = field_flag[i];
    // }

    // int n_basic_fields = basic_field_flag.count();
    int n_fields = field_flag.count();

    int rows_number = n_fields * N_obs;
    assert(G2.rows() == n_fields * N_obs);
    assert(Nm > 0);
    assert(N_obs > 0);

    vector<int> basic_field_index;
    for (int i = 0; i < 10; i++)
    {
        // bool status = flag[i];
        if (field_flag[i])
        {
            basic_field_index.push_back(i);
        }
    }

    assert(basic_field_index.size() == n_fields);
    map<unsigned int, Cell *>::const_iterator map_it = split_cells.begin();
    G.block(0, 0, rows_number, map_it->first) =
        G2.block(0, 0, rows_number, map_it->first);
    Cell *parent_cell;
    for (int i_cell = 0; i_cell < num_newcells; i_cell++)
    {
        int index = map_it->first;
        Cell *parent_cell = map_it->second;
        for (int j = 0; j < 8; j++)
        {
            assert(parent_cell->child_cells[j]->isleaf);
        }
        for (int i = 0; i < N_obs; i++)
        {
#pragma omp parallel for
            for (int j = 0; j < 8; j++)
            {
                GravFormula gra;
                vector<double> field;

                // gra.field_for_a_tesseroid(ob(i),
                // parent_cell->child_cells[j], 1.0,field, basic_field_flag);
                gra.field_caused_by_single_prism(ob(i),
                                                 *(parent_cell->child_cells[j]),
                                                 1.0, field, field_flag);

                for (int k = 0; k < n_fields; k++)
                {
                    temp(i + k * N_obs, j) = field[basic_field_index[k]];
                }
            }
        }
        G.block(0, index + i_cell * 7, rows_number, 8) = temp;
        map_it++;
        if (map_it != split_cells.end())
        {
            G.block(0, index + i_cell * 7 + 8, rows_number,
                    map_it->first - index - 1) =
                G2.block(0, index + 1, rows_number, map_it->first - index - 1);
        }
    }
    map_it--;
    int index = map_it->first;
    G.block(0, index + (num_newcells - 1) * 7 + 8, rows_number,
            G2.cols() - 1 - index) =
        G2.block(0, index + 1, rows_number, G2.cols() - 1 - index);
}

void AdaptiveInversion::indicator_calculator(VectorXd &indicator)
{
    indicator.resize(Nm);
    indicator.setZero();
    double dx, dy, dz;

    for (int i = 0; i < Nm; i++)
    {
        Cell *c = this->mesh.leaf_cells[i];
        double mc = m(c->get_id());
        c->get_size(dx, dy, dz);
        for (int j = 0; j < 2; j++)
        {
            Face *fx = c->external_faces_x[j];
            if (fx->isleaf)
            {
                if (fx->neigh_cells[0] != NULL && fx->neigh_cells[1] != NULL)
                {
                    Cell *neigh = fx->neigh_cells[j];
                    int id = neigh->get_id();
                    double area = (c->_z[1] - c->_z[0]) * (c->_y[1] - c->_y[0]);
                    indicator(i) += (m(id) - mc) * (m(id) - mc) * area;
                    // double vol = neigh->get_volumn();
                    // indicator(i) += (m(id) - mc) * (m(id) - mc) * vol;
                }
            }
            else
            {
                for (int k = 0; k < 4; k++)
                {
                    Face *f = fx->child_faces[k];
                    Cell *neigh = f->neigh_cells[j];
                    int id = neigh->get_id();
                    double area = (neigh->_z[1] - neigh->_z[0]) *
                                  (neigh->_y[1] - neigh->_y[0]);
                    indicator(i) += (m(id) - mc) * (m(id) - mc) * area;
                    // double vol = neigh->get_volumn();
                    // indicator(i) += (m(id) - mc) * (m(id) - mc) * vol;
                }
            }

            Face *fy = c->external_faces_y[j];
            if (fy->isleaf)
            {
                if (fy->neigh_cells[0] != NULL && fy->neigh_cells[1] != NULL)
                {
                    Cell *neigh = fy->neigh_cells[j];
                    int id = neigh->get_id();
                    double area = (c->_z[1] - c->_z[0]) * (c->_x[1] - c->_x[0]);
                    indicator(i) += (m(id) - mc) * (m(id) - mc) * area;
                    // double vol = neigh->get_volumn();
                    // indicator(i) += (m(id) - mc) * (m(id) - mc) * vol;
                }
            }
            else
            {
                for (int k = 0; k < 4; k++)
                {
                    Face *f = fy->child_faces[k];
                    Cell *neigh = f->neigh_cells[j];
                    int id = neigh->get_id();
                    double area = (neigh->_z[1] - neigh->_z[0]) *
                                  (neigh->_x[1] - neigh->_x[0]);
                    indicator(i) += (m(id) - mc) * (m(id) - mc) * area;
                    // double vol = neigh->get_volumn();
                    // indicator(i) += (m(id) - mc) * (m(id) - mc) * vol;
                }
            }

            Face *fz = c->external_faces_z[j];
            if (fz->isleaf)
            {
                if (fz->neigh_cells[0] != NULL && fz->neigh_cells[1] != NULL)
                {
                    Cell *neigh = fz->neigh_cells[j];
                    int id = neigh->get_id();
                    double area = (c->_x[1] - c->_x[0]) * (c->_y[1] - c->_y[0]);
                    indicator(i) += (m(id) - mc) * (m(id) - mc) * area;
                    // double vol = neigh->get_volumn();
                    // indicator(i) += (m(id) - mc) * (m(id) - mc) * vol;
                }
            }
            else
            {
                for (int k = 0; k < 4; k++)
                {
                    Face *f = fz->child_faces[k];
                    Cell *neigh = f->neigh_cells[j];
                    int id = neigh->get_id();
                    double area = (neigh->_x[1] - neigh->_x[0]) *
                                  (neigh->_y[1] - neigh->_y[0]);
                    indicator(i) += (m(id) - mc) * (m(id) - mc) * area;
                }
            }
        }
    }
}

void AdaptiveInversion::refine_mesh(double a,
                                    InterpMultilinear<3, double> &interp,
                                    string flag)
{
    this->set_density_to_mesh();
    this->set_reference_model_to_mesh();
    this->set_min_max_to_mesh();

    VectorXd indicator(Nm);
    assert(Nm == G.cols());
    assert(Nm == m.rows());

    // indicator = ((D_theta1 * m).cwiseAbs2() + (D_phi1 * m).cwiseAbs2() +
    // (D_r1
    // * m).cwiseAbs2()).cwiseSqrt();
    this->indicator_calculator(indicator);
    double max_idt = indicator.maxCoeff();
    double min_idt = indicator.minCoeff();

    // indicator =(indicator - VectorXd::Constant(Nm, min_idt)) / (max_idt -
    // min_idt);

    VectorXd sorted_indicator;
    VectorXi sorted_index;
    sort_vec(indicator, sorted_indicator, sorted_index);

    vector<Cell *> cells_marked(0);
    cells_marked.clear();
    assert(Nm == mesh.leaf_cells.size());

    double tol = sorted_indicator(int(a * Nm));
    for (int i = 0; i < Nm; i++)
    {
        if (indicator(i) > tol)
        {
            Cell *c = this->mesh.leaf_cells[i];
            double c_dx;
            double c_dy;
            double c_dz;
            c->get_size(c_dx, c_dy, c_dz);
            if (c_dx < min_dx || std::abs(c_dx - min_dx) < 1e-7 ||
                c_dy < min_dy || std::abs(c_dy - min_dy) < 1e-7 ||
                c_dz < min_dz || std::abs(c_dz - min_dz) < 1e-7)
            {
                // If the size of the cell is less than the limit of minimum
                // size, do nothing
            }
            else
            {
                cells_marked.push_back(c);
            }
        }
    }

    for (int i = 1; i < cells_marked.size(); i++)
    {
        assert(cells_marked[i - 1]->id < cells_marked[i]->id);
    }

    int counter = 0;
    map<unsigned int, Cell *> split_cells;
    for (int i = 0; i < cells_marked.size(); i++)
    {
        if (cells_marked[i]->isleaf)
        {
            map<unsigned int, Cell *> split_cells0 =
                mesh.refinement(cells_marked[i]);
            split_cells.insert(split_cells0.begin(), split_cells0.end());
        }
    }

    if (split_cells.size() > 0)
    {
        expand_G(split_cells);

        map<unsigned int, Cell *>::const_iterator map_it = split_cells.begin();
        Cell *parent_cell;
        for (int i = 0; i < split_cells.size(); i++)
        {
            int index = map_it->first;
            parent_cell = map_it->second;
            for (int j = 0; j < 8; j++)
            {
                assert(parent_cell->child_cells[j]->isleaf);
                double xc, yc, zc;
                parent_cell->child_cells[j]->get_center(xc, yc, zc);
                array<double, 3> args = {xc, yc, zc};

                if (flag == "crg")
                {
                    parent_cell->child_cells[j]->set_parameter(
                        interp.interp(args.begin()), 2);
                }
                else if (flag == "pet")
                {
                    parent_cell->child_cells[j]->set_parameter(
                        interp.interp(args.begin()), 1);
                }
                else if (flag == "both")
                {
                    parent_cell->child_cells[j]->set_parameter(
                        interp.interp(args.begin()), 2);
                    parent_cell->child_cells[j]->set_parameter(
                        interp.interp(args.begin()), 1);
                }
                else
                {
                    cout << "flag must be \"crg\" or \"pet\" or\"both\""
                         << endl;
                    abort();
                }
            }
            map_it++;
        }
    }
    // cout << G.rows() << ", " << G.cols() << endl;
    mesh.get_model_parameter_from_mesh(m, 0);
    mesh.get_model_parameter_from_mesh(m0, 1);
    mesh.get_model_parameter_from_mesh(m0_s, 2);
    mesh.get_model_parameter_from_mesh(m_min, 3);
    mesh.get_model_parameter_from_mesh(m_max, 4);
    mesh.rearrange_id();
    this->init_matrices();
    // cout << G.rows() << ", " << G.cols() << endl;
    // this->compute_G();
}
void AdaptiveInversion::refine_mesh(
    double a, InterpMultilinear<3, double> &interp_ML_crg,
    InterpMultilinear<3, double> &interp_ML_m0)
{
    this->set_density_to_mesh();
    this->set_reference_model_to_mesh();
    this->set_min_max_to_mesh();

    VectorXd indicator(Nm);
    assert(Nm == G.cols());
    assert(Nm == m.rows());

    // indicator = ((D_theta1 * m).cwiseAbs2() + (D_phi1 * m).cwiseAbs2() +
    // (D_r1
    // * m).cwiseAbs2()).cwiseSqrt();
    this->indicator_calculator(indicator);
    double max_idt = indicator.maxCoeff();
    double min_idt = indicator.minCoeff();

    // indicator =
    //     (indicator - VectorXd::Constant(Nm, min_idt)) / (max_idt - min_idt);

    VectorXd sorted_indicator;
    VectorXi sorted_index;
    sort_vec(indicator, sorted_indicator, sorted_index);

    vector<Cell *> cells_marked(0);
    cells_marked.clear();
    assert(Nm == mesh.leaf_cells.size());

    double tol = sorted_indicator(int(a * Nm));
    for (int i = 0; i < Nm; i++)
    {
        if (indicator(i) > tol || abs(indicator(i) - tol) < 1e-6)
        {
            Cell *c = this->mesh.leaf_cells[i];
            double c_dx;
            double c_dy;
            double c_dz;
            c->get_size(c_dx, c_dy, c_dz);
            if (c_dx < min_dx || std::abs(c_dx - min_dx) < 1e-7 ||
                c_dy < min_dy || std::abs(c_dy - min_dy) < 1e-7 ||
                c_dz < min_dz || std::abs(c_dz - min_dz) < 1e-7)
            {
                // If the size of the cell is less than the limit of minimum
                // size, do nothing
            }
            else
            {
                cells_marked.push_back(c);
            }
        }
    }

    for (int i = 1; i < cells_marked.size(); i++)
    {
        assert(cells_marked[i - 1]->id < cells_marked[i]->id);
    }

    int counter = 0;
    map<unsigned int, Cell *> split_cells;
    for (int i = 0; i < cells_marked.size(); i++)
    {
        if (cells_marked[i]->isleaf)
        {
            map<unsigned int, Cell *> split_cells0 =
                mesh.refinement(cells_marked[i]);
            split_cells.insert(split_cells0.begin(), split_cells0.end());
        }
    }

    if (split_cells.size() > 0)
    {
        expand_G(split_cells);

        map<unsigned int, Cell *>::const_iterator map_it = split_cells.begin();
        Cell *parent_cell;
        for (int i = 0; i < split_cells.size(); i++)
        {
            int index = map_it->first;
            parent_cell = map_it->second;
            for (int j = 0; j < 8; j++)
            {
                assert(parent_cell->child_cells[j]->isleaf);
                double xc, yc, zc;
                parent_cell->child_cells[j]->get_center(xc, yc, zc);

                array<double, 3> args = {xc, yc, zc};

                parent_cell->child_cells[j]->set_parameter(
                    interp_ML_m0.interp(args.begin()), 1);
                parent_cell->child_cells[j]->set_parameter(
                    interp_ML_crg.interp(args.begin()), 2);
            }
            map_it++;
        }
    }
    // cout << G.rows() << ", " << G.cols() << endl;
    mesh.get_model_parameter_from_mesh(m, 0);
    mesh.get_model_parameter_from_mesh(m0, 1);
    mesh.get_model_parameter_from_mesh(m0_s, 2);
    mesh.get_model_parameter_from_mesh(m_min, 3);
    mesh.get_model_parameter_from_mesh(m_max, 4);
    mesh.rearrange_id();
    this->init_matrices();

    // cout << G.rows() << ", " << G.cols() << endl;
    // this->compute_G();
}

// void AdaptiveInversion::refine_mesh(double a, InterpMultilinear<3, double>
// &interp_ML_crg, InterpMultilinear<3, double> &interp_ML_m0, double
// reference_surface)
// {
//     this->set_density_to_mesh();
//     this->set_reference_model_to_mesh();
//     this->set_min_max_to_mesh();

//     VectorXd indicator(Nm);
//     assert(Nm == G.cols());
//     assert(Nm == m.rows());

//     // indicator = ((D_theta1 * m).cwiseAbs2() + (D_phi1 * m).cwiseAbs2() +
//     (D_r1 * m).cwiseAbs2()).cwiseSqrt();
//     this->indicator_calculator(indicator);
//     double max_idt = indicator.maxCoeff();
//     double min_idt = indicator.minCoeff();

//     indicator = (indicator - VectorXd::Constant(Nm, min_idt)) / (max_idt -
//     min_idt);

//     VectorXd sorted_indicator;
//     VectorXi sorted_index;
//     sort_vec(indicator, sorted_indicator, sorted_index);

//     vector<Cell *> cells_marked(0);
//     cells_marked.clear();
//     assert(Nm == mesh.leaf_cells.size());

//     double tol = sorted_indicator(int(a * Nm));
//     for (int i = 0; i < Nm; i++)
//     {
//         if (indicator(i) > tol)
//         {
//             Cell *c = this->mesh.leaf_cells[i];
//             cells_marked.push_back(c);
//         }
//         // if (indicator(i) > a)
//         // {
//         //     Cell *c = this->mesh.leaf_cells[i];
//         //     cells_marked.push_back(c);
//         // }
//     }

//     for (int i = 1; i < cells_marked.size(); i++)
//     {
//         assert(cells_marked[i - 1]->id < cells_marked[i]->id);
//     }

//     int counter = 0;
//     map<unsigned int, Cell *> split_cells;
//     for (int i = 0; i < cells_marked.size(); i++)
//     {
//         if (cells_marked[i]->isleaf)
//         {
//             map<unsigned int, Cell *> split_cells0 =
//             mesh.refinement(cells_marked[i]);
//             split_cells.insert(split_cells0.begin(), split_cells0.end());
//         }
//     }

//     expand_G(split_cells);

//     map<unsigned int, Cell *>::const_iterator map_it = split_cells.begin();
//     Cell *parent_cell;
//     for (int i = 0; i < split_cells.size(); i++)
//     {

//         int index = map_it->first;
//         parent_cell = map_it->second;
//         for (int j = 0; j < 8; j++)
//         {
//             assert(parent_cell->child_cells[j]->isleaf);
//             double rc, thetac, phic;
//             parent_cell->child_cells[j]->get_center(rc, thetac, phic);

//             double latc = 90.0 - thetac * 180.0 / GS::PI;
//             double lonc = phic * 180.0 / GS::PI - 180.0;
//             double depthc = (reference_surface - rc) / 1000.0;
//             array<double, 3> args = {latc, lonc, depthc};

//             parent_cell->child_cells[j]->set_parameter(interp_ML_m0.interp(args.begin()),
//             1);
//             parent_cell->child_cells[j]->set_parameter(interp_ML_crg.interp(args.begin()),
//             2);
//         }
//         map_it++;
//     }
//     // cout << G.rows() << ", " << G.cols() << endl;
//     mesh.get_model_parameter_from_mesh(m, 0);
//     mesh.get_model_parameter_from_mesh(m0, 1);
//     mesh.get_model_parameter_from_mesh(m0_s, 2);
//     mesh.get_model_parameter_from_mesh(m_min, 3);
//     mesh.get_model_parameter_from_mesh(m_max, 4);
//     mesh.rearrange_id();
//     this->init_matrices();

//     // cout << G.rows() << ", " << G.cols() << endl;
//     // this->compute_G();
// }

void AdaptiveInversion::refine_mesh(double a)
{
    if (a < 0 || abs(a * Nm) < 1)
        return;
    this->set_density_to_mesh();
    this->set_reference_model_to_mesh();
    this->set_min_max_to_mesh();

    VectorXd indicator(Nm);
    assert(Nm == G.cols());
    assert(Nm == m.rows());

    this->indicator_calculator(indicator);
    // cout<<indicator<<endl;

    VectorXd sorted_indicator;
    VectorXi sorted_index;
    sort_vec(indicator, sorted_indicator, sorted_index);

    vector<Cell *> cells_marked(0);
    cells_marked.clear();
    assert(Nm == mesh.leaf_cells.size());

    double tol = sorted_indicator(int(a * Nm));
    // cout<<tol<<endl;
    // cout<<a * Nm<<endl;
    for (int i = 0; i < Nm; i++)
    {
        if (indicator(i) > tol)
        {
            Cell *c = this->mesh.leaf_cells[i];
            double c_dx;
            double c_dy;
            double c_dz;
            c->get_size(c_dx, c_dy, c_dz);
            if (c_dx < min_dx || std::abs(c_dx - min_dx) < 1e-7 ||
                c_dy < min_dy || std::abs(c_dy - min_dy) < 1e-7 ||
                c_dz < min_dz || std::abs(c_dz - min_dz) < 1e-7)
            {
                // If the size of the cell is less than the limit of minimum
                // size, do nothing
            }
            else
            {
                cells_marked.push_back(c);
            }
        }
    }
    // cout<<cells_marked.size()<<endl;

    for (int i = 1; i < cells_marked.size(); i++)
    {
        assert(cells_marked[i - 1]->id < cells_marked[i]->id);
    }

    int counter = 0;
    map<unsigned int, Cell *> split_cells;
    for (int i = 0; i < cells_marked.size(); i++)
    {
        if (cells_marked[i]->isleaf)
        {
            map<unsigned int, Cell *> split_cells0 =
                mesh.refinement(cells_marked[i]);
            split_cells.insert(split_cells0.begin(), split_cells0.end());
        }
    }
    if (split_cells.size() > 0)
    {
        expand_G(split_cells);
    }

    mesh.get_model_parameter_from_mesh(m, 0);
    mesh.get_model_parameter_from_mesh(m0, 1);
    mesh.get_model_parameter_from_mesh(m0_s, 2);
    mesh.get_model_parameter_from_mesh(m_min, 3);
    mesh.get_model_parameter_from_mesh(m_max, 4);
    mesh.rearrange_id();
    this->init_matrices();
}
