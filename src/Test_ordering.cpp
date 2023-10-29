#include "GaussNewtonInversion.h"
#include "Mesh.h"

int main()
{
    double x_lim[2] = {200, 500};
    double y_lim[2] = {900, 1300};
    double z_lim[2] = {0, 400};

    Mesh mesh;
    mesh.generate_regular_mesh(x_lim, 3, y_lim, 4, z_lim, 4);
    mesh.show_ordering();

    mesh.set_block_parameter(1, 1, 1, 2, 1, 1, 100);
    VectorXd rho;
    mesh.get_model_parameter_from_mesh(rho, 0);
    mesh.out_model_vtk("test_ordering_model.vtk");

    // test ordering after refinement
    int nlat = 3;
    int nlon = 3;
    double observation_r = 6428137;
    Observation ob;
    int nx = 3;
    int ny = 3;
    ob.generate(nx, x_lim[0], x_lim[1], ny, y_lim[0], y_lim[1], -0.1);

    cout << mesh.leaf_cells[0]->get_ordering_forward() << endl;

    GaussNewtonInversion inv(mesh, ob, Compute_g_z);
    Mesh &invmesh0 = inv.get_mesh();
    cout << invmesh0.leaf_cells[0]->get_ordering_forward() << endl;
    inv.compute_G();
    inv.set_m(rho);
    inv.set_density_to_mesh();
    inv.refine_mesh(0.05);
    inv.result2vtk("test_ordering_model_refined");
    Mesh &invmesh = inv.get_mesh();
    invmesh.show_ordering();

    return 0;
}
