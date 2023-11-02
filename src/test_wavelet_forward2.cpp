#include <chrono>
#include <fstream>
#include <random>
// #include<function>

#include "GaussNewtonInversion.h"
#include "timer.h"
void func1();
void func2(double);
int main(int argc, char **argv)
{
    if (argc != 2)
    {
        cout << "A thresholding parameter should follow the program name. "
                "Usage:\n";
        cout << "program_name [relative threshold]" << endl;
        return 1;
    }
    double relative_threshold = atof(argv[1]);
    func1();
    func2(relative_threshold);
    return 0;
}
void func1()
{
    double x_lim[2] = {0, 2000};
    double y_lim[2] = {0, 2000};
    double z_lim[2] = {0, 1000};

    Mesh mesh;
    mesh.generate_regular_mesh(x_lim, 40, y_lim, 40, z_lim, 20);

    double block_z[2] = {200, 500};
    double block_x[2] = {600, 900};
    double block_y[2] = {550, 850};
    mesh.set_parameter_in_a_region(block_x, block_y, block_z, 500);

    double block_z2[2] = {200, 500};
    double block_x2[2] = {1200, 1500};
    double block_y2[2] = {550, 850};
    mesh.set_parameter_in_a_region(block_x2, block_y2, block_z2, -500);

    double block_z3[2] = {200, 500};
    double block_x3[2] = {900, 1500};
    double block_y3[2] = {1300, 1500};
    mesh.set_parameter_in_a_region(block_x3, block_y3, block_z3, 500);

    //     mesh.out_model_vtk("test_model.vtk");
    // #ifdef USE_NETCDF
    //     mesh.out_model_netcdf("test_model.nc");
    // #endif
    //     mesh.out_model_txt("test_model.txt");

    /******************Foward modelling*****************/
    VectorXd rho;
    mesh.get_model_parameter_from_mesh(rho, 0);
    Observation ob;

    int nx = 41;
    int ny = 41;
    ob.generate(nx, x_lim[0], x_lim[1], ny, y_lim[0], y_lim[1], -0.1);
    ofstream os("sites");
    os << ob;

    GaussNewtonInversion temp_var(mesh, ob, Compute_g_z);
    Mesh &temp_mesh0 = temp_var.get_mesh();
    temp_var.compute_G();
    temp_var.set_m(rho);
    temp_var.set_density_to_mesh();
    temp_var.refine_mesh(0.05);
    // temp_var.out_model_vtk("test_model.vtk");
    Mesh &mesh1 = temp_var.get_mesh();
    mesh1.out_model_vtk("test_model.vtk");

    mesh1.get_model_parameter_from_mesh(rho, 0);
    // temp_mesh.show_ordering();

    Fwd forward(mesh1, ob, Compute_g_z);
    // Fwd forward(mesh, ob,
    // Compute_T_rr|Compute_T_rtheta|Compute_T_rphi|Compute_T_thetatheta|Compute_T_thetaphi|Compute_T_phiphi,8);
    // Fwd forward(&mesh, &ob, Compute_T_rr|Compute_T_rtheta, 8);

    Timer timer;
    timer.start();
    cout << "Generating synthetic data..." << endl;
    forward.compute_G();
    const Eigen::MatrixXd &G = forward.get_G();
    timer.stop();
    cout << "Time for calculating sensitivity: " << timer.getElapsedTimeInSec()
         << " s" << endl;
    timer.start();
    Eigen::VectorXd d_obs = G * rho;
    cout << "Calculation of synthetic data completed" << endl;
    timer.stop();
    cout << "Multiplication time: " << timer.getElapsedTimeInSec() << " s"
         << endl;

    VectorXd test_vec(d_obs.size());
    for (int i = 0; i < test_vec.size(); ++i)
    {
        test_vec(i) = 1;
    }
    timer.start();
    VectorXd GTv1 = G.transpose() * test_vec;
    timer.stop();
    cout << "Multiplication time of GT*v: " << timer.getElapsedTimeInSec() << "s"
         << endl;

    ofstream ofs("GTv1.txt");
    Mesh m = forward.get_mesh();
    for (int i = 0; i < GTv1.size(); i++)
    {
        int id = m.get_reordered_id(i);
        ofs << setw(23) << scientific << GTv1(id);
        ofs << endl;
    }

    // cout << "||d_obs||=" << d_obs.norm() << endl;
    // cout << "0.001*d_obs.norm()=" << 0.001 * d_obs.norm() << endl;

    // GaussNewtonInversion inv(mesh, ob,
    // Compute_T_rr|Compute_T_rtheta|Compute_T_rphi|Compute_T_thetatheta|Compute_T_thetaphi|Compute_T_phiphi);
    GaussNewtonInversion inv(mesh1, ob, Compute_g_z);

    inv.set_dobs(d_obs);
    inv.output_obs_data("dobs");
}

void func2(double relative_threshold)
{
    double x_lim[2] = {0, 2000};
    double y_lim[2] = {0, 2000};
    double z_lim[2] = {0, 1000};
    double reference_surface = 6378137;

    Mesh mesh;
    mesh.generate_regular_mesh(x_lim, 40, y_lim, 40, z_lim, 20);

    double block_z[2] = {200, 500};
    double block_x[2] = {600, 900};
    double block_y[2] = {550, 850};
    mesh.set_parameter_in_a_region(block_x, block_y, block_z, 500);

    double block_z2[2] = {200, 500};
    double block_x2[2] = {1200, 1500};
    double block_y2[2] = {550, 850};
    mesh.set_parameter_in_a_region(block_x2, block_y2, block_z2, -500);

    double block_z3[2] = {200, 500};
    double block_x3[2] = {900, 1500};
    double block_y3[2] = {1300, 1500};
    mesh.set_parameter_in_a_region(block_x3, block_y3, block_z3, 500);

    //     mesh.out_model_vtk("test_model.vtk");
    // #ifdef USE_NETCDF
    //     mesh.out_model_netcdf("test_model.nc");
    // #endif
    //     mesh.out_model_txt("test_model.txt");

    VectorXd rho;
    mesh.get_model_parameter_from_mesh(rho, 0);
    Observation ob;

    int nx = 41;
    int ny = 41;
    ob.generate(nx, x_lim[0], x_lim[1], ny, y_lim[0], y_lim[1], -0.1);
    ofstream os("sites");
    os << ob;

    //***********************************WAVELET
    GaussNewtonInversion temp_var(mesh, ob, Compute_g_z);
    Mesh &temp_mesh0 = temp_var.get_mesh();
    temp_var.compute_G();
    temp_var.set_m(rho);
    temp_var.set_density_to_mesh();
    temp_var.refine_mesh(0.05);
    // temp_var.out_model_vtk("test_model.vtk");
    Mesh &mesh1 = temp_var.get_mesh();
    mesh1.out_model_vtk("test_model.vtk");

    mesh1.get_model_parameter_from_mesh(rho, 0);

    Fwd forward_wl(mesh1, ob, Compute_g_z);
    forward_wl.set_use_wavelet(true);
    Timer timer2;
    timer2.start();
    forward_wl.set_compression_threshold(relative_threshold);
    cout << endl;
    cout
        << "=====================Using wavelet transfrom======================="
        << endl;
    cout << "Using wavelet transform\nGenerating synthetic data..." << endl;
    forward_wl.compute_G_wavelet();
    timer2.stop();
    cout << "Time for calculating sensitivity: " << timer2.getElapsedTimeInSec()
         << " s" << endl;

    VectorXd d_obs_wl;
    timer2.start();
    // cout << rho << endl;
    forward_wl.G_vec_mul(rho, d_obs_wl);
    timer2.stop();
    cout << "Multiplication time: " << timer2.getElapsedTimeInSec() << " s"
         << endl;

    VectorXd test_vec(ob.get_n_obs());
    for (int i = 0; i < test_vec.size(); ++i)
    {
        test_vec(i) = 1;
    }
    timer2.start();
    VectorXd GTv2;
    forward_wl.GT_vec_mul(test_vec, GTv2);
    timer2.stop();
    cout << "Multiplication time of GT*v using wavelet compression: "
         << timer2.getElapsedTimeInSec() << endl;

    ofstream ofs("GTv2.txt");
    Mesh m = forward_wl.get_mesh();
    for (int i = 0; i < GTv2.size(); i++)
    {
        int id = m.get_reordered_id(i);
        ofs << setw(23) << scientific << GTv2(id);
        ofs << endl;
    }

    GaussNewtonInversion inv2(mesh, ob, Compute_g_z);
    inv2.set_dobs(d_obs_wl);
    inv2.output_obs_data("dobs_wavelet");
}
