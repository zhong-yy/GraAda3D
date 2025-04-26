#include <chrono>
#include <fstream>
#include <random>
// #include<function>

// #include "GaussNewtonInversion.h"
// #include "timer.h"

#include <boost/program_options.hpp>
#include <iostream>
#include <iterator>
#include "Mesh.h"
using namespace std;

namespace po = boost::program_options;

// A helper function to simplify the main part.
template <class T>
ostream &operator<<(ostream &os, const vector<T> &v)
{
    copy(v.begin(), v.end(), ostream_iterator<T>(os, " "));
    return os;
}

std::vector<double> split_to_double(const std::string &s, char delim)
{
    std::vector<double> result;
    std::stringstream ss(s);
    std::string item;

    while (getline(ss, item, delim))
    {
        result.push_back(std::stod(item));
    }

    return result;
}

void line_process(std::string &line, const std::string comment_str = "#")
{
    for (char &c : line) // C++11以上版本的语法
    {
        // 制表符 tab，逗号，分号都当作有效的分隔符，统一转成空格
        // 为了避免错误，回车符和换行符也转为空格（否则无法处理空行）
        if (c == '\t' || c == ',' || c == ';' || c == '\r' || c == '\n')
            c = ' ';
    }

    // 查找注释符所在位置，如果不存在，则得到string::npos
    int n_comment_start = line.find_first_of(comment_str);
    if (n_comment_start != std::string::npos) // 这一句必须的
        line.erase(n_comment_start);          // 删除注释

    line.erase(0, line.find_first_not_of(" ")); // 删除行首空格
    line.erase(line.find_last_not_of(" ") + 1); // 删除行末空格
    // 调用的string& erase (size_t pos = 0, size_t len = npos);
    //  len为默认参数,
    //  size_t可以当作无符号整数，npos是string内置的静态常量，为size_t的最大值

    /****************************************************************
     *处理完毕。如果这一行只有空格，制表符 tab，注释，那么处理后line为空；
     *如果行首有多个空格(或者空格和tab交错)，行尾为注释，如
     *“   a b c#坐标”
     *那么处理后字符串line的行首多个空格(和tab)和行尾注释被删掉，得到
     *“a b c”
     ****************************************************************/
}

int main(int ac, char *av[])
{
    try
    {
        int nx, ny, nz;

        int opt;
        string region;
        po::options_description desc("Allowed options");
        desc.add_options()("help,h", "show help message")("region,R", po::value<string>(&region), "specify the region of interet. For example, -R 0/4000/0/2000/0/1000 defines a region with x ranging from 0 to 4000, y from 0 to 2000 and z from 0 to 1000.")
            // ("region_x,x",po::value<string>(), "specify the region of interet in the x direction. For example, 0/6000 specify a range from 0 to 6000.")
            // ("region_y,y",po::value<string>(), "specify the region of interet in the y direction")
            // ("region_z,z",po::value<string>(), "specify the region of interet in the z direction")
            ("nx", po::value<int>(&nx), "number of cells in the x direction.")("ny", po::value<int>(&ny), "number of cells in the y direction.")("nz", po::value<int>(&nz), "number of cells in the z direction.")("anomalies,a", po::value<vector<string>>(), "specify anomalies from command line. For example, -a 200/500/600/900/550/850/-300 means adding an anomalous body, with x dimension [200, 500], y dimension [600, 900] and z dimension [550, 850]. The physical property is -300.")("anomaly_file,A", po::value<string>(), "specify the file describing anomalous bodies")("output,o", po::value<string>(), "specify the file to save the model")
            // ("calculate,c",po::value<string>(),"calculate gravity")
            // ("observation,o",po::value<string>(),"specify observation points")
            // ("input-file", po::value< vector<string> >(), "input file")
            ;

        po::positional_options_description p;

        po::variables_map vm;
        po::store(po::command_line_parser(ac, av).options(desc).positional(p).run(), vm);
        po::notify(vm);

        if (vm.count("help"))
        {
            cout << "Usage: makeModel [options]\n";
            cout << desc;
            return 0;
        }

        double x_lim[2] = {0, 0};
        double y_lim[2] = {0, 0};
        double z_lim[2] = {0, 0};
        Mesh mesh;

        if (vm.count("region") && vm.count("nx") && vm.count("ny") && vm.count("nz"))
        {
            vector<double> region_params = split_to_double(region, '/');
            x_lim[0] = region_params[0];
            x_lim[1] = region_params[1];
            y_lim[0] = region_params[2];
            y_lim[1] = region_params[3];
            z_lim[0] = region_params[4];
            z_lim[1] = region_params[5];
            cout << "X range: [" << x_lim[0] << ", " << x_lim[1] << "]" << ", " << nx << " cells" << endl;
            cout << "Y range: [" << y_lim[0] << ", " << y_lim[1] << "]" << ", " << ny << " cells" << endl;
            cout << "Z range: [" << z_lim[0] << ", " << z_lim[1] << "]" << ", " << nz << " cells" << endl;

            mesh.generate_regular_mesh(x_lim, nx, y_lim, ny, z_lim, nz);
        }
        if (vm.count("anomalies") && vm.count("anomaly_file"))
        {
            cout << "Offending options: --anomaly_file[-a] and --anomaly_file[-A] can not be used together" << endl;
            return 1;
        }
        if (vm.count("anomalies") == 0 && vm.count("anomaly_file") == 0)
        {
            cout << "You need to specify at least one anomaly either from command line using option --anomaly_file[-a] or from a file using --anomaly_file[-A]" << endl;
            return 1;
        }
        else if (vm.count("anomaly_file"))
        {
            string anomaly_file = vm["anomaly_file"].as<string>();
            ifstream input_stream(anomaly_file.c_str());
            if (!input_stream.good())
            {
                cout << "Please check whether file " << anomaly_file << " exists." << endl;
                return 1;
            }
            assert(input_stream.good());
            string line;
            while (std::getline(input_stream, line))
            {
                line_process(line);
                if (line.empty())
                {
                    continue;
                }
                else
                {
                    double x0, y0, z0, x1, y1, z1, rho;
                    std::istringstream iss(line);
                    iss >> x0 >> x1 >> y0 >> y1 >> z0 >> z1 >> rho;
                    double block_z[2] = {0, 0};
                    double block_x[2] = {0, 0};
                    double block_y[2] = {0, 0};
                    block_x[0] = x0;
                    block_x[1] = x1;
                    block_y[0] = y0;
                    block_y[1] = y1;
                    block_z[0] = z0;
                    block_z[1] = z1;
                    double value = rho;
                    mesh.set_parameter_in_a_region(block_x, block_y, block_z, value);
                }
            }
        }
        else if (vm.count("anomalies"))
        {
            vector<string> anomalies = vm["anomalies"].as<vector<string>>();
            for (auto ano : anomalies)
            {
                double block_z[2] = {0, 0};
                double block_x[2] = {0, 0};
                double block_y[2] = {0, 0};

                vector<double> ano_params = split_to_double(ano, '/');
                block_x[0] = ano_params[0];
                block_x[1] = ano_params[1];
                block_y[0] = ano_params[2];
                block_y[1] = ano_params[3];
                block_z[0] = ano_params[4];
                block_z[1] = ano_params[5];
                double value = ano_params[6];
                mesh.set_parameter_in_a_region(block_x, block_y, block_z, value);
            }
        }
        string output_filename = vm["output"].as<string>();
        mesh.out_model_vtk(output_filename + ".vtk");
        mesh.out_model_txt(output_filename + ".xyz");
#ifdef USE_NETCDF
        mesh.out_model_netcdf(output_filename+".nc");
#endif
    }
    catch (exception &e)
    {
        cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch (...)
    {
        cerr << "Exception of unknown type!\n";
    }
    return 0;
}