#include <chrono>
#include <fstream>
#include <random>
// #include<function>

// #include "GaussNewtonInversion.h"
// #include "timer.h"

#include <boost/program_options.hpp>
#include <iostream>
#include <iterator>
using namespace std;

#include "Fwd.h"

namespace po = boost::program_options;

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

std::vector<std::string> split_to_string(const std::string &s, char delim)
{
    std::vector<std::string> result;
    std::stringstream ss(s);
    std::string item;

    while (getline(ss, item, delim))
    {
        result.push_back(item);
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
        int opt;
        int portnum;
        po::options_description desc("Allowed options");
        desc.add_options()("help,h", "produce help message")("model,m", po::value<string>(), "Specify the model file")("observation,p", po::value<string>(), "specify a file containing locations of observation points. Each line of this file should be three floating-point numbers corresponding to x y z coordinates of an observation point.")("component,c", po::value<string>()->default_value("V/gz/gx/gy/Txx/Txy/Txz/Tyy/Tyz/Tzz"), "specify the components of gravitational potential/vector/gradient tensor to be components. The possible values are: \n\tV: gravitational potential\n\tgx/gy/gz: x/y/z component of gravity filed\n\tTxx/Txy/Txz/Tyy/Tyz/Tzz.\nFor example, the option '-c gx/gz/Txx/Tyz/Tzz' makes the program to calculate gx, gz, Txx, Tyz and Tzz.")("output,o", po::value<string>(), "output file");

        po::positional_options_description p;
        // p.add("input-file", -1);

        po::variables_map vm;
        po::store(po::command_line_parser(ac, av).options(desc).positional(p).run(), vm);
        po::notify(vm);

        if (vm.count("help"))
        {
            cout << "Usage: options_description [options]\n";
            cout << desc;
            return 0;
        }
        vector<RectPrism> rects;
        vector<double> density;
        if (vm.count("model"))
        {
            string model_file = vm["model"].as<string>();
            ifstream input_stream(model_file.c_str());
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
                    double x0, y0, z0, x1, y1, z1, xc, yc, zc, rho;
                    std::istringstream iss(line);
                    iss >> x0 >> x1 >> y0 >> y1 >> z0 >> z1 >> xc >> yc >> zc >> rho;

                    rects.push_back(RectPrism(x0, y0, z0, x1, y1, z1));
                    density.push_back(rho);
                }
            }
        }
        else
        {
        }

        Observation ob;
        if (vm.count("observation"))
        {
            string observation_file = vm["observation"].as<string>();
            ifstream input_stream(observation_file.c_str());
            assert(input_stream.good());
            string line;

            int n_obs = 0;
            while (std::getline(input_stream, line))
            {
                line_process(line);
                if (line.empty())
                {
                    continue;
                }
                else
                {
                    n_obs++;
                    double x0, y0, z0;
                    std::istringstream iss(line);
                    iss >> x0 >> y0 >> z0;
                    ob.add_point(x0, y0, z0);
                }
            }
        }

        unsigned long long field_flag_ULL;

        if (vm.count("component"))
        {
            string field_str = vm["component"].as<string>();
            vector<string> component_strings = split_to_string(field_str, '/');
            cout << "The following field components will be computed: " << endl;
            for (int i = 0; i < component_strings.size(); i++)
            {
                cout << component_strings[i] << " ";
            }
            cout << endl;
            for (int i = 0; i < component_strings.size(); i++)
            {
                if (component_strings[i] == "V")
                {
                    field_flag_ULL = field_flag_ULL | Compute_V;
                }
                else if (component_strings[i] == "gz")
                {
                    field_flag_ULL = field_flag_ULL | Compute_g_z;
                }
                else if (component_strings[i] == "gx")
                {
                    field_flag_ULL = field_flag_ULL | Compute_g_x;
                }
                else if (component_strings[i] == "gy")
                {
                    field_flag_ULL = field_flag_ULL | Compute_g_y;
                }
                else if (component_strings[i] == "Txx" || component_strings[i] == "gxx")
                {
                    field_flag_ULL = field_flag_ULL | Compute_T_xx;
                }
                else if (component_strings[i] == "Txy" || component_strings[i] == "Tyx" || component_strings[i] == "gxy" || component_strings[i] == "gyx")
                {
                    field_flag_ULL = field_flag_ULL | Compute_T_xy;
                }
                else if (component_strings[i] == "Txz" || component_strings[i] == "Tzx" || component_strings[i] == "gxz" || component_strings[i] == "gzx")
                {
                    field_flag_ULL = field_flag_ULL | Compute_T_zx;
                }
                else if (component_strings[i] == "Tyy" || component_strings[i] == "gyy")
                {
                    field_flag_ULL = field_flag_ULL | Compute_T_yy;
                }
                else if (component_strings[i] == "Tyz" || component_strings[i] == "Tzy" || component_strings[i] == "gyz" || component_strings[i] == "gzy")
                {
                    field_flag_ULL = field_flag_ULL | Compute_T_zy;
                }
                else if (component_strings[i] == "Tzz" || component_strings[i] == "gzz")
                {
                    field_flag_ULL = field_flag_ULL | Compute_T_zz;
                }
                else
                {
                    cout << "Invalid component name" << endl;
                    return 1;
                }
            }
        }
        bitset<10> field_flag(field_flag_ULL);

        vector<unsigned int> field_index;
        for (unsigned int i = 0; i < 10; i++)
        {
            if (field_flag[i])
            {
                field_index.push_back(i);
            }
        }
        vector<vector<double>> results;
        results.resize(ob.get_n_obs());
        unsigned int n_components = field_index.size();
        for (unsigned int i = 0; i < ob.get_n_obs(); i++)
        {

            results[i].resize(n_components);
            for (unsigned int k = 0; k < n_components; k++)
            {
                results[i][k] == 0.0;
            }
        }

#pragma omp parallel for
        for (unsigned int i = 0; i < ob.get_n_obs(); i++)
        {
            for (unsigned int j = 0; j < rects.size(); j++)
            {
                GravFormula gra;
                vector<double> field;
                gra.field_caused_by_single_prism(ob(i), rects[j], density[j], field,
                                                 field_flag);
                for (unsigned int k = 0; k < n_components; k++)
                {
                    results[i][k] += field[field_index[k]];
                }
            }
        }

        string outfile = vm["output"].as<string>();
        ofstream out_s(outfile);
        vector<string> field_names_in_header = {"V (m2/s2)", "g_z (mGal)", "g_x (mGal)", "g_y (mGal)", "T_zz (E)",
                                                "T_xz (E)", "T_yz (E)", "T_xx (E)", "T_xy (E)", "T_yy (E)"};
        out_s << setw(23) << left << "# x" << setw(23) << left << "y" << setw(23) << left << "z";
        for (int i = 0; i < n_components; i++)
        {
            out_s << setw(23) << left << field_names_in_header[field_index[i]];
        }
        out_s << endl;

        for (int i = 0; i < ob.get_n_obs(); i++)
        {
            const Point &p = ob(i);
            out_s << scientific;
            out_s << setw(23) << setprecision(15) << left << p.x() << setw(23) << left
                  << p.y() << setw(23) << setprecision(15) << left << p.z();
            out_s << scientific;
            for (unsigned int k = 0; k < n_components; k++)
            {
                out_s << setw(23) << setprecision(15) << left << results[i][k];
            }
            out_s << endl;
        }
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