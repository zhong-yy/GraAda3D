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
        result.push_back(std::stod(item));
    }

    return result;
}

int main(int ac, char *av[])
{
    try
    {
        int opt;
        int portnum;
        po::options_description desc("Allowed options");
        desc.add_options()
        ("help,h", "produce help message")
        ("model,m", po::value<string>(), "Specify the model file")
        ("observation,p", po::value<string>(), "specify a file containing locations of observation points. Each line of this file should be three floating-point numbers corresponding to x y z coordinates of an observation point.")
        ("component,c", po::value<string>(&field_str)->default_value("V/gz/gx/gy/Txx/Txy/Txz/Tyy/Tyz/Tzz"), "specify the components of gravitational potential/vector/gradient tensor to be components. The possible values are: \n\tV: gravitational potential\n\tgx/gy/gz: x/y/z component of gravity filed\n\tTxx/Txy/Txz/Tyy/Tyz/Tzz.\nFor example, the option '-c gx/gz/Txx/Tyz/Tzz' makes the program to calculate gx, gz, Txx, Tyz and Tzz.")("output,o", po::value<string>(), "output file");

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

                }
            }
        }
        else
        {
            
        }

        Observation ob;
        if (vm.count("observation")){
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
                    iss >> x0 >> y0>>z0;
                    ob.add_point(x0, y0, z0);
                }
            }            
        }

        unsigned long field_flag;
        if (vm.count("component"))
        {
            // string field_str = vm["component"].as<string>();
            vector<string> component_strings = split_to_string(field_str, '/');
            cout << "The following field components will be computed: " << endl;
            for (int i = 0; i < component_strings.size(); i++)
            {
                cout << component_strings[i] << " ";
            }
            cout << endl;
            for (int i = 0; i < component_strings.size(); i++)
            {
                switch (component_strings[i])
                {
                case "V":
                    field_flag = field_flag | Compute_V;
                    break;
                case "gz":
                    field_flag = field_flag | Compute_g_z;
                    break;
                case "gx":
                    field_flag = field_flag | Compute_g_x;
                    break;
                case "gy":
                    field_flag = field_flag | Compute_g_y;
                    break;                    
                case "Txx":
                    field_flag = field_flag | Compute_T_xx;
                    break;
                case "Txy":
                    field_flag = field_flag | Compute_T_xy;
                    break;
                case "Tyx":
                    field_flag = field_flag | Compute_T_xy;
                    break;
                case "Txz":
                    field_flag = field_flag | Compute_T_zx;
                    break;
                case "Tzx":
                    field_flag = field_flag | Compute_T_zx;
                    break;
                case "Tyy":
                    field_flag = field_flag | Compute_T_yy;
                    break;
                case "Tyz":
                    field_flag = field_flag | Compute_T_zy;
                    break;
                case "Tzy":
                    field_flag = field_flag | Compute_T_zy;
                    break;                    
                case "Tzz":
                    field_flag = field_flag | Compute_T_zz;
                    break;
                default:
                    cout<<"Invalid component name"<<endl;
                    return 1;
                    break;
                }
            }
        }

        GravFormula gra;
        vector<double> field;
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