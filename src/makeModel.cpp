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

namespace po = boost::program_options;

// A helper function to simplify the main part.
template<class T>
ostream& operator<<(ostream& os, const vector<T>& v)
{
    copy(v.begin(), v.end(), ostream_iterator<T>(os, " "));
    return os;
}

int main(int ac, char *av[])
{
    try
    {
        double x_lim[2] = {0, 0};
        double y_lim[2] = {0, 0};
        double z_lim[2] = {0, 0};
        int nx,ny,nz;

        int opt;
        int portnum;
        vector<string> include_path;
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message") 
            ("region_x,x",po::value<string>(), "specify the region of interet in the x direction. For example, 0/6000/300 specify a range from 0 to 6000 with a spacing of 300.")
            ("region_y,y",po::value<string>(), "specify the region of interet in the y direction")
            ("region_z,z",po::value<string>(), "specify the region of interet in the z direction")
            ("nx",po::value<int>(&nx), "number of cell in the x direction.")
            ("ny",po::value<int>(), "number of cell in the y direction.")            
            ("nz",po::value<int>(), "number of cell in the z direction.")            
            ("anomaly_file,a",po::value<string>(),"specify anomalies from command line")
            ("Anomaly_file,A",po::value<string>(),"specify the file describing anomalous bodies")
            ("output,o",po::value<string>(), "specify the file to save the model")
            // ("calculate,c",po::value<string>(),"calculate gravity")
            // ("observation,o",po::value<string>(),"specify observation points")
            ("optimization", po::value<int>(&opt)->default_value(10),
                  "optimization level")
            ("verbose,v", po::value<int>()->implicit_value(1),
                  "enable verbosity (optionally specify level)")
            ("listen,l", po::value<int>(&portnum)->implicit_value(1001)
                  ->default_value(0,"no"),
                  "listen on a port.")
            ("include-path,I", po::value< vector<string> >(&include_path),
                  "include path")
            // ("input-file", po::value< vector<string> >(), "input file")
        ;

        po::positional_options_description p;
        // p.add("input-file", -1);

        po::variables_map vm;
        po::store(po::command_line_parser(ac, av).
                  options(desc).positional(p).run(), vm);
        po::notify(vm);

        if (vm.count("help")) {
            cout << "Usage: makeModel [options]\n";
            cout << desc;
            return 0;
        }
        if (vm.count("anomaly_file")&&vm.count("Anomaly_file")){
            cout << "Offending options: --anomaly_file[-a] and --Anomaly_file[-A] can not be used together"<<endl;
            return 1;
        }
        if (vm.count("anomaly_file")==0 && vm.count("Anomaly_file")==0){
            cout << "You need to specify at least one anomaly either from command line using option --anomaly_file[-a] or from a file using --Anomaly_file[-A]"<<endl;
            return 1;
        }
        else if(vm.count("anomaly_file")){

        }

        if (vm.count("include-path"))
        {
            cout << "Include paths are: "
                 << vm["include-path"].as< vector<string> >() << "\n";
            // include_path=vm["include-path"].as<vector<string>>();
            cout <<include_path.size()<<"\n";
            cout <<include_path<<"\n";
        }

        // if (vm.count("input-file"))
        // {
        //     cout << "Input files are: "
        //          << vm["input-file"].as< vector<string> >() << "\n";
        // }

        if (vm.count("verbose")) {
            cout << "Verbosity enabled.  Level is " << vm["verbose"].as<int>()
                 << "\n";
        }

        cout << "Optimization level is " << opt << "\n";

        cout << "Listen port is " << portnum << "\n";
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