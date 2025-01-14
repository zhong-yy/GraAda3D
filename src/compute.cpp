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

int main(int ac, char* av[]){
   try
    {
        int opt;
        int portnum;
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message")
            ("model,r",po::value<string>(), "Specify the model file")
            ("observation,p",po::value<string>(),"specify observation points")
            ("fiedl,f",po::value<string>(),"specify observation points")
            ("outpu,f",po::value<string>(),"specify observation points")
            ("optimization", po::value<int>(&opt)->default_value(10),
                  "optimization level")
            ("verbose,v", po::value<int>()->implicit_value(1),
                  "enable verbosity (optionally specify level)")
            ("listen,l", po::value<int>(&portnum)->implicit_value(1001)
                  ->default_value(0,"no"),
                  "listen on a port.")
            ("include-path,I", po::value< vector<string> >(),
                  "include path")
            ("input-file", po::value< vector<string> >(), "input file")
        ;

        po::positional_options_description p;
        p.add("input-file", -1);

        po::variables_map vm;
        po::store(po::command_line_parser(ac, av).
                  options(desc).positional(p).run(), vm);
        po::notify(vm);

        if (vm.count("help")) {
            cout << "Usage: options_description [options]\n";
            cout << desc;
            return 0;
        }

        if (vm.count("include-path"))
        {
            cout << "Include paths are: "
                 << vm["include-path"].as< vector<string> >() << "\n";
        }

        if (vm.count("input-file"))
        {
            cout << "Input files are: "
                 << vm["input-file"].as< vector<string> >() << "\n";
        }

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