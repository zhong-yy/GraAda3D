#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <cstdlib>
using namespace std;
#define EPS 1e-14

#include <boost/program_options.hpp>
namespace po = boost::program_options;

int main(int ac, char *av[])
{
	// if(argc <5) {
	// 	printf("Usage: %s input_paramters_filename\n", argv[0]);
	// 	return 1;
	// }
	// first string: x
	// second string: y
	// third string: z
	try
	{
		po::options_description desc("Allowed options");
		desc.add_options()("help,h", "produce help message")
		("xcoordinates,x", po::value<string>(), "x coordinates. For example, '-x 0:2.5:10' will generate a list of x coordinates starting at 0 and ending at 10, with an interval of 2.5. If only one number is given, e.g '-x 5', there will be only one value for x. ")
		("ycoordinates,y", po::value<string>(), "y coordinates. The usage is similar to x.")
		("zcoordinates,z", po::value<string>(), "z coordinates. The usage is similar to z.")
		("output,o", po::value<string>(), "output file");

		po::positional_options_description p;

		po::variables_map vm;
		po::store(po::command_line_parser(ac, av).options(desc).positional(p).run(), vm);
		po::notify(vm);
        if (vm.count("help")) {
            cout << "Usage: makeModel [options]\n";
            cout << desc;
            return 0;
        }
		
		typedef vector<double>::iterator iter;
		int num = 0;
		string xstring = vm["xcoordinates"].as<string>();
		string ystring = vm["ycoordinates"].as<string>();
		string zstring = vm["zcoordinates"].as<string>();

		// parse x y z
		double start_x, xspace, end_x;
		double start_y, yspace, end_y;
		double start_z, zspace, end_z;
		cout<<"x"<<endl;
		if (vm.count("xcoordinates")==0){
			cout<<"You need to specify the x coordinates;"<<endl;
			return 1;
		}

		if (vm.count("ycoordinates")==0){
			cout<<"You need to specify the y coordinates;"<<endl;
			return 1;
		}

		if (vm.count("zcoordinates")==0){
			cout<<"You need to specify the z coordinates;"<<endl;
			return 1;
		}
		cout<<"Generating spatial points ..."<<endl;

		if (xstring.find(':') == string::npos)
		{
			start_x = std::stod(xstring);
			end_x = std::stod(xstring);
			xspace = 1;
			cout << setw(30) << left << "x coordinate:" << start_x << endl;
		}
		else
		{
			if (xstring.find(':') == xstring.rfind(':'))
			{
				start_x = std::stod(xstring.substr(0, xstring.find(':')));
				xspace = 1;
				end_x = std::stod(xstring.substr(xstring.rfind(':') + 1));
			}
			else
			{
				start_x = std::stod(xstring.substr(0, xstring.find(':')));
				xspace = std::stod(xstring.substr(xstring.find(':') + 1, xstring.rfind(':') - xstring.find(':')));
				end_x = std::stod(xstring.substr(xstring.rfind(':') + 1));
			}
			cout << setw(30) << left << "Starting x coordinate:" << start_x << endl;
			cout << setw(30) << left << "x interval:" << xspace << endl;
			cout << setw(30) << left << "Ending x coordinate:" << end_x << endl;
		}

		cout << endl;

		if (ystring.find(':') == string::npos)
		{
			start_y = stod(ystring);
			end_y = stod(ystring);
			yspace = 1;
			cout << setw(30) << left << "y coordinate:" << start_y << endl;
		}
		else
		{
			if (ystring.find(':') == ystring.rfind(':'))
			{
				start_y = std::stod(ystring.substr(0, ystring.find(':')));
				yspace = 1;
				end_y = std::stod(ystring.substr(ystring.rfind(':') + 1));
			}
			else
			{
				start_y = std::stod(ystring.substr(0, ystring.find(':')));
				yspace = std::stod(ystring.substr(ystring.find(':') + 1, ystring.rfind(':') - ystring.find(':')));
				end_y = std::stod(ystring.substr(ystring.rfind(':') + 1));
			}
			cout << setw(30) << left << "Starting y coordinate:" << start_y << endl;
			cout << setw(30) << left << "y interval:" << yspace << endl;
			cout << setw(30) << left << "Ending y coordinate:" << end_y << endl;
		}

		cout << endl;

		if (zstring.find(':') == string::npos)
		{
			start_z = std::stod(zstring);
			end_z = std::stod(zstring);
			zspace = 1;
			cout << setw(30) << left << "z coordinate:" << start_y << endl;
		}
		else
		{
			if (zstring.find(':') == zstring.rfind(':'))
			{
				start_z = std::stod(zstring.substr(0, zstring.find(':')));
				zspace = 1;
				end_z = std::stod(zstring.substr(zstring.rfind(':') + 1));
			}
			else
			{
				start_z = std::stod(zstring.substr(0, zstring.find(':')));
				zspace = std::stod(zstring.substr(zstring.find(':') + 1, zstring.rfind(':') - zstring.find(':')));
				end_z = std::stod(zstring.substr(zstring.rfind(':') + 1));
			}
			cout << setw(30) << left << "Starting z coordinate:" << start_z << endl;
			cout << setw(30) << left << "z interval:" << zspace << endl;
			cout << setw(30) << left << "Ending z coordinate:" << end_z << endl;
		}

		// double xspace=0.05,yspace=0.05,start_x=-3,start_y=-3.5,end_x=6,end_y=3.5;
		//	double z=atof(argv[3]);
		//	cout <<'\n'<< setw(30) << left << "z coordinate:" << z << endl;
		vector<double> x;
		vector<double> y;
		vector<double> z;

		ofstream os;
		string outfile=vm["output"].as<string>();
		os.open(outfile.c_str());

		int cou = 0;
		double temp = start_x + cou * xspace;
		do
		{
			x.push_back(temp);
			cou++;
			temp = start_x + cou * xspace;
		} while ((temp < end_x) || (abs(temp - end_x) < EPS));

		cou = 0;
		temp = start_y + cou * yspace;
		do
		{
			y.push_back(temp);
			cou++;
			temp = start_y + cou * yspace;
		} while ((temp < end_y) || abs(temp - end_y) < EPS);

		cou = 0;
		temp = start_z + cou * zspace;
		do
		{
			z.push_back(temp);
			cou++;
			temp = start_z + cou * zspace;
		} while ((temp < end_z) || abs(temp - end_z) < EPS);

		num = y.size() * x.size() * z.size();
		// os << num << endl;
		for (iter iz = z.begin(); iz != z.end(); ++iz)
		{
			for (iter iy = y.begin(); iy != y.end(); ++iy)
			{
				for (iter ix = x.begin(); ix != x.end(); ++ix)
				{
					os << setw(15) << left << (*ix)
					   << setw(15) << left << (*iy)
					   << setw(15) << left << (*iz)
					   << '\n';
				}
			}
		}
	}
	catch (std::exception &e)
	{
		cout << e.what() << "\n";
		return 1;
	}
	return 0;
}
