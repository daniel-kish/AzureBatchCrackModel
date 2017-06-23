#include <vector>
#include <fstream>
#include <utility>
#include <string>
#include <sstream>
#include <iostream>
#include <exception>
using namespace std;


void reSize(string in, string out)
{
	ifstream from( in );

	if (!from)
	{
		cout << "from error " << in << '\n';
		throw std::exception(in.c_str());
	}

	ofstream to( out );
	if (!to)
	{
		cout << "to error\n" << out << '\n';
		throw std::exception(out.c_str());
	}

	double a{}, b{};
	int N{ 30 };
	int i{ N };
	int total{ 0 };

	while (from >> a >> b) 
	{
		//if (total < 15) {
		//	to << b*6250 << ' ' << a * 0.006324 << '\n'; // L*6250 T*0.0006324
		//	total++;
		//}
		//else
		//{
		//	if (i == 0) {
		//		to << b*0.001238 << ' ' << a * 2026 << '\n';
		//		i = N;
		//	}
		//	else
		//		i--;
		//}
		to << a * 6250 << ' ' << b * 0.000275 << '\n';
	}
}

int main(int argc, char*argv[])
{
	string fromDir = argv[1];
	string toDir = argv[2];
	int nr_files = stoi(argv[3]);
	ostringstream fromFilename;
	ostringstream toFilename;
	int i{1};

	while (i != nr_files)
	{
		ostringstream fromFilename;
		ostringstream toFilename;

		fromFilename << fromDir << '\\' << "out_" << i << ".txt";
		toFilename << toDir << '\\' << "out_" << i << ".txt";
		reSize( fromFilename.str(), toFilename.str() );

		fromFilename.clear(); cout << fromFilename.str() << '\n';
		toFilename.clear(); cout << toFilename.str() << '\n';
		i++;
	}
}