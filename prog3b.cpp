#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLImageDataWriter.h>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

int main(int argc, char** argv) {
	using namespace std;
	try {
		namespace po = boost::program_options;
		po::options_description desc("Allowed options");
		desc.add_options()
			("help,h", "Produce this help message.")
			("resolution,r", po::value<double>()->default_value(0.6), "Resolution of output mesh lattices (mm).")
			("input,i", po::value<vector<string> >()->multitoken(), "Input vti files.")
			("output,o", po::value<string>()->default_value("."), "Output folder for the DDSCAT shape files.")
			("clobber,c", "Overwrite existing files.")
			;
		po::positional_options_description pos;
		pos.add("input", -1);
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).
			options(desc).positional(pos).run(), vm);
		po::notify(vm);
		auto doHelp = [&]()
		{
			std::cerr << desc << std::endl;
			exit(2);
		};
		double resolution_mm = vm["resolution"].as<double>();

		if (vm.count("help")) doHelp();
		if (!vm.count("input")) doHelp();
		bool clobber = false;
		if (vm.count("clobber")) clobber = true;
		string sOutFolder = vm["output"].as<string>();
		using namespace boost::filesystem;
		path pOutFolder(sOutFolder);
		vector<string> vInputs = vm["input"].as<vector<string> >();
		vector<path> veInputs;
		{ // Validate and expand inputs.
			while (is_symlink(pOutFolder)) pOutFolder = read_symlink(pOutFolder);
			if (!exists(pOutFolder)) {
				if (!create_directory(pOutFolder)) {
					cerr << "Cannot create output directory at " << pOutFolder << endl;
					exit(3);
				}
			}
			else if (!is_directory(pOutFolder)) {
				cerr << "Output path " << pOutFolder << " exists and is not a folder." << endl;
				exit(3);
			}
			for (const auto &i : vInputs) {
				path p(i);
				while (is_symlink(p)) p = read_symlink(p);
				if (!exists(p)) {
					cerr << "Input file or directory " << p << " does not exist." << endl;
					exit(3);
				}
				if (is_directory(p)) {
					// Iterate over one level of directory entries to get files.
					// Useful when there are too many files for shell expansion, or on Windows.
					vector<path> dircontents;
					copy(directory_iterator(p), directory_iterator(),
						back_inserter(dircontents));
					for (const auto &di : dircontents) {
						path dp(di);
						while (is_symlink(dp)) dp = read_symlink(dp);
						if (!exists(dp)) continue;
						if (!is_regular_file(dp)) continue;
						if (!dp.has_extension()) continue;
						if (dp.extension().string() == ".vti") veInputs.push_back(dp);
					}
				}
				else if (is_regular_file(p)) {
					if (!p.has_extension()) {
						cerr << "Input file " << p << " is not an vti file." << endl;
						exit(3);
					}
					if (p.extension().string() != ".vti") {
						cerr << "Input file " << p << " is not a vti file." << endl;
						exit(3);
					}
					veInputs.push_back(p);
				}
				else {
					cerr << "Input file " << p << " is not a regular file." << endl;
					exit(3);
				}
			}
		}
		// Input is validated. Now, to loop over each file and process.

		for (const boost::filesystem::path &p : veInputs) {
			path pOut = pOutFolder / p.filename().replace_extension(path(".shp"));
			cout << p << " to " << pOut << endl;

			vtkSmartPointer<vtkXMLImageDataReader> reader =
				vtkSmartPointer<vtkXMLImageDataReader>::New();
			reader->SetFileName(p.string().c_str());
			reader->Update();

			// Get number of filled lattice sites. Calculate the center of mass.
			double c[3] = { 0, 0, 0 };

			vtkSmartPointer<vtkImageData> imgd = reader->GetOutput();
			size_t numIceLatticeSites = 0;
			int* dims = imgd->GetDimensions();
			for (int z = 0; z < dims[2]; z++)
			{
				for (int y = 0; y < dims[1]; y++)
				{
					for (int x = 0; x < dims[0]; x++)
					{
						double pixel = (imgd->GetScalarComponentAsDouble(x, y, z, 0));
						if (pixel > 0.1) numIceLatticeSites++;
						c[0] += x;
						c[1] += y;
						c[2] += z;
					}
				}
			}
			c[0] /= (double)numIceLatticeSites;
			c[1] /= (double)numIceLatticeSites;
			c[2] /= (double)numIceLatticeSites;

			// Open the output file, write the header and then write the points.
			std::ofstream out(pOut.string().c_str());
			out << "From " << p.filename() << endl;
			out << "\t" << numIceLatticeSites << "\t= NAT" << endl;
			out << "\t1.0\t0.0\t0.0 =\tvector" << endl
				<< "\t0.0\t1.0\t0.0 =\tvector" << endl
				<< "\t1.0\t1.0\t1.0 =\tlattice spacings (d_x,d_y,d_z)/d" << endl
				<< "\t" << c[0] << "\t" << c[1] << "\t" << c[2] << " =\tlattice offset" << endl
				<< "\tJA\tIX\tIY\tIZ\tICOMP(x, y, z)" << endl;
			size_t i = 0;
			for (int z = 0; z < dims[2]; z++)
			{
				for (int y = 0; y < dims[1]; y++)
				{
					for (int x = 0; x < dims[0]; x++)
					{
						double pixel = (imgd->GetScalarComponentAsDouble(x, y, z, 0));
						if (pixel > 0.1) {
							++i;
							out << "\t" << i << "\t" << x << "\t" << y << "\t" << z << "\t\t1\t1\t1\n";
						}
					}
				}
			}
		}

	}
	catch (std::exception &e) {
		cerr << e.what() << endl;
		return 1;
	}
	return 0;
}