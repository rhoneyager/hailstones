#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
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
#include "refract.h"

/// Program 3b takes vti shape files and generates DDSCAT runs based on the provided temperature and frequency.

int main(int argc, char** argv) {
	using namespace std;
	try {
		namespace po = boost::program_options;
		po::options_description desc("Allowed options");
		desc.add_options()
			("help,h", "Produce this help message.")
			("resolution,r", po::value<string>()->default_value("0.6"), "Resolution of output mesh lattices (mm).")
			("input,i", po::value<vector<string> >()->multitoken(), "Input vti files.")
			("output,o", po::value<string>()->default_value("."), "Output base folder for the DDSCAT runs. Each run is placed in a uniquely-named subfolder.")
			("clobber,c", "Overwrite existing files.")
			("temperature,t", po::value<string>()->default_value("270"), "Temperature (in K) for refractive index calculations.")
			("frequency,f", po::value<string>()->default_value("13.6"), "Frequency, in GHz")
			("num-betas", po::value<unsigned short>()->default_value(12), "Number of orientation in Beta direction.")
			("num-thetas", po::value<unsigned short>()->default_value(13), "Number of orientation in Theta direction.")
			("num-phis", po::value<unsigned short>()->default_value(12), "Number of orientation in Phi direction.")
			;
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).
			options(desc).run(), vm);
		po::notify(vm);
		auto doHelp = [&]()
		{
			std::cerr << "This program takes a volumetric shape file (.vti), a frequency (GHz) and "
				"a temperature (K), and generates the files needed to start a DDSCAT run. DDSCAT is "
				"used to deterine the radiative properties of the particle." << std::endl;
			std::cerr << desc << std::endl;
			exit(2);
		};

		if (vm.count("help")) doHelp();
		if (!vm.count("input")) doHelp();
		bool clobber = false;
		if (vm.count("clobber")) clobber = true;
		string sOutFolder = vm["output"].as<string>();
		unsigned short nb = vm["num-betas"].as<unsigned short>();
		unsigned short nt = vm["num-thetas"].as<unsigned short>();
		unsigned short np = vm["num-phis"].as<unsigned short>();
		string stempK = vm["temperature"].as<string>();
		string sfreqGHz = vm["frequency"].as<string>();
		string sres_mm = vm["resolution"].as<string>();
		double tempK = boost::lexical_cast<double>(stempK), freqGHz = boost::lexical_cast<double>(sfreqGHz),
			res_mm = boost::lexical_cast<double>(sres_mm);
		complex<double> m;
		mIceMatzler(freqGHz, tempK, m);
		const double c = 2.9979e8;
		const double wvlen_mm = c / (freqGHz * 1000 * 1000);
		const double pi = 3.141592654;

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
			string sRun = p.filename().replace_extension(path("")).string();
			sRun = sRun + "_" + sfreqGHz + "_" + stempK
				+ "_" + boost::lexical_cast<string>(nb)
				+ "_" + boost::lexical_cast<string>(nt)
				+ "_" + boost::lexical_cast<string>(np);
			path pOut = pOutFolder / path(sRun);
			path pOutShp = pOut / path("shape.dat");
			path pOutDiel = pOut / path("diel.tab");
			path pOutPar = pOut / path("ddscat.par");
			cout << p << " to " << pOut << endl;
			if (exists(pOut)) {
				path pSym = pOut;
				while (is_symlink(pSym)) pSym = read_symlink(pSym);
				if (exists(pSym)) {
					if (is_directory(pSym)) {
						if (!clobber) {
							cerr << "Directory " << pSym << " already exists, and not clobbering." << endl;
							exit(4);
						}
					}
				} else {
					bool res = boost::filesystem::create_directory(pSym);
					if (!res) {
						cerr << "Could not create " << pSym << endl;
						exit(4);
					}
				}
			} else {
				bool res = boost::filesystem::create_directory(pOut);
				if (!res) {
					cerr << "Could not create " << pOut << endl;
					exit(4);
				}
			}

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
						if (pixel > 0.1) {
							numIceLatticeSites++;
							c[0] += x;
							c[1] += y;
							c[2] += z;
						}
					}
				}
			}
			c[0] /= (double)numIceLatticeSites;
			c[1] /= (double)numIceLatticeSites;
			c[2] /= (double)numIceLatticeSites;

			double V_mm3 = std::pow<double,double>(res_mm, 3.) * (double) numIceLatticeSites;
			double aeff_mm = std::pow<double,double>(3.*V_mm3 / (4. * pi), 1. / 3.);

			writeShp(pOutShp.string().c_str(), 
				p.filename().replace_extension(path("")).string().c_str(), 
				imgd, c, numIceLatticeSites, sres_mm.c_str());
			writeDiel(pOutDiel.string().c_str(), m, sfreqGHz.c_str(), stempK.c_str());
			writePar(pOutPar.string().c_str(), wvlen_mm, aeff_mm, nb, nt, np, dims);
		}

	}
	catch (std::exception &e) {
		cerr << e.what() << endl;
		return 1;
	}
	return 0;
}