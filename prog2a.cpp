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

void readMasses(const boost::filesystem::path &p, std::map<std::string, double> &res) {
	using namespace std;
	ifstream in(p.c_str());
	while (in.good()) {
		string lin;
		getline(in, lin);
		if (!lin.size()) continue;
		if (lin.at(0) == '\r' || lin.at(0) == '\n' || lin.at(0) == '#') continue;
		if (lin.find("ID") == 0) continue;
		// Split according to the presence of a comma. Reading only the first two columns of each line.
		size_t copos = lin.find_first_of(',');
		if (copos == string::npos) continue; // Bad line
		if (copos < 2) continue; // Bad line
		if (copos >= ((lin.length() > 2) ? lin.length() - 2 : 0)) continue;
		string id = lin.substr(0, copos);
		size_t endnum = lin.find_first_not_of("0123456789.", copos+1);
		string smass = lin.substr(copos+1,endnum);
		istringstream iconv(smass);
		double val = 0;
		iconv >> val;

		// Reformat ID to be in lower case, and to stop after the first space
		string idt = id.substr(0, id.find_first_of(' '));
		string idtl;
		idtl.resize(idt.length());
		transform(idt.begin(), idt.end(), idtl.begin(), tolower);
		
		res[idtl] = val;
	}
}

int main(int argc, char** argv) {
	using namespace std;
	try {
		namespace po = boost::program_options;
		po::options_description desc("Allowed options");
		desc.add_options()
			("help,h", "Produce this help message.")
			("masses,m", po::value<string>(), "Path to csv file containing shape id,mass(g) of each particle")
			("resolution,r", po::value<double>()->default_value(0.6), "Resolution of output mesh lattices (mm).")
			("input,i", po::value<vector<string> >()->multitoken(), "Input vti files.")
			("output,o", po::value<string>()->default_value("."), "Output folder for the decimated hailstones.")
			("clobber,c", "Overwrite existing files.")
			;
		po::positional_options_description p;
		p.add("input", -1);
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).
			options(desc).positional(p).run(), vm);
		po::notify(vm);
		auto doHelp = [&]()
		{
			std::cerr << desc << std::endl;
			exit(2);
		};
		double resolution_mm = vm["resolution"].as<double>();
		unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
		mt19937 rndgen(seed);  // std::mt19937 is a standard mersenne_twister_engine

		if (vm.count("help")) doHelp();
		if (!vm.count("input")) doHelp();
		if (!vm.count("masses")) doHelp();
		bool clobber = false;
		if (vm.count("clobber")) clobber = true;
		string sOutFolder = vm["output"].as<string>();
		using namespace boost::filesystem;
		path pOutFolder(sOutFolder);
		string sMasses = vm["masses"].as<string>();
		path pMasses(sMasses);
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
			while (is_symlink(pMasses)) pMasses = read_symlink(pMasses);
			if (!exists(pMasses)) {
				cerr << "Mass information file does not exist: " << pMasses << endl;
				exit(3);
			}
			else if (!is_regular_file(pMasses)) {
				cerr << "Mass information file " << pMasses << " exists and is not a regular file." << endl;
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
						cerr << "Input file " << p << " is not an stl file." << endl;
						exit(3);
					}
					if (p.extension().string() != ".stl") {
						cerr << "Input file " << p << " is not an stl file." << endl;
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
		
		const double den_ice_g_cm3 = 0.9167;
		map<string, double> m_ids_masses;
		readMasses(pMasses, m_ids_masses);
		if (m_ids_masses.size() != veInputs.size())
			cout << "Warning: veInputs and m_ids_masses do not have the same number of entries." << endl
			<< "veInputs has " << veInputs.size() << " file"
			<< (veInputs.size() == 1 ? "" : "s")
			<< ", and m_ids_masses has " << m_ids_masses.size()
			<< (m_ids_masses.size() == 1 ? " entry." : " entries.") << endl;
		// Take file names and make sure that they all match an entry in the m_ids_masses.

		map<path, double> mFileMasses;
		for (const auto &p : veInputs) {
			string preid = p.filename().replace_extension("").string();
			string preidt = preid.substr(0, preid.find_first_of(' '));
			string id;
			id.resize(preidt.length());
			transform(preidt.begin(), preidt.end(), id.begin(), tolower);

			if (!m_ids_masses.count(id)) {
				cerr << "Cannot find a mass for the file " << p << endl;
				continue; // exit(5);
			}
			mFileMasses[p] = m_ids_masses.at(id);
		}

		struct analysis_data {
			double contouredMass_g, actualMass_g, volumeFraction;
			std::string sFileId;
			static void writeHeader(std::ostream &out) {
				out << "ID,contouredMass_g,actualMass_g,volumeFraction" << endl;
			}
			void write(std::ostream &out) {
				out << sFileId << "," << contouredMass_g << "," << actualMass_g << "," << volumeFraction << endl;
			}
		};

		path pOutAnalysis = pOutFolder / path("analysis-2a.csv");
		if (exists(pOutAnalysis) && !clobber) {
			cout << "\t\tAnalysis file already exists." << endl;
			exit(4);
		}
		std::ofstream canal(pOutAnalysis.c_str());
		analysis_data::writeHeader(canal);
		analysis_data data;


		for (const auto &i : mFileMasses) {
			double mass_g = i.second;
			path p = i.first;
			path pOut = pOutFolder / p.filename().replace_extension(path(".vti"));
			cout << p << " to " << pOut << endl;
			if (mass_g < 0) {
				cout << "\tSkipping as there is no valid mass." << endl;
				continue;
			}


			vtkSmartPointer<vtkXMLImageDataReader> reader =
				vtkSmartPointer<vtkXMLImageDataReader>::New();
			reader->SetFileName(p.string().c_str());
			reader->Update();

			// Get number of filled sites, and calculate filled mass.

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
					}
				}
			}

			// Compare filled mass to actual mass to get the retention fraction
			data.contouredMass_g = (double)(numIceLatticeSites) * pow(resolution_mm/10.,3.) * den_ice_g_cm3;
			data.volumeFraction = data.actualMass_g / data.contouredMass_g;
			data.actualMass_g = mass_g;
			cout << "\tmass(g): " << mass_g << "\tfilled mass (g): " << data.contouredMass_g 
				<< "\tRetention fraction: " << data.volumeFraction << endl;
			if (data.volumeFraction > 1.) data.volumeFraction = 1.;

			// Iterate over points and randomly decimate
			for (int z = 0; z < dims[2]; z++)
			{
				for (int y = 0; y < dims[1]; y++)
				{
					for (int x = 0; x < dims[0]; x++)
					{
						double pixel = (imgd->GetScalarComponentAsDouble(x, y, z, 0));
						if (pixel > 0.1) {
							double rndval = (double)rndgen() / (double)pow(2., rndgen.word_size);
							if (rndval > data.volumeFraction) {
								imgd->SetScalarComponentFromDouble(x, y, z, 0, 0.);
							}
						}
					}
				}
			}

			// Write output file
			vtkSmartPointer<vtkXMLImageDataWriter> writer =
				vtkSmartPointer<vtkXMLImageDataWriter>::New();
			writer->SetInputData(imgd);
			writer->SetFileName(pOut.string().c_str());
			writer->Write();
			
			data.write(canal);
		}

		/*
		path pOutAnalysis = pOutFolder / path("analysis.csv");
		if (exists(pOutAnalysis) && !clobber) {
			cout << "\t\tAnalysis file already exists." << endl;
			exit(4);
		}
		std::ofstream canal(pOutAnalysis.c_str());
		analysis_data::writeHeader(canal);
		analysis_data data;

		for (const auto &p : veInputs) {
			path pOut = pOutFolder / p.filename().replace_extension(path(".vti"));
			cout << "\t" << p << " to " << pOut << endl;
			if (exists(pOut) && !clobber) {
				cout << "\t\tFile already exists. Skipping." << endl;
				continue;
			}
			data.sIn = p.string();
			data.sOut = pOut.string();
			data.sFileId = p.filename().replace_extension().string();

			vtkSmartPointer<vtkSTLReader> reader =
				vtkSmartPointer<vtkSTLReader>::New();
			reader->SetFileName(p.string().c_str());
			reader->Update();

			vtkSmartPointer<vtkPolyData> pd = reader->GetOutput();
			pd->GetBounds(data.actual_bounds);

			// Determine contoured mesh volume, solid ice mass and volume fraction

			vtkSmartPointer<vtkTriangleFilter> triangleFilter =
				vtkSmartPointer<vtkTriangleFilter>::New();
			triangleFilter->SetInputConnection(reader->GetOutputPort());
			triangleFilter->Update();

			vtkSmartPointer<vtkMassProperties> massProps =
				vtkSmartPointer<vtkMassProperties>::New();
			massProps->SetInputConnection(triangleFilter->GetOutputPort());
			massProps->Update();

			data.contouredVol_cm3 = massProps->GetVolume();
			data.contouredVol_cm3 /= std::pow(10., 3.);
			data.contouredMass_g = data.contouredVol_cm3 / den_ice_g_cm3;

			// Translate mesh so that the min point is located at (2,2,2). Needed for ddscat.

			vtkSmartPointer<vtkTransform> translation =
				vtkSmartPointer<vtkTransform>::New();
			translation->Translate(2. - data.actual_bounds[0], 2. - data.actual_bounds[2], 2. - data.actual_bounds[4]);
			vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter =
				vtkSmartPointer<vtkTransformPolyDataFilter>::New();
			transformFilter->SetInputConnection(reader->GetOutputPort());
			transformFilter->SetTransform(translation);
			transformFilter->Update();
			vtkSmartPointer<vtkPolyData> pt = transformFilter->GetOutput();
			pt->GetBounds(data.actual_bounds);

			data.resolution_mm = resolution_mm;

			// Convert contoured surface to volume representation
			{
				// These initial bounds are in mm
				data.mesh_bounds[0] = data.actual_bounds[0] - 1;
				data.mesh_bounds[1] = data.actual_bounds[1] + 1;
				data.mesh_bounds[2] = data.actual_bounds[2] - 1;
				data.mesh_bounds[3] = data.actual_bounds[3] + 1;
				data.mesh_bounds[4] = data.actual_bounds[4] - 1;
				data.mesh_bounds[5] = data.actual_bounds[5] + 1;
				data.actual_bound_volume_cm3 = std::pow(1. / 10., 3.)
					* (data.actual_bounds[1] - data.actual_bounds[0])
					* (data.actual_bounds[3] - data.actual_bounds[2])
					* (data.actual_bounds[5] - data.actual_bounds[4]);
				data.mesh_bound_volume_cm3 = std::pow(1. / 10., 3.)
					* (data.mesh_bounds[1] - data.mesh_bounds[0])
					* (data.mesh_bounds[3] - data.mesh_bounds[2])
					* (data.mesh_bounds[5] - data.mesh_bounds[4]);
				data.vf_contoured_boxed = data.contouredVol_cm3 / data.mesh_bound_volume_cm3;

				double spacing[3]; // desired volume spacing
				spacing[0] = resolution_mm;
				spacing[1] = resolution_mm;
				spacing[2] = resolution_mm;
				vtkSmartPointer<vtkImageData> whiteImage =
					vtkSmartPointer<vtkImageData>::New();
				whiteImage->SetSpacing(spacing);

				// compute dimensions
				int dim[3];
				for (int i = 0; i < 3; i++)
				{
					dim[i] = static_cast<int>(ceil((data.mesh_bounds[i * 2 + 1] - data.mesh_bounds[i * 2]) / spacing[i]));
				}
				whiteImage->SetDimensions(dim);
				whiteImage->SetExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1);

				double origin[3];
				origin[0] = data.mesh_bounds[0] + spacing[0] / 2;
				origin[1] = data.mesh_bounds[2] + spacing[1] / 2;
				origin[2] = data.mesh_bounds[4] + spacing[2] / 2;
				whiteImage->SetOrigin(origin);

				whiteImage->AllocateScalars(VTK_UNSIGNED_CHAR, 1);

				// fill the image with foreground voxels:
				unsigned char inval = 1;
				unsigned char outval = 0;
				vtkIdType count = whiteImage->GetNumberOfPoints();
				for (vtkIdType i = 0; i < count; ++i)
				{
					whiteImage->GetPointData()->GetScalars()->SetTuple1(i, inval);
				}
				data.numLatticeSites = (size_t)count;

				// polygonal data --> image stencil:
				vtkSmartPointer<vtkPolyDataToImageStencil> pol2stenc =
					vtkSmartPointer<vtkPolyDataToImageStencil>::New();
				pol2stenc->SetInputData(pt);
				pol2stenc->SetOutputOrigin(origin);
				pol2stenc->SetOutputSpacing(spacing);
				pol2stenc->SetOutputWholeExtent(whiteImage->GetExtent());
				pol2stenc->Update();

				// cut the corresponding white image and set the background:
				vtkSmartPointer<vtkImageStencil> imgstenc =
					vtkSmartPointer<vtkImageStencil>::New();
				imgstenc->SetInputData(whiteImage);
				imgstenc->SetStencilConnection(pol2stenc->GetOutputPort());
				imgstenc->ReverseStencilOff();
				imgstenc->SetBackgroundValue(outval);
				imgstenc->Update();

				// Determine the number of filled points on the grid
				vtkSmartPointer<vtkImageData> imgd = imgstenc->GetOutput();
				data.numIceLatticeSites = 0;
				for (int z = 0; z < dim[2]; z++)
				{
					for (int y = 0; y < dim[1]; y++)
					{
						for (int x = 0; x < dim[0]; x++)
						{
							double pixel = (imgd->GetScalarComponentAsDouble(x, y, z, 0));
							if (pixel > 0.1) data.numIceLatticeSites++;
						}
					}
				}

				vtkSmartPointer<vtkXMLImageDataWriter> writer =
					vtkSmartPointer<vtkXMLImageDataWriter>::New();
				writer->SetInputConnection(imgstenc->GetOutputPort());
				writer->SetFileName(pOut.string().c_str());
				writer->Write();
			}

			// V = Nd^3
			// Calculate the discretized volume and effective radius
			data.discretized_volume_cm3 = data.numIceLatticeSites * std::pow(data.resolution_mm / 10., 3.);
			//double dv_um3 = data.discretized_volume_cm3 * std::pow(10000., 3);
			const double pi = 3.141592654;
			double aeff_cm = 3. * data.discretized_volume_cm3 / 4. * pi;
			data.aeff_um = aeff_cm * 10000.;
			data.write(canal);
		}
		*/
		


	}
	catch (std::exception &e) {
		cerr << e.what() << endl;
		return 1;
	}
	return 0;
}