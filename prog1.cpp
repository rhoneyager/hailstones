#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <vtkImageData.h>
#include <vtkImageDataToPointSet.h>
#include <vtkImageStencil.h>
#include <vtkMassProperties.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataToImageStencil.h>
#include <vtkSmartPointer.h>
#include <vtkSTLReader.h>
#include <vtkStructuredGrid.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTriangleFilter.h>
#include <vtkXMLImageDataWriter.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char** argv) {
	using namespace std;
	try {
		namespace po = boost::program_options;
		po::options_description desc("Allowed options");
		desc.add_options()
			("help,h", "Produce this help message.")
			("resolution,r", po::value<double>()->default_value(1.0), "Resolution of output mesh lattices (mm).")
			("input,i", po::value<vector<string> >()->multitoken(), "Input stl files.")
			("output,o", po::value<string>()->default_value("."), "Output folder.")
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

		if (vm.count("help")) doHelp();
		if (!vm.count("input")) doHelp();
		bool clobber = false;
		if (vm.count("clobber")) clobber = true;
		string sOutFolder = vm["output"].as<string>();
		using namespace boost::filesystem;
		path pOutFolder(sOutFolder);
		vector<string> vInputs = vm["input"].as<vector<string> >();
		{
			while (is_symlink(pOutFolder)) pOutFolder = read_symlink(pOutFolder);
			if (!exists(pOutFolder)) {
				if (!create_directory(pOutFolder)) {
					cerr << "Cannot create output directory at " << pOutFolder << endl;
					exit(3);
				}
			} else if (!is_directory(pOutFolder)) {
				cerr << "Output path " << pOutFolder << " exists and is not a folder." << endl;
				exit(3);
			}

			for (const auto &i : vInputs) {
				path p(i);
				while (is_symlink(p)) p = read_symlink(p);
				if (!exists(p)) {
					cerr << "Input file " << p << " does not exist." << endl;
					exit(3);
				}
				if (!is_regular_file(p)) {
					cerr << "Input file " << p << " is not a regular file." << endl;
					exit(3);
				}
				if (!p.has_extension()) {
					cerr << "Input file " << p << " is not an stl file." << endl;
					exit(3);
				}
				if (p.extension().string() != ".stl") {
					cerr << "Input file " << p << " is not an stl file." << endl;
					exit(3);
				}
			}
		}
		// Input is validated. Now, to loop over each file and process.
		double resolution_mm = vm["resolution"].as<double>();
		cout << "Converting files at resolution of " << resolution_mm << " mm." << endl;

		const double den_ice_g_cm3 = 0.9167;
		struct analysis_data {
			double contouredVol_cm3, contouredMass_g, discretized_volume_cm3;
			double actual_bounds[6], mesh_bounds[6], actual_bound_volume_cm3;
			double vf_contoured_boxed, resolution_mm, aeff_um;
			size_t numIceLatticeSites, numLatticeSites;
			std::string sIn, sOut, sFileId;
			static void writeHeader(std::ostream &out) {
				out << "ID,contouredVol_cm3,contouredMass_g,discretized_volume_cm3,"
					"actual_bounds[x_min](mm),actual_bounds[x_max](mm),actual_bounds[y_min](mm),"
					"actual_bounds[y_max](mm),actual_bounds[z_min](mm),actual_bounds[z_max](mm),"
					"mesh_bounds[x_min](mm),mesh_bounds[x_max](mm),mesh_bounds[y_min](mm),"
					"mesh_bounds[y_max](mm),mesh_bounds[z_min](mm),mesh_bounds[z_max](mm),"
					"actual_bound_volume_cm3,vf_contoured_boxed,resolution_mm,aeff_um,"
					"numLatticeSites,numIceLatticeSites" << endl;
			}
			void write(std::ostream &out) {
				out << sFileId << "," << contouredVol_cm3 << "," << contouredMass_g << "," << discretized_volume_cm3 << ","
					<< actual_bounds[0] << "," << actual_bounds[1] << "," << actual_bounds[2] << ","
					<< actual_bounds[3] << "," << actual_bounds[4] << "," << actual_bounds[5] << ","
					<< mesh_bounds[0] << "," << mesh_bounds[1] << "," << mesh_bounds[2] << ","
					<< mesh_bounds[3] << "," << mesh_bounds[4] << "," << mesh_bounds[5] << ","
					<< actual_bound_volume_cm3 << "," << vf_contoured_boxed << "," << resolution_mm << ","
					<< aeff_um << "," << numLatticeSites << "," << numIceLatticeSites << endl;
			}
		};
		path pOutAnalysis = pOutFolder / path("analysis.csv");
		if (exists(pOutAnalysis) && !clobber) {
			cout << "\t\tAnalysis file already exists." << endl;
			exit(4);
		}
		std::ofstream canal(pOutAnalysis.c_str());
		analysis_data::writeHeader(canal);
		analysis_data data;

		for (const auto &i : vInputs) {
			path p(i);
			path pOut = pOutFolder / p.filename().replace_extension(path(".vti"));
			cout << "\t" << p << " to " << pOut << endl;
			if (exists(pOut) && !clobber) {
				cout << "\t\tFile already exists. Skipping." << endl;
				continue;
			}
			data.sIn = i;
			data.sOut = pOut.string();
			data.sFileId = p.filename().replace_extension().string();

			vtkSmartPointer<vtkSTLReader> reader =
				vtkSmartPointer<vtkSTLReader>::New();
			reader->SetFileName(i.c_str());
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
			translation->Translate(2.-data.actual_bounds[0], 2. - data.actual_bounds[2], 2. - data.actual_bounds[4]);
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
				data.vf_contoured_boxed = data.contouredVol_cm3 / data.actual_bound_volume_cm3;

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
					dim[i] = static_cast<int>(ceil((data.actual_bounds[i * 2 + 1] - data.actual_bounds[i * 2]) / spacing[i]));
				}
				whiteImage->SetDimensions(dim);
				whiteImage->SetExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1);

				double origin[3];
				origin[0] = data.actual_bounds[0] + spacing[0] / 2;
				origin[1] = data.actual_bounds[2] + spacing[1] / 2;
				origin[2] = data.actual_bounds[4] + spacing[2] / 2;
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
				//int* dims = imgd->GetDimensions();
				// int dims[3]; // can't do this
				data.numIceLatticeSites = 0;
				//std::cout << "Dims: " << " x: " << dims[0] << " y: " << dims[1] << " z: " << dims[2] << std::endl;

				//std::cout << "Number of points: " << imgd->GetNumberOfPoints() << std::endl;
				//std::cout << "Number of cells: " << imgd->GetNumberOfCells() << std::endl;

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
		

	} catch (std::exception &e) {
		cerr << e.what() << endl;
		return 1;
	}
	return 0;
}