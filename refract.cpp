#include <complex>

#include "refract.h"

/// Needs frequency in GHz and temp in K
/// Valid for 0 to 1000 GHz, and and ice temp
void mIceMatzler(double f, double t, std::complex<double> &m)
{
	double er = 0;
	if (t>243.0)
		er = 3.1884 + 9.1e-4*(t - 273.0);
	else
		er = 3.1611 + 4.3e-4*(t - 243.0);
	// Imaginary part
	double theta = 300.0 / (t)-1.0;
	double alpha = (0.00504 + 0.0062*theta)*exp(-22.1*theta);
	double dbeta = exp(-9.963 + 0.0372*(t - 273.16));
	const double B1 = 0.0207;
	const double B2 = 1.16e-11;
	const double b = 335;
	double betam = B1 / t*exp(b / t) / pow((exp(b / t) - 1.0), 2) + B2*pow(f, 2.0);
	double beta = betam + dbeta;
	double ei = alpha / f + beta*f;
	std::complex<double> e(er, -ei);
	m = sqrt(e);
}

void writeShp(const char* file, const char* desc, vtkSmartPointer<vtkImageData> imgd, double c[3], size_t numIceLatticeSites) {
	// Open the output file, write the header and then write the points.
	std::ofstream out(file);
	if (desc)
		out << desc << endl;
	else throw;
	out << "\t" << numIceLatticeSites << "\t= NAT" << endl;
	out << "\t1.0\t0.0\t0.0 =\tvector" << endl
		<< "\t0.0\t1.0\t0.0 =\tvector" << endl
		<< "\t1.0\t1.0\t1.0 =\tlattice spacings (d_x,d_y,d_z)/d" << endl
		<< "\t" << c[0] << "\t" << c[1] << "\t" << c[2] << " =\tlattice offset" << endl
		<< "\tJA\tIX\tIY\tIZ\tICOMP(x, y, z)" << endl;
	size_t i = 0;
	int* dims = imgd->GetDimensions();
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

void writeDiel(const char* file, std::complex<double>& r) {
	using namespace std;
	ofstream out(file);
	//out.setf( ios::scientific, ios::floatfield);
	//out.precision(7);
	//out.unsetf(ios_base::floatfield);
	double rr = r.real();
	double ri = (-1.0 * abs(r.imag()));
	out << " m = " << rr << " + " << ri << " i" << endl;
	out << " 0 1 2" << endl;
	out << " = columns for wave, Re(m), Im(m)" << endl;
	out << " 0.000001    " << rr << "      " << ri << endl;
	out << " 1.000000    " << rr << "      " << ri << endl;
	out << " 100000.0    " << rr << "      " << ri << endl;
}

void writePar(const char* file, double wvlen_mm, double aeff_mm, unsigned short nb, unsigned short nt, unsigned short np) {
	throw;
}
