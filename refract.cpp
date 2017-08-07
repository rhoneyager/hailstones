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

void writeShp(const char* file, const char* desc, 
	vtkSmartPointer<vtkImageData> imgd, double c[3], size_t numIceLatticeSites,
	const char* res
) {
	// Open the output file, write the header and then write the points.
	std::ofstream out(file);
	if (desc)
		out << desc << ", at resolution of " << res << " mm" << endl;
	else throw;
	out << "\t" << numIceLatticeSites << "\t= NAT" << endl;
	out << "\t1.0\t0.0\t0.0 =\tvector" << endl
		<< "\t0.0\t1.0\t0.0 =\tvector" << endl
		<< "\t1.0\t1.0\t1.0 =\tlattice spacings (d_x,d_y,d_z)/d" << endl
		<< "\t" << c[0] << "\t" << c[1] << "\t" << c[2] << " =\tlattice offset" << endl
		<< "JA\tIX\tIY\tIZ\tICOMP(x, y, z)" << endl;
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
					out << i << "\t" << x << "\t" << y << "\t" << z << "\t\t1\t1\t1\n";
				}
			}
		}
	}
}

void writeDiel(const char* file, std::complex<double>& r, const char* f, const char* t) {
	using namespace std;
	ofstream out(file);
	//out.setf( ios::scientific, ios::floatfield);
	//out.precision(7);
	//out.unsetf(ios_base::floatfield);
	double rr = r.real();
	double ri = (1.0 * abs(r.imag()));
	out << " Diels for f = " << f << " GHz, T = " << t << " K (m = " << rr << " + " << ri << " i)" << endl;
	out << " 1 2 3 0 0 = columns for wave, Re(m), Im(m)" << endl;
	out << " LAMBDA\tRe(M)\tIm(M)" << endl;
	out << " 0.000001    " << rr << "      " << ri << endl;
	out << " 1.000000    " << rr << "      " << ri << endl;
	out << " 100000.0    " << rr << "      " << ri << endl;
}

void writePar(const char* file, double wvlen_mm, double aeff_mm, unsigned short nb, unsigned short nt, unsigned short np, const int* ndims) {
	using namespace std;
	ofstream out(file);
	out <<
		"' ========= Parameter file for v7.3 =================== '\n"
		"'**** Preliminaries ****'\n"
		"'NOTORQ' = CMTORQ * 6 (NOTORQ, DOTORQ) --either do or skip torque calculations\n"
		"'PBCGS2' = CMDSOL * 6 (PBCGS2, PBCGST, PETRKP) --select solution method\n"
		"'GPFAFT' = CMDFFT * 6 (GPFAFT, FFTMKL)-- - FFT method\n"
		"'GKDLDR' = CALPHA * 6 (GKDLDR, LATTDR)\n"
		"'NOTBIN' = CBINFLAG(NOTBIN, ORIBIN, ALLBIN)\n"
		"'**** Initial Memory Allocation ****'\n"
		<< ndims[0] + 2 << " " << ndims[1] + 2 << " " << ndims[2] + 2 << " = dimension\n"
		"'**** Target Geometry and Composition ****'\n"
		"'FROM_FILE' = CSHAPE * 9 shape directive\n"
		"50 50 924 = shape parameters 1 - 3\n"
		"1 = NCOMP = number of dielectric materials\n"
		"'diel.tab' = file with refractive index 1\n"
		"'**** Additional Nearfield calculation? ****'\n"
		"0 = NRFLD(= 0 to skip nearfield calc., =1 to calculate nearfield E)\n"
		"0 0 0 0 0 0 = (fract.extens.of calc.vol.in - x, +x, -y, +y, -z, +z)\n"
		"'**** Error Tolerance ****'\n"
		"1e-05 = TOL = MAX ALLOWED(NORM OF | G >= AC | E>-ACA | X>) / (NORM OF AC | E>)\n"
		"'**** Maximum number of iterations ****'\n"
		"300 = MXITER\n"
		"'**** Integration cutoff parameter for PBC calculations ****'\n"
		"0.005 = GAMMA(1e-2 is normal, 3e-3 for greater accuracy)\n"
		"'**** Angular resolution for calculation of <cos>, etc. ****'\n"
		"0.5 = ETASCA(number of angles is proportional to[(3 + x) / ETASCA] ^ 2)\n"
		"'**** Vacuum wavelengths (mm) ****'\n"
		<< wvlen_mm << " " << wvlen_mm << " 1 'LIN' = wavelengths\n"
		"'**** Refractive index of ambient medium'\n"
		"1 = NAMBIENT\n"
		"'**** Effective Radii (mm) **** '\n"
		<< aeff_mm << " " << aeff_mm << " 1 'LIN' = aeff\n"
		"'**** Define Incident Polarizations ****'\n"
		"(0, 0) (1, 0) (0, 0) = Polarization state e01(k along x axis)\n"
		"2 = IORTH(= 1 to do only pol.state e01; = 2 to also do orth.pol.state)\n"
		"'**** Specify which output files to write ****'\n"
		"1 = IWRKSC(= 0 to suppress, =1 to write \".sca\" file for each target orient.\n"
		"'**** Specify Target Rotations ****'\n"
		"0 360 " << nb << "  = BETAMI, BETAMX, NBETA  (beta=rotation around a1)\n"
		"0 180 " << nt << "  = THETMI, THETMX, NTHETA (theta=angle between a1 and k)\n"
		"0 360 " << np << "  = PHIMIN, PHIMAX, NPHI (phi=rotation angle of a1 around k)\n"
		"'**** Specify first IWAV, IRAD, IORI (normally 0 0 0) ****'\n"
		"0 0 0 = first IWAV, first IRAD, first IORI(0 0 0 to begin fresh)\n"
		"'**** Select Elements of S_ij Matrix to Print ****'\n"
		"6 = NSMELTS = number of elements of S_ij to print(not more than 9)\n"
		"11 12 21 22 31 41 = indices ij of elements to print\n"
		"'**** Specify Scattered Directions ****'\n"
		"'LFRAME' = CMDFRM(LFRAME, TFRAME for Lab Frame or Target Frame)\n"
		"2 = NPLANES = number of scattering planes\n"
		"0 0 180 10 = phi, thetan_min, thetan_max, dtheta(in deg) for plane 1\n"
		"90 0 180 10 = phi, thetan_min, thetan_max, dtheta(in deg) for plane 2"
		;
}
