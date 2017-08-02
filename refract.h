#pragma once
#include <complex>

#include <vtkSmartPointer.h>
#include <vtkImageData.h>

void mIceMatzler(double fGHz, double tK, std::complex<double> &m);
void writeDiel(const char* file, std::complex<double>& r);
void writePar(const char* file, double wvlen_mm, double aeff_mm, unsigned short nb, unsigned short nt, unsigned short np);
void writeShp(const char* file, const char* desc, vtkSmartPointer<vtkImageData> imgd, double c[3], size_t numPoints);
