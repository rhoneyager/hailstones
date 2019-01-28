# hailstones

Codes for modeling hailstone structural and electromagnetic scattering properties from stereolithographic measurements.

![Picture of a modeled hailstone](./Hailstone%20example.PNG?raw=true)

## Building

Requirements:
- A C++14 compiler
- CMake
- Boost
- VTK 6

1. Run CMake and generate the build scripts.
2. Build the project. Usually, this involves typing "make"
3. Run the programs. 

## The Programs

### Prog1 - Read STL surface meshes and generate volumetric meshes
### Prog2 - Take volumetric meshes and scale them to that they match the measured hailstone masses / densities.
### Prog3a - Convert the volumetric meshes into [DDSCAT](http://ddscat.wikidot.com/) shape files
### Prog3b - Convert the volumetric meshes into [DDSCAT](http://ddscat.wikidot.com/) runs at desired spectral frequency and temperature.
