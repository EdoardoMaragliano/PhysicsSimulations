
# Na22 decay

## Overview
This program simulates the interaction of Na22 decay products with a detector material. It calculates the energy deposition and tracks the trajectories of neutrons within the detector. The simulation includes different scenarios such as biased beam, isotropic radiation, and variable bias conditions.

## Compilation
To compile the program, execute the following command in the terminal:
```
g++ bias_var.cpp -o bias_var.x `root-config --cflags --glibs`
```
This command compiles the code using the ROOT libraries.

## Usage
To run the simulation, execute the compiled executable with one of the following options:
- `./bias_var.x lin_beam`: Simulates a linear beam of photons.
- `./bias_var.x isotropic`: Simulates isotropic radiation.
- `./bias_var.x bias_var`: Simulates variable bias conditions.

## Dependencies
The program requires the ROOT framework for data analysis and visualization. Ensure that ROOT is properly installed on your system before compiling and running the code.

## Input Parameters
- `nev`: Number of events to simulate.
- `r`: Radius of the detector.
- `dz`: Thickness of the detector.
- `dist`: Distance between the source and the detector.
- `E`: Array containing the energies of the photons.

## Output
The program generates histograms displaying the energy deposition in the detector for each simulation scenario. Additionally, it calculates and prints the figure of merit (FOM) for each simulation, which measures the performance of the detector.

