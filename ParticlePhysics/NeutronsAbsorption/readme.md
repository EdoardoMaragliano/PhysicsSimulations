
# Neutron Interaction Simulation

## Overview
This program simulates the interaction of neutrons with a detector material. It calculates the energy deposition and tracks the trajectories of neutrons within the detector. The simulation includes variations in neutron flux and detection efficiency.

## Compilation
To compile the program, execute the following command in the terminal:
```
g++ SimNeuVar.cpp -o SimNeuVar.x `root-config --cflags --glibs`
```
This command compiles the code using the ROOT libraries.

## Usage
To run the simulation, execute the compiled executable:
```
./SimNeuVar.x
```
The program simulates the interaction of neutrons with the detector and displays the results.

## Dependencies
The program requires the ROOT framework for data analysis and visualization. Ensure that ROOT is properly installed on your system before compiling and running the code.

## Input Parameters
- `A`: Atomic mass number.
- `n`: Number density of detector material nuclei.
- `dz`: Thickness of the detector.
- `detr`: Radius of the detector.
- `E`: Initial energy of neutrons.
- `nneu`: Number of neutrons to simulate.
- `ngr`: Number of events to display graphically.

## Output
The program generates histograms displaying the energy deposition in the detector and calculates the fraction of neutrons detected. It also calculates the figure of merit (FOM) to evaluate the performance of the detection system.
