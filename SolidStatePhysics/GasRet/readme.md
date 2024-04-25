### README

#### Introduction
This C++ program simulates the evolution of a lattice system using the Metropolis algorithm. The lattice represents a two-dimensional grid with particles occupying some of the grid sites. The program evolves the system over a certain number of Monte Carlo steps (MCS) and tracks the order parameter of the lattice.

#### Requirements
- C++ compiler
- ROOT framework (for plotting)

#### Files
- **Main.cpp**: Contains the main code for the simulation.
- **Matrice.hpp**: Header file for the Matrice class, which defines the lattice and its operations.
- **TRandom3.h**: Header file for random number generation.
- **TGraph.h**, **TApplication.h**, **TAxis.h**, **TCanvas.h**, **TStyle.h**, **TSystem.h**: ROOT framework headers for plotting.
- **TH2I.h**: Header file for 2D histograms in ROOT.

#### Compilation
To compile the program, you need to link against the ROOT library. Use the following command:

```bash
g++ Main.cpp -o Main `root-config --cflags --libs`
```

#### Usage
The program prompts the user to input parameters interactively:
- **L**: Size of the lattice grid (must be even).
- **theta**: Fill percentage of the lattice with particles.
- **N_mcs**: Number of Monte Carlo steps for the simulation.
- **BetaJ**: Parameter for energy calculation.
- **dyn**: Option to enable dynamic evolution plotting (1 for yes, 0 for no).

#### Output
- The program generates output files named "output0.dat", "output1.dat", and "output2.dat", containing the state of the lattice at different stages of the simulation.
- Plots of the lattice at the initial, intermediate, and final stages are displayed using ROOT.
- A graph showing the evolution of the order parameter is also displayed.

#### Notes
- The simulation time increases significantly with larger lattice sizes.
- Adjust the number of MCS for convergence based on the lattice size and fill percentage.

#### Matrice Class
- **GetMatrix()**: Returns the matrix representing the lattice.
- **ContaParticelle()**: Counts the number of particles in the lattice.
- **Print()**: Prints the lattice configuration to the console.
- **Print(string nomefile)**: Writes the lattice configuration to a file with the specified name.
- **PrintGr(TH2I &grP)**: Updates a 2D histogram (TH2I) with the lattice configuration.
- **Init()**: Initializes the lattice with all sites empty.
- **RandomFill(int Nfilled)**: Randomly fills the lattice with a specified number of particles.
- **NBounds(int x, int y)**: Computes the number of neighboring particles for a given lattice site.
- **RandomMove(int &successi)**: Performs a random move of particles according to the Metropolis algorithm.
- **ComputeParOrdine()**: Computes the order parameter of the lattice.
- **GetParOrdine()**: Returns the computed order parameter.
