### Lattice Growth Simulation

## Contents

1. [Main](#main)
2. [Reticolo class](#reticolo-class)
3. [Island class](#islands-class)

## Main
#### Overview
This simulation code models the growth dynamics of a lattice system using the `Reticolo` class. It incorporates processes such as diffusion and deposition to mimic real-world phenomena, enabling the study of various physical systems.

#### Key Features
- **Lattice Structure**: Defines a square lattice with a user-specified size `L`.
- **Parameters**: Allows customization of parameters such as diffusion rate (`D/F`), filling percentage (`theta`), and interaction strength (`nu`).
- **Visualization**: Utilizes ROOT histograms to visualize the evolution of the lattice over time.
- **Dynamic Evolution**: Implements dynamic growth evolution by continuously updating the lattice until a maximum filling percentage is reached.
- **Post-Growth Diffusion**: Conducts additional diffusion iterations after reaching maximum filling to simulate further lattice relaxation.

#### Public Functions

##### Main Function
- **int main(int argc, char **argv)**: The entry point of the program. Initializes the lattice, performs dynamic growth, and visualizes the evolution using ROOT histograms.

#### Utility Functions
- **int approx(float a)**: Approximates a floating-point number to the nearest integer.

#### How to Use
1. **Compilation**: Compile the code using a C++ compiler, ensuring that ROOT libraries are linked.
2. **Execution**: Run the compiled executable, providing the desired `D/F` ratio as a command-line argument.
3. **Visualization**: View the evolution of the lattice dynamics in the ROOT canvas window. The simulation continues until the lattice reaches the maximum filling percentage, followed by additional diffusion iterations.

#### Example
```
./lattice_simulation 1000
```
This command executes the simulation with a `D/F` ratio of 1000.

#### Dependencies
- ROOT libraries for visualization (`TRandom3.h`, `TGraph.h`, `TApplication.h`, etc.).
- `Reticolo.hpp` header file defining the `Reticolo` class for lattice operations.

# Reticolo Class 

#### Overview
The `Reticolo` class represents a lattice structure used for modeling various physical processes such as diffusion and growth. This version of the class introduces several enhancements and modifications for improved functionality and flexibility.

#### Key Updates
- **Single Random Number Generator**: Now uses a single random number generator initialized with the seed `time(0)`, ensuring consistent and reproducible results across different instances.
- **Diffusion Process**: Implements a diffusion process only when a lattice site has 0 neighboring particles, providing more realistic modeling of diffusion phenomena.
- **Periodic Boundary Conditions**: Introduces periodic boundary conditions for lattice operations, enabling simulations of systems with cyclic boundaries.

#### Public Methods

##### Constructors
- **Reticolo(int rows = 10, int cols = 10, double piene = 0)**: Constructs a lattice grid with the specified number of rows and columns, optionally filling it randomly based on the given filling percentage.

##### Printing Methods
- **Print()**: Prints the current lattice configuration to the terminal.
- **Print(string nomefile)**: Prints the current lattice configuration to a file specified by `nomefile`.
- **PrintGr(TH2I &grP)**: Updates a `TH2I` histogram `grP` to visualize the current lattice configuration.

##### Initialization and Modification Methods
- **Init()**: Initializes the lattice grid with zeros.
- **RandomFill(int Nfilled)**: Randomly fills the lattice with a specified number of particles.
- **Deposizione()**: Deposits a particle in an empty site of the lattice.
- **Diffusione(int x, int y)**: Performs a diffusion move from the specified position `(x, y)`.
- **CreaClassi()**: Divides lattice sites into classes based on certain conditions.
- **Crescita(double nu, double DovF, double pDep, bool flag)**: Simulates lattice growth with specified parameters and a flag indicating dynamic evolution.

##### Computation and Analysis Methods
- **ContaParticelle()**: Counts the number of filled sites in the lattice.
- **NBounds(int x, int y)**: Counts the number of neighboring particles at a given site `(x, y)`.
- **ComputeParOrdine()**: Computes the order parameter of the lattice.

##### Getters
- **GetMatrix()**: Retrieves the current lattice matrix.
- **GetNPart()**: Retrieves the total number of particles in the lattice.
- **GetParOrdine()**: Retrieves the order parameter of the lattice.
- **GetClassi()**: Retrieves the classes of lattice sites.
- **GetRandom()**: Generates a random integer within the column dimension of the lattice.

#### Note
This version of the `Reticolo` class is designed to provide enhanced functionality for simulating lattice-based systems with improved accuracy and flexibility.

# Islands class

This class contains an implementation of an algorithm to count the number of islands in a given boolean 2D matrix. An "island" here refers to a group of adjacent cells containing '1's (indicating land), where adjacency includes horizontally, vertically, and diagonally neighboring cells.

## Contents

1. [Introduction](#introduction)
2. [Dependencies](#dependencies)
3. [Usage](#usage)
4. [Algorithm Overview](#algorithm-overview)
5. [Example](#example)

## Introduction

The provided code implements the algorithm to count islands in a given boolean 2D matrix. It utilizes Depth-First Search (DFS) to explore connected cells within the matrix and identify separate islands.

## Usage

To use the provided code:

1. Include the "Islands.hpp" header file in your project.
2. Create an instance of the `Islands` class.
3. Call the `countIslands` function with the boolean 2D matrix as input to get the count of islands.

## Algorithm Overview

The algorithm follows these main steps:

1. Initialize a 2D boolean matrix to keep track of visited cells.
2. Iterate through each cell of the input matrix:
    - If a cell is land (denoted by '1') and has not been visited yet, initiate a DFS traversal to mark all connected cells of that island.
    - Increment the count each time a new island is encountered during the traversal.
3. Return the count of islands found in the matrix.

## Example

```cpp
#include "Islands.hpp"
#include <vector>
#include <iostream>

int main() {
    // Example boolean 2D matrix representing land (1) and water (0)
    std::vector<std::vector<int>> matrix = {
        {1, 0, 0, 0, 1},
        {1, 1, 0, 0, 0},
        {0, 0, 0, 1, 1},
        {0, 0, 0, 0, 0},
        {1, 1, 0, 1, 1}
    };

    // Create an instance of the Islands class
    Islands islandCounter;

    // Get the count of islands
    int islandCount = islandCounter.countIslands(matrix);

    // Output the count of islands
    std::cout << "Number of islands: " << islandCount << std::endl;

    return 0;
}
