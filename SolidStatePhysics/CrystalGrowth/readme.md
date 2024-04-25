

# Reticolo Class

The `Reticolo` class represents a lattice structure and provides methods to manipulate and analyze it. It is designed to simulate lattice growth processes and perform computations related to lattice properties.

## Contents

1. [Introduction](#introduction)
2. [Class Members](#class-members)
3. [Methods](#methods)
4. [Usage](#usage)
5. [Example](#example)

## Introduction

The `Reticolo` class is used to model and analyze lattice structures, commonly encountered in various scientific and engineering applications. It provides functionality to initialize, fill, manipulate, and analyze the lattice, as well as simulate growth processes.

## Class Members

- `m_rows`: Integer representing the number of rows in the lattice.
- `m_cols`: Integer representing the number of columns in the lattice.
- `m_Npart`: Double representing the number of filled sites in the lattice.
- `m_matr`: 2D vector representing the lattice structure.
- `m_classi`: Array of vectors representing different classes of lattice points based on their neighboring sites.
- `m_rand`: Random number generator object.
- `m_p`: Double representing the order parameter of the lattice.

## Methods

- `Init()`: Initializes the lattice to all empty sites.
- `RandomFill(int)`: Fills the lattice with a given percentage of filled sites randomly.
- `Diffusione(int, int)`: Performs a diffusion move from a specified position in the lattice.
- `Deposizione()`: Fills an empty site of the lattice.
- `ComputeParOrdine()`: Computes the order parameter of the lattice.
- `CreaClassi()`: Divides lattice points into classes based on their neighboring sites.
- `Crescita(double, double, double, double, double)`: Simulates lattice growth using a Monte Carlo algorithm.
- `NBounds(int, int)`: Counts the number of bonds at a specified position in the lattice.
- `ContaParticelle()`: Counts the filled sites of the lattice.
- `GetNPart()`: Returns the number of filled sites in the lattice.
- `GetParOrdine()`: Returns the order parameter of the lattice.
- `GetClassi()`: Returns the array of vectors representing different classes of lattice points.

## Usage

1. Include the "Reticolo.hpp" header file in your project.
2. Create an instance of the `Reticolo` class with the desired lattice dimensions.
3. Use the provided methods to initialize, fill, manipulate, and analyze the lattice.
