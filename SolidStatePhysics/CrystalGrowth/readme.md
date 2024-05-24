

# Lattice Class

The `Lattice` class represents a lattice structure and provides methods to manipulate and analyze it. It is designed to simulate lattice growth processes and perform computations related to lattice properties.

## Contents

1. [Introduction](#introduction)
2. [Class Members](#class-members)
3. [Methods](#methods)
4. [Usage](#usage)
5. [Example](#example)

## Introduction

The `Lattice` class is used to model and analyze lattice structures, commonly encountered in various scientific and engineering applications. It provides functionality to initialize, fill, manipulate, and analyze the lattice, as well as simulate growth processes.

## Class Members

- `m_rows`: Integer representing the number of rows in the lattice.
- `m_cols`: Integer representing the number of columns in the lattice.
- `m_Npart`: Double representing the number of filled sites in the lattice.
- `m_matr`: 2D vector representing the lattice structure.
- `m_classes`: Array of vectors representing different classes of lattice points based on their neighboring sites.
- `m_rand`: Random number generator object.
- `m_p`: Double representing the order parameter of the lattice.

## Methods

- `Init()`: Initializes the lattice to all empty sites.
- `RandomFill(int)`: Fills the lattice with a given percentage of filled sites randomly.
- `Diffusion(int, int)`: Performs a diffusion move from a specified position in the lattice.
- `Deposition()`: Fills an empty site of the lattice.
- `ComputeOrderPar()`: Computes the order parameter of the lattice.
- `CreateClasses()`: Divides lattice points into classes based on their neighboring sites.
- `Growth(double, double, double, double, double)`: Simulates lattice growth using a Monte Carlo algorithm.
- `NBounds(int, int)`: Counts the number of bonds at a specified position in the lattice.
- `CountNparticles()`: Counts the filled sites of the lattice.
- `GetNPart()`: Returns the number of filled sites in the lattice.
- `GetOrderPar()`: Returns the order parameter of the lattice.
- `GetClasses()`: Returns the array of vectors representing different classes of lattice points.

## Usage

1. Include the "Lattice.hpp" header file in your project.
2. Create an instance of the `Lattice` class with the desired lattice dimensions.
3. Use the provided methods to initialize, fill, manipulate, and analyze the lattice.
