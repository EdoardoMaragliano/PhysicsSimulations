# dE/dx + TOF Simulation

## Introduction

This C++ program simulates the specific energy loss (dE/dx) and Time of Flight (TOF) measurements for charged particles passing through a material. It computes the dE/dx and TOF at a fixed particle momentum, studying truncated means.

## Features

- **Bethe-Bloch Formula**: Calculates the specific energy loss (dE/dx) for charged particles passing through a material.
- **Time of Flight (TOF) Simulation**: Simulates the TOF measurements for the particles.
- **Truncated Means Analysis**: Computes truncated means of the dE/dx measurements.
- **Visualization**: Plots the dE/dx distributions, RMS, and difference between reconstructed and generated energy.
- **Customizable Settings**: Allows customization of material properties (atomic number, mass density) and particle type.

## Dependencies

- **ROOT Framework**: Required for graph plotting and manipulation.
- **C++ Compiler**: Compile using a C++11 compatible compiler.

## Compilation

Compile the program with the following command:

```bash
g++ --std=c++11 -o dedx_tof_simulation.x dedx_tof_simulation.cpp `root-config --cflags --glibs`
```

# Particle Identification Simulation

## Introduction

This C++ program simulates particle identification (PID) in high-energy physics experiments. It generates events corresponding to particle decays, simulates detector effects, and attempts to identify particles based on energy loss (dE/dx) and time-of-flight (TOF) measurements.

## Features

- **Particle Generation**: Generates events for particle decays (e.g., D*, D0, kaons, pions).
- **Detector Simulation**: Simulates detector resolution effects (e.g., momentum, angle, energy loss).
- **Particle Identification**: Identifies particles using measured dE/dx and TOF values.
- **Customizable Settings**: Customize event count, D0 decay probability, and enable/disable TOF and dE/dx PID.
- **Output Analysis**: Produces histograms showing invariant mass distribution with and without PID, and provides statistics on particle identification performance.

## Dependencies

- **ROOT Framework**: Required for histogram plotting and manipulation.
- **C++ Compiler**: Compile using a C++11 compatible compiler.

## Compilation

Compile the program using the provided Makefile or manually with the following command:

```bash
g++ --std=c++11 -o particle_id.x particle_id.cpp `root-config --cflags --glibs`
```

## Usage

Execute the program with default or custom settings:

```bash
./particle_id.x [nEv] [pD0] [TOF] [dEdx]
```

- `nEv`: Number of events to simulate.
- `pD0`: Probability of D0 decay.
- `TOF`: Enable/disable TOF PID (on/off).
- `dEdx`: Enable/disable dE/dx PID (on/off).

## Output

Generates histograms displaying invariant mass distribution with and without PID. Prints statistics on particle identification performance.

