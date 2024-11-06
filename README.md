
# Exact Diagonalization of the Fermionic Hubbard Model

This package implements Exact Diagonalization (ED) for the Fermionic Hubbard Model. It consists of several Fortran90 source files designed to generate the basis set, construct the Hamiltonian matrix, perform Lanczos diagonalization, and compute physical observables like the ground state energy, spin correlations, and dynamical structure factor. Below is an overview of the code structure and usage.

## File Overview

- **gendet.f90**: This file generates the basis set for the Hubbard model, conserving both spin and particle number.
- **matppp.f90**: Constructs the Hamiltonian matrix in sparse (CSR) format. Originally designed for the Pariser-Parr-Pople (PPP) model, it has been adapted for the Hubbard model.
- **makeinit.f90**: Creates the initial state vector for the Lanczos diagonalization method.
- **diagonly.f90**: Performs Lanczos diagonalization to obtain the ground state energy and expectation values of physical observables such as spin (`<Sz>`) per site and spin correlations.
- **lanc_dsf.f90**: Performs Lanczos diagonalization and extends the Lanczos method to calculate the dynamical structure factor using continued fractions.
- **paramod.f90**: Defines the variables and parameters used in the simulation.
- **interface.f90**: Contains Fortran interfaces for necessary variables and functions.
- **compile**: Shell script to compile all the Fortran90 files.
- **initial.inp**: input file for the parameters


## Requirements

- Fortran90 compiler (e.g., `gfortran`)
- A Unix-like operating system (Linux/MacOS) for the `compile` script (Windows users may need a Unix-like environment such as WSL or Cygwin).

## Installation

1. Clone or download the repository.
2. Navigate to the directory containing the Fortran90 files.
3. Run the `compile` script to compile the source code:
   ```bash
   ./compile

## Usage

1. Input file: 
2. Generating the Basis Set
Run gendet.f90 to generate the basis set for the Fermionic Hubbard Model. This step ensures that the particle number and spin are conserved in your simulations.

3. Constructing the Hamiltonian Matrix
Use matppp.f90 to build the Hamiltonian matrix in sparse format. This matrix will be used in the diagonalization steps to find the ground state energy and compute observables.

4. Performing Lanczos Diagonalization
The Lanczos method is implemented in diagonly.f90. Running this code will compute the ground state energy and observables like spin expectation values per site and spin correlations.

bash
./xlanc
4. Computing the Dynamical Structure Factor
The lanc_dsf.f90 file performs continued fractions in the Lanczos method to compute the dynamical structure factor. This can be used to explore the frequency-dependent properties of the system.

bash
Copy code
./a.out
5. Setting Parameters
The paramod.f90 file allows you to set the physical parameters (e.g., interaction strength, lattice size) for the model. Modify this file before running the simulations.

Parameters
Lattice size: Set in paramod.f90
Interaction strength (U): Set in paramod.f90
Particle number: Defined in gendet.f90
Example Output
Ground state energy: Printed by diagonly.f90.
Spin correlation: Computed and printed by diagonly.f90.
Dynamical structure factor: Output by lanc_dsf.f90.


