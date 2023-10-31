#include "cubic_lattice.h"

Atoms create_cubic_lattice(double lattice_constant, int atoms_per_side) {
    // Create atoms
    Atoms atoms(pow(atoms_per_side, 3));

    // Set atom positions
    for (int i = 0; i < atoms_per_side; i++) {
        for (int j = 0; j < atoms_per_side; j++) {
            for (int k = 0; k < atoms_per_side; k++) {
                int index = i * atoms_per_side * atoms_per_side + j * atoms_per_side + k;
                atoms.positions.col(index)(0) = i * lattice_constant;
                atoms.positions.col(index)(1) = j * lattice_constant;
                atoms.positions.col(index)(2) = k * lattice_constant;
            }
        }
    }

    return atoms;
}