#ifndef CUBICLATTICE_H
#define CUBICLATTICE_H

#include "atoms.h"

/**
 * Create a cubic lattice
 * 
 * @param lattice_constant Lattice constant
 * @param atoms_per_side Atoms per side
 * @return Atoms in a cubic lattice with the given lattice constant and number of atoms per side
*/
Atoms create_cubic_lattice(double lattice_constant, int atoms_per_side);

#endif // CUBICLATTICE_H