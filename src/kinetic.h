#ifndef KINETIC_H
#define KINETIC_H

#include "types.h"

/**
* Calculate the kinetic energy given velocities and mass
* 
* @param velocities Velocities of the atoms in (Å/fs) or dimensionless
* @param m Mass of the atoms
* @return Kinetic energy in eV or dimensionless
*/
double kinetic_energy(const Velocities_t &velocities, double m = 1.0);

/**
 * Calculate kinetic energy for first n atoms
 * 
 * @param velocities Velocities of the atoms in (Å/fs) or dimensionless
 * @param m Mass of the atoms
 * @param n Number of atoms to calculate kinetic energy for
 * @return Kinetic energy in eV or dimensionless
 */
double kinetic_energy_subset(const Velocities_t &velocities, double m, int n);

#endif // KINETIC_H