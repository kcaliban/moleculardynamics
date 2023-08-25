#include "kinetic.h"
#include <iostream>

/**
* Calculate the kinetic energy given velocities and mass
* 
* @param velocities Velocities of the atoms in (Å/fs)
* @param m Mass of the atoms
* @return Kinetic energy in eV
*/
double kinetic_energy(const Velocities_t &velocities, double m) {
    return 0.5 * m * velocities.square().sum();
}