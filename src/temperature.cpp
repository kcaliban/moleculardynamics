#include "temperature.h"
#include "kinetic.h"

/**
* Calculate the temperature for given atoms
* 
* @param atoms Atoms to calculate the temperature for
* @return Temperature in Kelvin
*/
double temperature(const Atoms &atoms) {
    double boltzmann_constant = 8.617333262e-5; // eV/K
    double temperature = kinetic_energy(atoms.velocities, atoms.m) * 2 / (3 * atoms.nb_atoms() * boltzmann_constant);
    return temperature;
}