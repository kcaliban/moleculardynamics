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

/**
 * Calculate the temperature for kinetic energy
 * 
 * @param kinetic_energy Kinetic energy to calculate the temperature for
 * @param nb_atoms Number of atoms
 * @return Temperature in Kelvin
 */
double temperature(const double kinetic_energy, const int nb_atoms) {
    double boltzmann_constant = 8.617333262e-5; // eV/K
    double temperature = kinetic_energy * 2 / (3 * nb_atoms * boltzmann_constant);
    return temperature;
}