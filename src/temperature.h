#ifndef TEMPERATURE_H
#define TEMPERATURE_H

#include "atoms.h"

/**
* Calculate the temperature for given atoms
* 
* @param atoms Atoms to calculate the temperature for
* @param dimensionless Whether not to calculate temperature using Boltzmann constant
* @return Temperature in Kelvin or dimensionless
*/
double temperature(const Atoms &atoms, bool dimensionless = false);

/**
 * Calculate the temperature for given kinetic energy and number of atoms
 * 
 * @param kinetic_energy Kinetic energy to calculate the temperature for
 * @param nb_atoms Number of atoms
 * @param dimensionless Whether not to calculate temperature using Boltzmann constant
 * @return Temperature in Kelvin or dimensionless
 */
double temperature(const double kinetic_energy, const int nb_atoms, bool dimensionless = false);

#endif // TEMPERATURE_H