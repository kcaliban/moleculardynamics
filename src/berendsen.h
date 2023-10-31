#ifndef BERENDSEN_H
#define BERENDSEN_H

#include "atoms.h"

/**
 * Berendsen thermostat
 * 
 * @param atoms Atoms to apply thermostat to
 * @param target_temperature Target temperature
 * @param timestep Timestep
 * @param relaxation_time Relaxation time
 * @param dimensionless Whether not to calculate temperature using Boltzmann constant
*/
void berendsen_thermostat(Atoms &atoms, double target_temperature, double timestep, double relaxation_time, bool dimensionless = false);

/**
 * Berendsen thermostat
 * 
 * @param atoms Atoms to apply thermostat to
 * @param kinetic_energy Kinetic energy
 * @param target_temperature Target temperature
 * @param timestep Timestep
 * @param relaxation_time Relaxation time
 * @param total_number_atoms Total number of atoms
 * @param dimensionless Whether not to calculate temperature using Boltzmann constant
*/
void berendsen_thermostat(Atoms &atoms, double kinetic_energy, double target_temperature, double timestep, 
                          double relaxation_time, int total_number_atoms, bool dimensionless = false);

#endif // BERENDSEN_H