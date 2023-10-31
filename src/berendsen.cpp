#include "berendsen.h"
#include "temperature.h"
#include <iostream>

void berendsen_thermostat(Atoms &atoms, double target_temperature, double timestep, double relaxation_time, bool dimensionless) {
    double current_temperature{temperature(atoms, dimensionless)};

    if (current_temperature) {
        atoms.velocities *= sqrt(1 + timestep / relaxation_time * (target_temperature / current_temperature - 1));
    }
}

void berendsen_thermostat(Atoms &atoms, double kinetic_energy, double target_temperature, 
                          double timestep, double relaxation_time, int total_number_atoms, bool dimensionless) {
    double current_temperature{temperature(kinetic_energy, total_number_atoms, dimensionless)};

    if (current_temperature) {
        atoms.velocities *= sqrt(1 + timestep / relaxation_time * (target_temperature / current_temperature - 1));
    }
}