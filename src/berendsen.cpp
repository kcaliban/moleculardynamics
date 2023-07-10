#include "berendsen.h"
#include "temperature.h"

void berendsen_thermostat(Atoms &atoms, double target_temperature, double timestep, double relaxation_time) {
    double current_temperature = temperature(atoms);
    double lambda = sqrt(1 + timestep / relaxation_time * (target_temperature / current_temperature - 1));
    atoms.velocities *= lambda;
}
