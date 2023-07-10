#ifndef BERENDSEN_H
#define BERENDSEN_H

#include "atoms.h"

void berendsen_thermostat(Atoms &atoms, double target_temperature, double timestep, double relaxation_time);

#endif // BERENDSEN_H