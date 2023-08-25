#ifndef TEMPERATURE_H
#define TEMPERATURE_H

#include "atoms.h"

double temperature(const Atoms &atoms);
double temperature(const double kinetic_energy, const int nb_atoms);

#endif // TEMPERATURE_H