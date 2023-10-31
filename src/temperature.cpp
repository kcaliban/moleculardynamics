#include "temperature.h"
#include "kinetic.h"

double temperature(const Atoms &atoms, bool dimensionless) {
    double boltzmann_constant{dimensionless ? 1 : 8.617'333'262e-5}; // dimensionless or eV/K
    return kinetic_energy(atoms.velocities, atoms.m) * 2 / (3 * atoms.nb_atoms() * boltzmann_constant);
}

double temperature(const double kinetic_energy, const int nb_atoms, bool dimensionless) {
    double boltzmann_constant{dimensionless ? 1 : 8.617'333'262e-5}; // dimensionless or eV/K
    return kinetic_energy * 2 / (3 * nb_atoms * boltzmann_constant);
}