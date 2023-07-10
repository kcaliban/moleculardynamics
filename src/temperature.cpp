#include "temperature.h"
#include "kinetic.h"

double temperature(const Atoms &atoms) {
    // double boltzmann_constant = 1.380649e-23;
    double boltzmann_constant = 1;
    double temperature = kinetic_energy(atoms.velocities, atoms.m) * 2 / (3 * atoms.nb_atoms() * boltzmann_constant);
    return temperature;
}