#include "kinetic.h"

double kinetic_energy(const Velocities_t &velocities, double m) {
    // Kinetic energy is:
    // E_kin = 0.5 * m * v^2
    // where m is the mass of the atom and v is the velocity
    return 0.5 * m * velocities.square().sum();
}