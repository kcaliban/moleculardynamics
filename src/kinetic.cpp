#include "kinetic.h"
#include <iostream>

double kinetic_energy(const Velocities_t &velocities, double m) {
    return 0.5 * m * velocities.square().sum();
}

double kinetic_energy_subset(const Velocities_t &velocities, double m, int n) {
    return 0.5 * m * velocities.leftCols(n).square().sum();
}