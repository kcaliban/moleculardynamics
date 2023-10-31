#include "lj_direct_summation.h"
#include <iostream>

constexpr double lj(double r, double epsilon, double sigma) {
    // Lennard-Jones potential is:
    // V(r) = 4 * epsilon * ((sigma / r)^12 - (sigma / r)^6)
    // where r is the distance between two atoms
    return 4 * epsilon * (pow(sigma / r, 12) - pow(sigma / r, 6));
}

constexpr double lj_derivative(double r, double epsilon, double sigma) {
    // Derivative of Lennard-Jones potential is:
    // dV/dr = 4 * epsilon * (12 * sigma^12 / r^13 - 6 * sigma^6 / r^7)
    // where r is the distance between two atoms
    return 4 * epsilon * ((12 * pow(sigma, 12) / pow(r, 13)) - (6 * pow(sigma, 6) / pow(r, 7)));
}

double lj_direct_summation(Atoms &atoms, double epsilon = 1.0, double sigma = 1.0) {
    // Initialize potential energy and forces to zero
    double potential_energy{0.0};
    atoms.forces.setZero();

    // Loop over all pairs of atoms
    for (size_t i = 0; i < atoms.nb_atoms(); i++) {
        for (size_t j = i + 1; j < atoms.nb_atoms(); j++) {
            // Calculate distance between atoms
            Eigen::Array3Xd diff{atoms.positions.col(i) - atoms.positions.col(j)};
            double r{sqrt(diff.square().sum())};

            // Calculate potential energy
            potential_energy += lj(r, epsilon, sigma);

            // Calculate forces
            double lj_derived{lj_derivative(r, epsilon, sigma)};
            atoms.forces.col(i) += lj_derived * (diff / r);
            atoms.forces.col(j) -= lj_derived * (diff / r);
        }
    }

    return potential_energy;
}
