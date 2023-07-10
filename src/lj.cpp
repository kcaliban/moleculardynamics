#include "lj.h"
#include <iostream>
#include "neighbors.h"

double lj_pot(double r, double epsilon, double sigma) {
    // Lennard-Jones potential is:
    // V(r) = 4 * epsilon * ((sigma / r)^12 - (sigma / r)^6)
    // where r is the distance between two atoms
    return 4 * epsilon * (pow(sigma / r, 12) - pow(sigma / r, 6));
}

double lj_derivative(double r, double epsilon, double sigma) {
    // Derivative of Lennard-Jones potential is:
    // dV/dr = 4 * epsilon * (12 * sigma^12 / r^13 - 6 * sigma^6 / r^7)
    // where r is the distance between two atoms
    return 4 * epsilon * ((12 * pow(sigma, 12) / pow(r, 13)) - (6 * pow(sigma, 6) / pow(r, 7)));
}

double lj(Atoms &atoms, double epsilon = 1.0, double sigma = 1.0, double cutoff = 1.0) {
    double potential_energy = 0.0;
    atoms.forces.setZero();

    NeighborList neighbor_list;
    neighbor_list.update(atoms, cutoff);

    for (auto[i, j]: neighbor_list) {
        if (i < j) {
            Eigen::Array3Xd diff = atoms.positions.col(i) - atoms.positions.col(j);
            double r = sqrt(diff.square().sum());
            potential_energy += lj_pot(r, epsilon, sigma);

            double lj_derived = lj_derivative(r, epsilon, sigma);
            atoms.forces.col(i) += lj_derived * (diff / r);
            atoms.forces.col(j) -= lj_derived * (diff / r);
        }
    }

    return potential_energy;
}
