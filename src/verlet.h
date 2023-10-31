#ifndef VERLET_H
#define VERLET_H

#include <Eigen/Dense>

/**
 * Perform the first step of the Velocity-Verlet algorithm
 * 
 * @param positions Positions of the atoms in (Å) or dimensionless
 * @param velocities Velocities of the atoms in (Å/fs) or dimensionless
 * @param forces Forces on the atoms in (eV/Å) or dimensionless
 * @param timestep Timestep in (fs) or dimensionless
 * @param mass Mass of the atoms
*/
void verlet_step1(Eigen::Array3Xd &positions, Eigen::Array3Xd &velocities,
                  const Eigen::Array3Xd &forces, double timestep, double mass = 1.0);

/**
 * Perform the second step of the Velocity-Verlet algorithm
 * 
 * @param velocities Velocities of the atoms in (Å/fs) or dimensionless
 * @param forces Forces on the atoms in (eV/Å) or dimensionless
 * @param timestep Timestep in (fs) or dimensionless
 * @param mass Mass of the atoms
*/
void verlet_step2(Eigen::Array3Xd &velocities, const Eigen::Array3Xd &forces, double timestep, double mass = 1.0);

#endif  // VERLET_H