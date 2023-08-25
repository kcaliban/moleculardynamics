#ifndef VERLET_H
#define VERLET_H

#include <Eigen/Dense>

void verlet_step1(Eigen::Array3Xd &positions, Eigen::Array3Xd &velocities,
                  const Eigen::Array3Xd &forces, double timestep, double mass = 1.0);
void verlet_step2(Eigen::Array3Xd &velocities, const Eigen::Array3Xd &forces, double timestep, double mass = 1.0);

#endif  // VERLET_H