#include "verlet.h"

void verlet_step1(Eigen::Array3Xd &positions, Eigen::Array3Xd &velocities,
                  const Eigen::Array3Xd &forces, double timestep) {
    velocities += 0.5 * forces * timestep;
    positions += velocities * timestep;
}

void verlet_step2(Eigen::Array3Xd &velocities, const Eigen::Array3Xd &forces, double timestep) {
    velocities += 0.5 * forces * timestep;
}