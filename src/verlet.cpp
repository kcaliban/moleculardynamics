#include "verlet.h"

void verlet_step1(Eigen::Array3Xd &positions, Eigen::Array3Xd &velocities,
                  const Eigen::Array3Xd &forces, double timestep) {
    velocities.row(0) += 0.5 * forces.row(0) * timestep;
    velocities.row(1) += 0.5 * forces.row(1) * timestep;
    velocities.row(2) += 0.5 * forces.row(2) * timestep;
    positions.row(0) += velocities.row(0) * timestep;
    positions.row(1) += velocities.row(1) * timestep;
    positions.row(2) += velocities.row(2) * timestep;

    // velocities += 0.5 * forces * timestep;
    // positions += velocities * timestep;
}

void verlet_step2(Eigen::Array3Xd &velocities, const Eigen::Array3Xd &forces, double timestep) {
    velocities.row(0) += 0.5 * forces.row(0) * timestep;
    velocities.row(1) += 0.5 * forces.row(1) * timestep;
    velocities.row(2) += 0.5 * forces.row(2) * timestep;

    // velocities += 0.5 * forces * timestep;
}

/*
void verlet_step1(double &x, double &y, double &z, double &vx, double &vy, double &vz,
                  double fx, double fy, double fz, double timestep) {
    vx += 0.5 * fx * timestep;
    vy += 0.5 * fy * timestep;
    vz += 0.5 * fz * timestep;
    x += vx * timestep;
    y += vy * timestep;
    z += vz * timestep;
}

void verlet_step2(double &vx, double &vy, double &vz, double fx, double fy, double fz,
                  double timestep) {
    vx += 0.5 * fx * timestep;
    vy += 0.5 * fy * timestep;
    vz += 0.5 * fz * timestep;
}
*/