#include "verlet.h"
#include "types.h"
#include <gtest/gtest.h>
#include <cmath>

// Simple molecular dynamics simulation of a single atom with constant force
TEST(VerletIntegrationTest, SingleAtomConstantForceTest) {
    Positions_t positions(3, 1);
    Velocities_t velocities(3, 1);
    Forces_t forces(3, 1);

    // Set initial conditions
    positions(0, 0) = 0; positions(1, 0) = 0; positions(2, 0) = 0;
    velocities(0, 0) = 0; velocities(1, 0) = 0; velocities(2, 0) = 0;
    forces(0, 0) = 1; forces(1, 0) = 1; forces(2, 0) = 1;

    // Set simulation params
    double timestep = 0.1;
    int n_steps = 10;

    // Run simulation for n_steps
    for (int i = 0; i < n_steps; i++) {
        verlet_step1(positions, velocities, forces, timestep);
        verlet_step2(velocities, forces, timestep);
    }

    // Compare results with analytical solution
    Velocities_t analytical_velocities = n_steps * timestep * forces;
    EXPECT_NEAR(velocities(0, 0), analytical_velocities(0, 0), 1e-6);
    EXPECT_NEAR(velocities(1, 0), analytical_velocities(1, 0), 1e-6);
    EXPECT_NEAR(velocities(2, 0), analytical_velocities(2, 0), 1e-6);

    Positions_t analytical_positions = 0.5 * pow(n_steps * timestep, 2) * forces;
    EXPECT_NEAR(positions(0, 0), analytical_positions(0, 0), 1e-6);
    EXPECT_NEAR(positions(1, 0), analytical_positions(1, 0), 1e-6);
    EXPECT_NEAR(positions(2, 0), analytical_positions(2, 0), 1e-6);
}

// Simple molecular dynamics simulation of multiple atoms with constant force
TEST(VerletIntegrationTest, MultipleAtomConstantForceTest) {
    Positions_t positions(3, 3);
    Velocities_t velocities(3, 3);
    Forces_t forces(3, 3);

    // Set initial conditions
    positions(0, 0) = 0; positions(1, 0) = 0; positions(2, 0) = 0;
    velocities(0, 0) = 1; velocities(1, 0) = 2; velocities(2, 0) = -2;
    forces(0, 0) = 1; forces(1, 0) = 1; forces(2, 0) = 2;

    positions(0, 1) = 1; positions(1, 1) = 3; positions(2, 1) = -2;
    velocities(0, 1) = 1; velocities(1, 1) = -1; velocities(2, 1) = 2;
    forces(0, 1) = -2; forces(1, 1) = 2; forces(2, 1) = -5;

    positions(0, 2) = 2; positions(1, 2) = 3; positions(2, 2) = -5;
    velocities(0, 2) = -1; velocities(1, 2) = 2; velocities(2, 2) = -3;
    forces(0, 2) = 4; forces(1, 2) = 1; forces(2, 2) = 1;

    // Set simulation params
    double timestep = 0.1;
    int n_steps = 10;

    // Calculate analytical solution
    Velocities_t analytical_velocities(3, 3);
    analytical_velocities = velocities + n_steps * timestep * forces;

    Positions_t analytical_positions(3, 3);
    analytical_positions = positions + velocities + 0.5 * pow(n_steps * timestep, 2) * forces;

    // Run simulation for n_steps
    for (int i = 0; i < n_steps; i++) {
        verlet_step1(positions, velocities, forces, timestep);
        verlet_step2(velocities, forces, timestep);
    }

    // Compare results with analytical solution
    for (int i = 0; i < 3; i++) {
        EXPECT_NEAR(velocities(0, i), analytical_velocities(0, i), 1e-6);
        EXPECT_NEAR(velocities(1, i), analytical_velocities(1, i), 1e-6);
        EXPECT_NEAR(velocities(2, i), analytical_velocities(2, i), 1e-6);

        EXPECT_NEAR(positions(0, i), analytical_positions(0, i), 1e-6);
        EXPECT_NEAR(positions(1, i), analytical_positions(1, i), 1e-6);
        EXPECT_NEAR(positions(2, i), analytical_positions(2, i), 1e-6);
    }
}