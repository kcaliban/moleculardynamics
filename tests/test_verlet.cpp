#include "verlet.h"
#include <gtest/gtest.h>
#include <cmath>

// Simple molecular dynamics simulation of a single atom with constant force
TEST(VerletIntegrationTest, SingleAtomConstantForceTest) {
    double x = 0, y = 0, z = 0;
    double vx = 0, vy = 0, vz = 0;
    double fx = 1, fy = 1, fz = 1;
    double timestep = 0.1;
    int n_steps = 10;

    for (int i = 0; i < n_steps; i++) {
        verlet_step1(x, y, z, vx, vy, vz, fx, fy, fz, timestep);
        verlet_step2(vx, vy, vz, fx, fy, fz, timestep);
    }

    double analytical_vx = n_steps * timestep * fx;
    double analytical_vy = n_steps * timestep * fy;
    double analytical_vz = n_steps * timestep * fz;
    EXPECT_NEAR(vx, analytical_vx, 1e-6);
    EXPECT_NEAR(vy, analytical_vy, 1e-6);
    EXPECT_NEAR(vz, analytical_vz, 1e-6);

    double analytical_x = 0.5 * pow(n_steps * timestep, 2) * fx;
    double analytical_y = 0.5 * pow(n_steps * timestep, 2) * fy;
    double analytical_z = 0.5 * pow(n_steps * timestep, 2) * fz;
    EXPECT_NEAR(x, analytical_x, 1e-6);
    EXPECT_NEAR(y, analytical_y, 1e-6);
    EXPECT_NEAR(z, analytical_z, 1e-6);
}
