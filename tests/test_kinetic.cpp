#include "kinetic.h"
#include <gtest/gtest.h>

// Test of calculation of kinetic energy (for regression)
TEST(KineticEnergyTest, Simple) {
    // Generate random velocities
    Velocities_t velocities(3, 1);
    velocities.setRandom();

    // Calculate kinetic energy
    double m = 1.0;
    double expected = 0.5 * m * velocities.square().sum();
    double actual = kinetic_energy(velocities, m);

    // Check
    EXPECT_NEAR(expected, actual, 1e-6);
}