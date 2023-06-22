#include "kinetic.h"
#include <gtest/gtest.h>

// Test of calculation of kinetic energy (pretty stupid test, for regression)
TEST(KineticEnergyTest, Simple) {
    Velocities_t velocities(3, 1);
    velocities.setRandom();

    double m = 1.0;
    double expected = 0.5 * m * velocities.square().sum();
    double actual = kinetic_energy(velocities, m);
    EXPECT_NEAR(expected, actual, 1e-6);
}