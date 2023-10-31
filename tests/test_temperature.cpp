#include "temperature.h"
#include <gtest/gtest.h>

// Test of calculation of temperature (for regression)
TEST(TemperatureTest, BasicTemperatureWithDimensions) {
    int nb_atoms = 10;
    double m = 1.0;

    // Generate random velocities
    Atoms atoms(nb_atoms);
    atoms.velocities.setRandom();

    // Calculate temperature
    double boltzmann_constant = 8.617333262e-5;
    double kinetic_energy = 0.5 * m * atoms.velocities.square().sum();
    double current_temperature = kinetic_energy * 2 / (3 * atoms.nb_atoms() * boltzmann_constant);

    // Check
    EXPECT_NEAR(temperature(atoms), current_temperature, 1e-10);
}

// Test of calculation of dimensionless temperature
TEST(TemperatureTest, BasicTemperatureDimensionless) {
    int nb_atoms = 10;
    double m = 1.0;

    // Generate random velocities
    Atoms atoms(nb_atoms);
    atoms.velocities.setRandom();

    // Calculate temperature
    double boltzmann_constant = 1;
    double kinetic_energy = 0.5 * m * atoms.velocities.square().sum();
    double current_temperature = kinetic_energy * 2 / (3 * atoms.nb_atoms() * boltzmann_constant);

    // Check
    EXPECT_NEAR(temperature(atoms, true), current_temperature, 1e-10);
}

