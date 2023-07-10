#include "temperature.h"
#include <gtest/gtest.h>

TEST(TemperatureTest, BasicTemperature) {
    constexpr int nb_atoms = 10;
    constexpr double m = 1.0;

    Atoms atoms(nb_atoms);
    atoms.velocities.setRandom();

    // double boltzmann_constant = 1.380649e-23;
    double boltzmann_constant = 1;
    double kinetic_energy = 0.5 * m * atoms.velocities.square().sum();
    double current_temperature = kinetic_energy * 2 / (3 * atoms.nb_atoms() * boltzmann_constant);

    EXPECT_NEAR(temperature(atoms), current_temperature, 1e-10);
}
