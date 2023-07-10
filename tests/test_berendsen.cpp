#include "berendsen.h"
#include "verlet.h"
#include "temperature.h"
#include "lj_direct_summation.h"
#include <gtest/gtest.h>

TEST(BerendsenTest, ReachesTemperatureSingleStep) {
    constexpr int nb_atoms = 10;
    constexpr double target_temperature = 300.0;

    Atoms atoms(nb_atoms);
    atoms.m = 1.0;
    atoms.positions.setRandom();
    atoms.velocities.setRandom();

    double timestep = 1;
    int n_steps = 1;

    for (int i = 0; i < n_steps; i++) {
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep);
        verlet_step2(atoms.velocities, atoms.forces, timestep);
        berendsen_thermostat(atoms, target_temperature, timestep, 1);
    }

    EXPECT_NEAR(temperature(atoms), target_temperature, 1e-1);
}

TEST(BerendsenTest, ReachesTemperatureMultipleSteps) {
    constexpr int nb_atoms = 10;
    constexpr double target_temperature = 300.0;

    Atoms atoms(nb_atoms);
    atoms.m = 1.0;
    atoms.positions.setRandom();
    atoms.velocities.setRandom();

    double timestep = 1;
    int n_steps = 20;

    for (int i = 0; i < n_steps; i++) {
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep);
        verlet_step2(atoms.velocities, atoms.forces, timestep);
        berendsen_thermostat(atoms, target_temperature, timestep, 2);
    }

    EXPECT_NEAR(temperature(atoms), target_temperature, 1e-1);
}
