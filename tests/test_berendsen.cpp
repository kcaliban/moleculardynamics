#include "berendsen.h"
#include "verlet.h"
#include "temperature.h"
#include "lj_direct_summation.h"
#include <gtest/gtest.h>

// Simple Test for Berendsen thermostat
TEST(BerendsenTest, ReachesTemperatureSingleStep) {
    int nb_atoms = 10;
    double target_temperature = 300.0;

    // Initialize atoms with random positions and velocities
    Atoms atoms(nb_atoms);
    atoms.m = 1.0;
    atoms.positions.setRandom();
    atoms.velocities.setRandom();

    // Set params
    double epsilon = 1.0;
    double sigma = 1.0;
    double timestep = 1;
    int n_steps = 1;

    // Run simulation for n_steps
    for (int i = 0; i < n_steps; i++) {
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep);
        lj_direct_summation(atoms, epsilon, sigma);
        verlet_step2(atoms.velocities, atoms.forces, timestep);
        berendsen_thermostat(atoms, target_temperature, timestep, 1, true);
    }

    // Check temperature
    EXPECT_NEAR(temperature(atoms), target_temperature, 1e-1);
}

// Another simple Test for Berendsen thermostat
TEST(BerendsenTest, ReachesTemperatureMultipleSteps) {
    int nb_atoms = 10;
    double target_temperature = 300.0;

    // Initialize atoms with random positions and velocities
    Atoms atoms(nb_atoms);
    atoms.m = 1.0;
    atoms.positions.setRandom();
    atoms.velocities.setRandom();

    // Set params
    double epsilon = 1.0;
    double sigma = 1.0;
    double timestep = 1;
    int n_steps = 100;

    // Run simulation for n_steps
    for (int i = 0; i < n_steps; i++) {
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep);
        lj_direct_summation(atoms, epsilon, sigma);
        verlet_step2(atoms.velocities, atoms.forces, timestep);
        berendsen_thermostat(atoms, target_temperature, timestep, 2, true);
    }

    // Check temperature
    EXPECT_NEAR(temperature(atoms), target_temperature, 1e-1);
}
