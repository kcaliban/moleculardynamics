#include "virial.h"
#include "ducastelle.h"
#include "types.h"
#include "neighbors.h"
#include <gtest/gtest.h>
#include <cmath>

// Test of calculation of virial stress (for regression)
TEST(VirialStressTest, SimpleTest) {
    Positions_t positions(3, 3);
    Velocities_t velocities(3, 3);
    Forces_t forces(3, 3);

    // Set initial conditions
    positions(0, 0) = 0; positions(1, 0) = 0; positions(2, 0) = 0;
    velocities(0, 0) = 1; velocities(1, 0) = 2; velocities(2, 0) = -2;

    positions(0, 1) = 1; positions(1, 1) = 3; positions(2, 1) = -2;
    velocities(0, 1) = 1; velocities(1, 1) = -1; velocities(2, 1) = 2;

    positions(0, 2) = 2; positions(1, 2) = 3; positions(2, 2) = -5;
    velocities(0, 2) = -1; velocities(1, 2) = 2; velocities(2, 2) = -3;

    Atoms atoms(positions, velocities);
    auto volume = 10;

    NeighborList neighborList;
    neighborList.update(atoms, 5);

    // Calculate virial stress
    auto stress_tensor = virial(atoms, volume, neighborList);

    // Check
    ASSERT_NEAR(stress_tensor.sum(), 0.447922, 1e-6);
}

// Test to make sure results of separated between ghosts and local atoms is the same
TEST(VirialStressTest, GhostAtomsTest) {
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

    Atoms atoms(positions, velocities);
    atoms.forces = forces;
    auto volume = 10;

    NeighborList neighborList;
    neighborList.update(atoms, 5);

    // Calculate ducastelle embeddings
    auto embeddings = ducastelle_embedding(atoms, neighborList);

    // Calculate virial stress from local and ghost atoms
    auto stress_tensor = virial(atoms, volume, neighborList);
    auto [stress_tensor_local, stress_tensor_ghosts] = virial(atoms, volume, neighborList, embeddings, 2);

    // Check
    ASSERT_NEAR(stress_tensor.sum(), stress_tensor_local.sum() + stress_tensor_ghosts.sum(), 1e-6);
}

// Test to make sure the right component is extracted when using stress_z (for regression)
TEST(VirialStressTest, StressZTest) {
        Positions_t positions(3, 3);
    Velocities_t velocities(3, 3);
    Forces_t forces(3, 3);

    // Set initial conditions
    positions(0, 0) = 0; positions(1, 0) = 0; positions(2, 0) = 0;
    velocities(0, 0) = 1; velocities(1, 0) = 2; velocities(2, 0) = -2;

    positions(0, 1) = 1; positions(1, 1) = 3; positions(2, 1) = -2;
    velocities(0, 1) = 1; velocities(1, 1) = -1; velocities(2, 1) = 2;

    positions(0, 2) = 2; positions(1, 2) = 3; positions(2, 2) = -5;
    velocities(0, 2) = -1; velocities(1, 2) = 2; velocities(2, 2) = -3;

    Atoms atoms(positions, velocities);
    auto volume = 10;

    NeighborList neighborList;
    neighborList.update(atoms, 5);

    // Calculate virial stress and only z component
    auto stress_tensor = virial(atoms, volume, neighborList);
    auto stressz = stress_z(atoms, volume, neighborList);

    // Check
    ASSERT_NEAR(stressz, stress_tensor.coeff(2, 2), 1e-6);
}

// Test to make sure results of separated between ghosts and local atoms is the same for z component (regression)
TEST(VirialStressTest, GhostAtomsStressZTest) {
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

    Atoms atoms(positions, velocities);
    atoms.forces = forces;
    auto volume = 10;

    NeighborList neighborList;
    neighborList.update(atoms, 5);

    // Calculate ducastelle embeddings
    auto embeddings = ducastelle_embedding(atoms, neighborList);

    // Calculate virial stress and only z component, for local and ghosts
    auto stress_tensor = virial(atoms, volume, neighborList);
    auto [stress_tensor_local, stress_tensor_ghosts] = virial(atoms, volume, neighborList, embeddings, 2);
    auto stressz = stress_z(atoms, volume, neighborList);
    auto [stressz_local, stressz_ghosts] = stress_z(atoms, volume, neighborList, embeddings, 2);

    // Check
    ASSERT_NEAR(stressz, stress_tensor.coeff(2, 2), 1e-6);
    ASSERT_NEAR(stressz_local, stress_tensor_local.coeff(2, 2), 1e-6);
    ASSERT_NEAR(stressz_ghosts, stress_tensor_ghosts.coeff(2, 2), 1e-6);
    ASSERT_NEAR(stressz, stressz_local + stressz_ghosts, 1e-6);
}