#include "virial.h"


Eigen::Matrix3d virial(const Atoms &atoms, double volume, const NeighborList &neighbor_list) {
    // Calculate the embedding energies
    Eigen::ArrayXd embeddings = ducastelle_embedding(atoms, neighbor_list);

    // Calculate the virial stress
    Eigen::Matrix3d virial_stress = Eigen::Matrix3d::Zero();
    for (auto [i, j] : neighbor_list) {
        if (i < j) {
            Eigen::Vector3d r = atoms.positions.col(j) - atoms.positions.col(i);
            Eigen::Vector3d f = ducastelle_force(atoms, embeddings, i, j);
            Eigen::Matrix3d prod = r * f.transpose();
            virial_stress += prod;
        }
    }

    // ((eV / Å) * Å) / Å^3 = eV / Å^3
    return virial_stress / volume;
}


std::pair<Eigen::Matrix3d,Eigen::Matrix3d> virial(const Atoms &atoms, double volume, const NeighborList &neighbor_list, const Eigen::ArrayXd &embeddings, int nb_local) {
    Eigen::Matrix3d virial_stress_local = Eigen::Matrix3d::Zero();
    // Stress for local atoms
    for (auto [i, j] : neighbor_list) {
        if (i < j && j < nb_local) {
            Eigen::Vector3d r = atoms.positions.col(j) - atoms.positions.col(i);
            Eigen::Vector3d f = ducastelle_force(atoms, embeddings, i, j);
            Eigen::Matrix3d prod = r * f.transpose();
            virial_stress_local += prod;
        }
    }
    // Stress for ghost atoms
    Eigen::Matrix3d virial_stress_ghost = Eigen::Matrix3d::Zero();
    for (auto [i, j] : neighbor_list) {
        if (i < j && i < nb_local && j >= nb_local) {
            Eigen::Vector3d r = atoms.positions.col(j) - atoms.positions.col(i);
            Eigen::Vector3d f = ducastelle_force(atoms, embeddings, i, j);
            Eigen::Matrix3d prod = r * f.transpose();
            virial_stress_ghost += prod;
        }
    }
    // ((eV / Å) * Å) / Å^3 = eV / Å^3
    return {virial_stress_local / volume, virial_stress_ghost / volume };
}

double stress_z(const Atoms &atoms, double volume, const NeighborList &neighbor_list) {
    Eigen::Matrix3d virial_stress = virial(atoms, volume, neighbor_list);
    return virial_stress.coeff(2, 2);
}

std::pair<double,double> stress_z(const Atoms &atoms, double volume, const NeighborList &neighbor_list, const Eigen::ArrayXd &embeddings, int nb_local) {
    auto [virial_stress_local, virial_stress_ghost] = virial(atoms, volume, neighbor_list, embeddings, nb_local);
    return {virial_stress_local.coeff(2, 2), virial_stress_ghost.coeff(2, 2)};
}
