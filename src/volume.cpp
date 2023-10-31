#include "volume.h"

double volume(const Atoms &atoms) {
    // Find the minimum and maximum coordinates in each direction
    Eigen::Vector3d min{atoms.positions.rowwise().minCoeff()};
    Eigen::Vector3d max{atoms.positions.rowwise().maxCoeff()};

    // Compute the volume of the box
    return (max - min).prod();
}