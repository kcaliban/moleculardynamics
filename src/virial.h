#ifndef VIRIAL_H
#define VIRIAL_H

#include <Eigen/Dense>
#include "atoms.h"
#include "ducastelle.h"

/**
 * Calculate the virial stress tensor for given atoms
 * 
 * @param atoms Atoms to calculate the virial stress tensor for
 * @param volume Volume of the simulation box
 * @param neighbor_list Neighbor list
 * @return Virial stress tensor
*/
Eigen::Matrix3d virial(const Atoms &atoms, double volume, const NeighborList &neighbor_list);

/**
 * Calculate the local and ghost virial stress tensor for given atoms and embeddings
 * 
 * @param atoms Atoms to calculate the virial stress tensor for
 * @param volume Volume of the simulation box
 * @param neighbor_list Neighbor list
 * @param embeddings Embedding energies
 * @param nb_local Number of local atoms
*/
std::pair<Eigen::Matrix3d,Eigen::Matrix3d> virial(const Atoms &atoms, double volume, const NeighborList &neighbor_list, const Eigen::ArrayXd &embeddings, int nb_local);

/**
 * Calculate the stress in z-direction for given atoms
 * 
 * @param atoms Atoms to calculate the stress for
 * @param volume Volume of the simulation box
 * @param neighbor_list Neighbor list
 * @return Stress in z-direction
*/
double stress_z(const Atoms &atoms, double volume, const NeighborList &neighbor_list);

/**
 * Calculate the local and ghost stress in z-direction for given atoms and embeddings
 * 
 * @param atoms Atoms to calculate the stress for
 * @param volume Volume of the simulation box
 * @param neighbor_list Neighbor list
 * @param embeddings Embedding energies
 * @param nb_local Number of local atoms
*/
std::pair<double,double> stress_z(const Atoms &atoms, double volume, const NeighborList &neighbor_list, const Eigen::ArrayXd &embeddings, int nb_local);

#endif // VIRIAL_H