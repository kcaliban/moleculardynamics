#ifndef TYPES_H
#define TYPES_H

#include <Eigen/Dense>

using Positions_t = Eigen::Array3Xd;
using Velocities_t = Eigen::Array3Xd;
using Accelerations_t = Eigen::Array3Xd;
using Masses_t = Eigen::VectorXd;
using Forces_t = Eigen::Array3Xd;

#endif // TYPES_H
