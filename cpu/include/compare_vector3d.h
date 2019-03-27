#ifndef COMPARE_VECTOR3D_H
#define COMPARE_VECTOR3D_H

#include <Eigen/Core>

struct comp_eigen_vector3d {
	bool operator()(const Eigen::Vector3d &lhs, const Eigen::Vector3d &rhs) const {
		return lhs.norm() < rhs.norm();
		//return lhs.x() < rhs.x() && lhs.y() < rhs.y() && lhs.z() < rhs.z();
	}
};

#endif
