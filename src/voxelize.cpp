#include "voxelize.h"

void voxelize_scene(const Eigen::MatrixXf &V,
	const Eigen::MatrixXi &F,
	Eigen::Vector3i &side,
	Eigen::MatrixXf &GV)
{
	// make voxel grid

	Eigen::AlignedBox<float, 3> box(V.colwise().minCoeff(), V.colwise().maxCoeff());

	//std::cout << V.colwise().minCoeff() << std::endl << V.colwise().maxCoeff() << std::endl;

	igl::voxel_grid(box, 100, 1, GV, side);
}