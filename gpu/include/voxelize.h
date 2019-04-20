#ifndef VOXELIZE_H
#define VOXELIZE_H

#include <Eigen/Core>
#include <igl/voxel_grid.h>
#include <igl/writeDMAT.h>
#include <igl/remove_unreferenced.h>

// Compute the voxelization of a boxed scene.
//
// Inputs:
//	 length: length of the longest side of the voxel grid
//   F  #F by 3 list of triangle indices into some vertex list V
//	 V	#v by 3 list of vertex positions in the mesh
inline void surround_scene_in_grid(
	const float length,
	const Eigen::MatrixXf &V, 
	const Eigen::MatrixXi &F, 
	Eigen::Vector3i &side,
	Eigen::MatrixXf &GV);


inline void surround_scene_in_grid(
	const float length,
	const Eigen::MatrixXf &V,
	const Eigen::MatrixXi &F,
	Eigen::Vector3i &side,
	Eigen::MatrixXf &GV)
{
	// make voxel grid

	Eigen::AlignedBox<float, 3> box(V.colwise().minCoeff(), V.colwise().maxCoeff());

	igl::voxel_grid(box, length, 3, GV, side);

	// int lowest = std::min(std::min(side(0),side(1)),side(2));
	// igl::voxel_grid(box, length, 2*(length-lowest)-1, GV, side);
}

#endif
