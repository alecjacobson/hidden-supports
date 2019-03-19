#ifndef VOXELIZE_H
#define VOXELIZE_H

#include <Eigen/Core>
#include <igl/voxel_grid.h>
#include <igl/winding_number.h>
#include "compare_vector3d.h"

// Compute the voxelization of a boxed scene.
//
// Inputs:
//   F  #F by 3 list of triangle indices into some vertex list V
//	 V	#v by 3 list of vertex positions in the mesh
void voxelize_scene(const Eigen::MatrixXf &V, 
	const Eigen::MatrixXi &F, 
	Eigen::Vector3i &side,
	Eigen::MatrixXf &GV);


void make_voxels_from_visibility(
	const Eigen::VectorXf &S, 
	const Eigen::MatrixXf &GV,
	const Eigen::Vector3i &side, 
	const int isovalue, 
	Eigen::MatrixXf &V_voxels, 
	Eigen::MatrixXi &F_voxels);

#endif
