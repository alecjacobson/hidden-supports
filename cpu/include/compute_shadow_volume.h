#ifndef COMPUTE_SHADOW_VOLUME_H
#define COMPUTE_SHADOW_VOLUME_H

#include <Eigen/Core>
#include <map>
#include <igl/ray_mesh_intersect.h>
#include <igl/embree/EmbreeIntersector.h>
#include <iostream>
#include "compare_vector3d.h"

// Compute the voxelization of a boxed scene.
//
// Inputs:
//   F  #F by 3 list of triangle indices into some vertex list V
//	 V	#v by 3 list of vertex positions in the mesh

void compute_shadow_volume(const Eigen::MatrixXf &V,
	const Eigen::MatrixXi &F,
	const Eigen::MatrixXf &GV,
	const Eigen::Vector3f &view,
	Eigen::VectorXi &S);

void compute_shadow_volume(
	const Eigen::MatrixXf &V,
	const Eigen::MatrixXi &F,
	const Eigen::MatrixXf &GV,
	const Eigen::MatrixXf &views,
	Eigen::VectorXf &S);

#endif
