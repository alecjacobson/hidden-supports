#ifndef GENERATE_VIEWS_H
#define GENERATE_VIEWS_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include <map>
#include <random>
#include "compare_vector3d.h"
// Given rectangle boundaries in 3D, generate a normal distribution of ray positions
//
// Inputs:
//   top_left size 3 vector position of top left of rectangle
//	 bottom_right length 3 vector position of bottom left of rectangle


void generate_views( 
	const Eigen::Vector3f &top_left,
	const Eigen::Vector3f &bottom_right,
	const float num_views,
	const std::vector<float> z_vals,
	Eigen::MatrixXf &views);
#endif
