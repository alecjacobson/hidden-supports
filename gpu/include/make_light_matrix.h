#ifndef MAKE_LIGHT_MATRIX_H
#define MAKE_LIGHT_MATRIX_H

#include <Eigen/Core>
#include "helpers.h"

inline void make_light_matrix(
	const scene &sc,
	const Eigen::Vector3f &viewpoint,
  const Eigen::Vector3f &centroid,
  Eigen::Affine3f &light_view);

inline void make_light_matrix(
	const Eigen::Vector3f &viewpoint,
  const Eigen::Vector3f &centroid,
  Eigen::Affine3f &light_view)
{
  Eigen::Vector3f l = centroid - viewpoint;
  l.normalize();
  Eigen::Vector3f n = -(viewpoint - Eigen::Vector3f(0,0,1));
  n.normalize();
  float cos_spin_angle = n.dot(l);
  // if x value is negative, then rotate positive angle
  float spin_angle = (viewpoint(0) < 0) ? acos(cos_spin_angle) : -1*acos(cos_spin_angle);

  light_view = Eigen::Affine3f::Identity() * 
      Eigen::Translation3f(viewpoint);
  light_view.rotate(
      Eigen::AngleAxisf(
      spin_angle,
      Eigen::Vector3f(0,1,0)));
}

#endif