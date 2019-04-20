#ifndef READ_SCENE_H
#define READ_SCENE_H

#include "helpers.h"

inline void read_scene(
  scene &sc,
  const std::string &filename
);

inline void read_scene(
  scene &sc,
  const std::string &filename
)
{
  // Read input scene from file
  igl::readSTL(filename, sc.V, sc.F, sc.N);

  Eigen::RowVector3f translation = sc.V.colwise().mean();
  sc.V.rowwise() -= translation;

  float scale_factor = (sc.V.colwise().maxCoeff()-sc.V.colwise().minCoeff()).maxCoeff();
  sc.V /= scale_factor;
}

#endif