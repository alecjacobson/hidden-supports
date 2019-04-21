#ifndef READ_SCENE_H
#define READ_SCENE_H

#include "helpers.h"
#include <igl/remove_duplicate_vertices.h>

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
  sc.N.resize(0,3);
  {
    Eigen::VectorXi _1,_2;
    igl::remove_duplicate_vertices(
      Eigen::Matrix< float,Eigen::Dynamic,3,Eigen::RowMajor>(sc.V),
      Eigen::Matrix<GLuint,Eigen::Dynamic,3,Eigen::RowMajor>(sc.F),
      0,
      sc.V,
      _1,_2,sc.F);
  }


  Eigen::RowVector3f translation = sc.V.colwise().mean();
  sc.V.rowwise() -= translation;

  float scale_factor = (sc.V.colwise().maxCoeff()-sc.V.colwise().minCoeff()).maxCoeff();
  sc.V /= scale_factor;
}

#endif
