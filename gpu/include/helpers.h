#ifndef HELPERS_H
#define HELPERS_H

// make sure the modern opengl headers are included before any others
#include "gl.h"

// #include <OpenGL/gl3.h>
// #define __gl_h_
#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>

#include <vector>
#include <string>

#include <igl/readSTL.h>
#include <igl/writeSTL.h>
#include <igl/writeOBJ.h>
#include <igl/writeDMAT.h>
#include <igl/get_seconds.h>
#include <igl/copyleft/marching_cubes.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/remove_unreferenced.h>
#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/centroid.h>

struct shader_path 
{
  std::vector<std::string> vertex_paths;
  std::vector<std::string> fragment_paths;
};

struct scene
{
  Eigen::Matrix< float,Eigen::Dynamic,3,Eigen::RowMajor> V;
  Eigen::Matrix<GLuint,Eigen::Dynamic,3,Eigen::RowMajor> F;
  Eigen::Matrix< float,Eigen::Dynamic,3,Eigen::RowMajor> N;
  Eigen::Matrix< float,Eigen::Dynamic,2,Eigen::RowMajor> TC;
};

#endif
