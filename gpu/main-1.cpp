#include <string>
#include <chrono>
#include <thread>
#include <iostream>
#include <climits>
#include <cstdlib>

#include "helpers.h"
#include "voxelize.h"
#include "generate_views.h"
#include "read_scene.h"
#include "calculate_visibilities.h"

int main(int argc, char * argv[])
{
  std::vector<std::string> paths = {"../data/shaders.json", 
                            "../data/quad_shaders.json", 
                            "../data/render.json"};

  std::string filename = argv[1];

  float num_views = std::atof (argv[2]);

  int resolution = 200;

  scene sc;
  read_scene(sc, filename);
  std::cout << "read scene" << std::endl;

  // make voxel grid
  Eigen::MatrixXf GV;
  Eigen::Vector3i side;
  Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor>  F_int = sc.F.cast <int> ();
  std::cout << "casted to int" << std::endl;
  surround_scene_in_grid(resolution, sc.V, F_int, side, GV);

  std::cout << "made voxel grid" << std::endl;

  // get view rays
  Eigen::Vector3f bottom_left = sc.V.colwise().minCoeff();
  // Eigen::Vector3f bottom_left(-0.5,0,-1);  
  Eigen::Vector3f top_right = sc.V.colwise().maxCoeff();
  // Eigen::Vector3f top_right(0.5,0.4,-1);
  Eigen::MatrixXf views;
  std::vector<float> z_vals{-1};//, -0.7, -1.5, -1.3};
  generate_views_on_plane(bottom_left, top_right, num_views, z_vals, views);

	std::cout << "generated viewing rays" << std::endl;
  std::cout << "views: " << std::endl;
  for(int i = 0; i < views.rows(); i++)
  {
    std::cout << views.row(i) << std::endl;
  }
  // igl::writeDMAT("../results/views.dmat", views, true);

  Eigen::VectorXi S;
  int ok = calculate_visibilities(side, views, paths, sc, S);    
  if(!ok)
  {
    std::cout << "not ok, exiting" << std::endl;
    return 0;
  }
  else
  {
    std::cout << "calculated visibilities" << std::endl;
  }
  
  std::cout << "size of S in main " << S.rows() << std::endl;
  igl::writeDMAT("../results/visibilities.dmat", S, true);

  std::cout << side << std::endl;

  for(int isovalue = 1; isovalue <= views.rows(); isovalue+=10)
  {
    std::cout << "isovalue " << isovalue << std::endl;
    std::stringstream ss;
    ss << std::setfill('0') << std::setw(3) << isovalue;

    // voxelize
    Eigen::VectorXi S_iso = (S.array() >= isovalue).select(isovalue, S.array()-S.array());
    igl::writeDMAT("../results/visibility_isovalue_"+std::to_string(isovalue)+".dmat", S_iso, false);

  }

// // transform back meshes
// Eigen::Vector3f m = V.colwise().maxCoeff();
// V *= scale_factor;
// V.rowwise() += translation;

// // std::cout << (V.colwise().maxCoeff()-V.colwise().minCoeff()).maxCoeff() << std::endl;
// V_voxels *= scale_factor;
// V_voxels.rowwise() += translation;

// GV *= scale_factor;
// GV.rowwise() += translation;

// views *= scale_factor;
// views.rowwise() += translation;
// // views.rowwise() += Eigen::RowVector3f(0,0,-translation(2)+scale_factor);

// // Create a libigl Viewer object
// igl::opengl::glfw::Viewer viewer;
// viewer.data().set_mesh(V_voxels.cast <double> (), F_voxels);
// viewer.data().set_face_based(true);
// viewer.data().set_colors(Eigen::RowVector3d(146,197,222)/255.);
// // viewer.append_mesh();
// // viewer.data().set_mesh(V.cast <double> (), F_int);
// // viewer.data().set_face_based(true);
// // viewer.data().set_colors(Eigen::RowVector3d(244,165,130)/255.);

// // viewer.data().set_points(GV.cast <double>(), Eigen::RowVector3d(0.0,0.0,0.0));
// // viewer.data().point_size = 10;

// viewer.data().set_points(views.cast <double>(), Eigen::RowVector3d(0.0,0.0,0.0));
// viewer.data().point_size = 10;
// viewer.launch();

  return EXIT_SUCCESS;
}
