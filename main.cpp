#include <igl/copyleft/marching_cubes.h>
#include "generate_views.h"
#include "compute_shadow_volume.h"
#include "voxelize.h"

#include <igl/readSTL.h>
#include <igl/readOBJ.h>
#include <igl/writeSTL.h>
#include <igl/opengl/glfw/Viewer.h>

int main(int argc, char *argv[])
{
	// read in stl file
	Eigen::MatrixXf V;
	Eigen::MatrixXi F;
	Eigen::MatrixXf CN;
	std::string filename = argc > 1 ? argv[1] : "../data/example-scene.stl";
	//igl::readOBJ(filename, V, F);
	igl::readSTL(filename, V, F, CN);

	std::cout << "read frame" << std::endl;

	// compute "area light" of possible views
	Eigen::Vector3f top_left(-9.0, 8.0, 5.0);
	//Eigen::Vector3d bottom_right(-1.0, 2.0, 5.0);
	Eigen::Vector3f bottom_right(-1.0, 5.0, 5.0);

	// get view rays
	Eigen::MatrixXf views;
	generate_views(top_left, bottom_right, 100, views);
	std::cout << "generated viewing rays" << std::endl;
	for (int i = 0; i < views.rows(); i++)
		std::cout << i << ": " << views.row(i) << std::endl;

	// voxelize
	Eigen::MatrixXf GV;
	Eigen::Vector3i side;
	voxelize_scene(V, F, side, GV);
	std::cout << "made voxel grid" << std::endl;

	// compute shadow volume
	Eigen::VectorXi S;
	//compute_shadow_volume(V, F, GV, top_left, S);
	compute_shadow_volume(V, F, GV, views, S);
	std::cout << "computed shadows" << std::endl;

	Eigen::MatrixXf V_voxels;
	Eigen::MatrixXi F_voxels;
	std::cout << side << std::endl;
	igl::copyleft::marching_cubes(S, GV, side(0), side(1), side(2), 0.0, V_voxels, F_voxels);
	
	// write a new stl
	std::cout << V_voxels.rows() << std::endl;
	
	igl::writeOBJ("voxels.obj", V_voxels, F_voxels);
	igl::writeSTL("voxelized_shadow_volumes.stl", V_voxels, F_voxels);
	
	// Create a libigl Viewer object 
	
	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(V_voxels.cast <double> (), F_voxels);
	viewer.data().set_points(views.cast <double>(), Eigen::RowVector3d(0.0,0.0,0.0));
	viewer.data().set_face_based(true);
	viewer.launch();
	
	return 0;
}

