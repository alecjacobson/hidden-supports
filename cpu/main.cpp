#include <igl/copyleft/marching_cubes.h>
#include "generate_views.h"
#include "compute_shadow_volume.h"
#include "voxelize.h"

#include <igl/readSTL.h>
#include <igl/readOBJ.h>
#include <igl/writeSTL.h>
#include <igl/writeOBJ.h>
#include <igl/opengl/glfw/Viewer.h>

int main(int argc, char *argv[])
{
	float num_views = 100;

	// compute "area light" of possible views
	Eigen::Vector3f top_left(-8.5, 7.5, 3.5);
	//Eigen::Vector3d bottom_right(-1.0, 2.0, 5.0);
	Eigen::Vector3f bottom_right(0.5, 0.0, 3.5);

	// get view rays
	Eigen::MatrixXf views;
	generate_views(top_left, bottom_right, num_views, views);
	std::cout << "generated viewing rays" << std::endl;
	//for (int i = 0; i < views.rows(); i++)
	//	std::cout << i << ": " << views.row(i) << std::endl;

	// read in stl file(s)
	// for(int frame = 0; frame < 12; frame++)
	// {
		int frame = 0;
		////
		Eigen::MatrixXf V;
		Eigen::MatrixXi F;
		Eigen::MatrixXf CN;
		std::stringstream ss;
		ss << std::setw(2) << std::setfill('0') << frame;
		std::string n = ss.str();
		std::string filename = argc > 1 ? argv[1] : "../data/canonical-scene-frame-00.stl";
		if(argc > 1)
			filename = filename + "-" + n + ".stl";
		std::cout << "filename " << filename << std::endl;
		igl::readSTL(filename, V, F, CN);

		std::cout << "read frame" << std::endl;

		////

		// make voxel grid
		Eigen::MatrixXf GV;
		Eigen::Vector3i side;
		voxelize_scene(V, F, side, GV);
		std::cout << "made voxel grid" << std::endl;

		// compute shadow volume
		Eigen::VectorXf S;
		//compute_shadow_volume(V, F, GV, top_left, S);
		compute_shadow_volume(V, F, GV, views, S);
		std::cout << "computed shadows" << std::endl;

		////

		for (int i = 0; i < S.rows(); i++)
		{
			// if(!(S.row(i).value() == 1))
			// 	if(!(S.row(i).value() == -1))
					std::cout << i << ": " << S.row(i) << std::endl;
		}

		Eigen::MatrixXf V_voxels;
		Eigen::MatrixXi F_voxels;
		std::cout << side << std::endl;

		////

		// voxelize

		igl::copyleft::marching_cubes(S, GV, side(0), side(1), side(2), -1.0, V_voxels, F_voxels);
		// write a new stl
		std::string output_file = "../data/voxelized-canonical-frame-"+n+".stl";
		std::cout << "writing file " << n << std::endl << "-----" << std::endl;
		igl::writeSTL(output_file, V_voxels, F_voxels);

		make_voxels_from_visibility(S, GV,side, -1, V_voxels, F_voxels);
		igl::writeOBJ("test.obj", V_voxels, F_voxels);

		// for(int i = 1; i <= 5; i++)
		// {
		// 	double iso = (double)i/5;
		// 	igl::copyleft::marching_cubes(S, GV, side(0), side(1), side(2), -iso, V_voxels, F_voxels);
		// 	// write a new stl
		// 	std::cout << "isovalue " << iso << std::endl;

		// 	std::stringstream ss;
		// 	ss << std::setw(3) << std::setfill('0') << i;
		// 	std::string n = ss.str();
		// 	std::string output_file = "../data/voxelized-canonical-frame-"+n+".stl";
		// 	std::cout << "writing file " << n << std::endl << "-----" << std::endl;
		// 	igl::writeSTL(output_file, V_voxels, F_voxels);

		// }

	// }
	
	/*
	// Create a libigl Viewer object 
	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(V_voxels.cast <double> (), F_voxels);
	viewer.data().set_points(views.cast <double>(), Eigen::RowVector3d(0.0,0.0,0.0));
	viewer.data().set_face_based(true);
	viewer.launch();
	*/

	return 0;
}

