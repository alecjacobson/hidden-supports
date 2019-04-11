#ifndef VOXELIZE_H
#define VOXELIZE_H

#include <Eigen/Core>
#include <igl/voxel_grid.h>
#include <igl/writeDMAT.h>

// #include "cc3d.hpp"

// // 3d array represented as 1d array
// int* labels = new int[512*512*512](); 

// int* cc_labels = cc3d::connected_components3d<int>(
//   labels, /*sx=*/512, /*sy=*/512, /*sz=*/512
// );

// Compute the voxelization of a boxed scene.
//
// Inputs:
//	 length: length of the longest side of the voxel grid
//   F  #F by 3 list of triangle indices into some vertex list V
//	 V	#v by 3 list of vertex positions in the mesh
inline void surround_scene_in_grid(
	const float length,
	const Eigen::MatrixXf &V, 
	const Eigen::MatrixXi &F, 
	Eigen::Vector3i &side,
	Eigen::MatrixXf &GV);

inline void	make_voxels_from_visibility(
	const Eigen::VectorXf S, 
	const Eigen::MatrixXf &GV,
	const Eigen::Vector3i &side, 
	const double isovalue, 
	Eigen::MatrixXf &V_voxels, 
	Eigen::MatrixXi &F_voxels);

inline void	make_hex_from_visibility(
	const Eigen::VectorXf &S, 
	const Eigen::MatrixXf &GV,
	const Eigen::Vector3i &side, 
	const double isovalue, 
	Eigen::MatrixXf &V_hex, 
	Eigen::MatrixXi &F_hex);

inline void surround_scene_in_grid(
	const float length,
	const Eigen::MatrixXf &V,
	const Eigen::MatrixXi &F,
	Eigen::Vector3i &side,
	Eigen::MatrixXf &GV)
{
	// make voxel grid

	Eigen::AlignedBox<float, 3> box(V.colwise().minCoeff(), V.colwise().maxCoeff());

	igl::voxel_grid(box, length, 1, GV, side);

	// int lowest = std::min(std::min(side(0),side(1)),side(2));
	// igl::voxel_grid(box, length, 2*(length-lowest)-1, GV, side);
}


inline void	make_voxels_from_visibility(
	const Eigen::VectorXf S, 
	const Eigen::MatrixXf &GV,
	const Eigen::Vector3i &side, 
	const double isovalue, 
	Eigen::MatrixXf &V_voxels, 
	Eigen::MatrixXi &F_voxels)
{
	int count = 0;

	int width = side(0);
    int height = side(1);
    int depth = side(2);
	// std::cout << side(0) << " * " << side(1) << " * " << side(2) << " = " << side(0)*side(1)*side(2)<< std::endl;

	Eigen::VectorXf S_iso = (S.array() >= isovalue).select(isovalue, S.array()-S.array());
	int filled_count = (S_iso.array() == isovalue).count();

	igl::writeDMAT("visibility_isovalue_"+std::to_string(isovalue)+".dmat", S_iso, false);

	float voxel_dims = std::abs(GV.row(0)(0) - GV.row(1)(0));
	//std::cout << voxel_dims << std::endl;

	// V_voxels.resize(GV.rows(),3);
	V_voxels.resize(filled_count*8,3); 
	F_voxels.resize(filled_count*12,3);

	int current_filled_v = 0;
	int current_filled_f = 0;

	for(int index = 0; index < GV.rows(); index ++)
	{
				// std::cout << index << std::endl;
				Eigen::VectorXf center = GV.row(index);
				auto make_cube = [&](){
					count ++;
					float diff = voxel_dims/2;
					// V_voxels.row(index) = center;

					// fill V matrix
					V_voxels.row(current_filled_v+0) = center + Eigen::Vector3f(-diff,-diff,-diff);
					V_voxels.row(current_filled_v+1) = center + Eigen::Vector3f(diff,-diff,-diff);
					V_voxels.row(current_filled_v+2) = center + Eigen::Vector3f(-diff,-diff,diff);
					V_voxels.row(current_filled_v+3) = center + Eigen::Vector3f(diff,-diff,diff);
					V_voxels.row(current_filled_v+4) = center + Eigen::Vector3f(-diff,diff,-diff);
					V_voxels.row(current_filled_v+5) = center + Eigen::Vector3f(diff,diff,-diff);
					V_voxels.row(current_filled_v+6) = center + Eigen::Vector3f(-diff,diff,diff);
					V_voxels.row(current_filled_v+7) = center + Eigen::Vector3f(diff,diff,diff);

					// fill F matrix
					// bottom (0,1,2,3)
					F_voxels.row(current_filled_f+0) = Eigen::Vector3i(
														current_filled_v+0,
														current_filled_v+1,
														current_filled_v+2);
					F_voxels.row(current_filled_f+1) = Eigen::Vector3i(
														current_filled_v+2,
														current_filled_v+1,
														current_filled_v+3);
					// left side (0,2,4,6)								
					F_voxels.row(current_filled_f+2) = Eigen::Vector3i(
														current_filled_v+2,
														current_filled_v+4,
														current_filled_v+0);
					F_voxels.row(current_filled_f+3) = Eigen::Vector3i(
														current_filled_v+2,
														current_filled_v+6,
														current_filled_v+4);
					// front side (2,3,6,7)								
					F_voxels.row(current_filled_f+4) = Eigen::Vector3i(
														current_filled_v+2,
														current_filled_v+7,
														current_filled_v+6);
					F_voxels.row(current_filled_f+5) = Eigen::Vector3i(
														current_filled_v+7,
														current_filled_v+2,
														current_filled_v+3);
					// right side (1,3,5,7)						
					F_voxels.row(current_filled_f+6) = Eigen::Vector3i(
														current_filled_v+3,
														current_filled_v+5,
														current_filled_v+7);
					F_voxels.row(current_filled_f+7) = Eigen::Vector3i(
														current_filled_v+3,
														current_filled_v+1,
														current_filled_v+5);
					// back side (0,1,4,5)							
					F_voxels.row(current_filled_f+8) = Eigen::Vector3i(
														current_filled_v+1,
														current_filled_v+0,
														current_filled_v+5);
					F_voxels.row(current_filled_f+9) = Eigen::Vector3i(
														current_filled_v+0,
														current_filled_v+4,
														current_filled_v+5);
					// top side (4,5,6,7)
					F_voxels.row(current_filled_f+10) = Eigen::Vector3i(
														current_filled_v+6,
														current_filled_v+7,
														current_filled_v+5);
					F_voxels.row(current_filled_f+11) = Eigen::Vector3i(
														current_filled_v+6,
														current_filled_v+5,
														current_filled_v+4);
					
														
					current_filled_v+=8;
					current_filled_f+=12;
					
					return 0;
				};
				// if this voxel is filled, put a cube there
				if(S_iso.row(index).value() == isovalue)
				{
					make_cube();

				}
	}
	std::cout << "number of voxels made " << count << std::endl;
}

inline void	make_hex_from_visibility(
	const Eigen::VectorXf &S, 
	const Eigen::MatrixXf &GV,
	const Eigen::Vector3i &side, 
	const double isovalue, 
	Eigen::MatrixXf &V_hex, 
	Eigen::MatrixXi &F_hex)
{
	int count = 0;

	int width = side(0);
    int height = side(1);
    int depth = side(2);
	std::cout << side(0) << " * " << side(1) << " * " << side(2) << " = " << side(0)*side(1)*side(2)<< std::endl;

	Eigen::VectorXf S_iso = (S.array() >= isovalue).select(1, S.array()-S.array());
	int filled_count = (S_iso.array() == 1).count();

	// igl::writeDMAT("visibility_isovalue_"+std::to_string(isovalue)+".dmat", S_iso, false);

	float voxel_dims = std::abs(GV.row(0)(0) - GV.row(1)(0));
	std::cout << voxel_dims << std::endl;

	// V_voxels.resize(GV.rows(),3);
	V_hex.resize(filled_count*8,3); 
	F_hex.resize(filled_count*12,3);

	int current_filled_v = 0;
	int current_filled_f = 0;

	for(int index = 0; index < GV.rows(); index ++)
	{
				std::cout << index << std::endl;
				Eigen::VectorXf center = GV.row(index);
				auto make_cube = [&](){
					count ++;
					float diff = voxel_dims/2;
					// V_hex.row(index) = center;

					// fill V matrix
					V_hex.row(current_filled_v+0) = center + Eigen::Vector3f(-diff,-diff,-diff);
					V_hex.row(current_filled_v+1) = center + Eigen::Vector3f(diff,-diff,-diff);
					V_hex.row(current_filled_v+2) = center + Eigen::Vector3f(-diff,-diff,diff);
					V_hex.row(current_filled_v+3) = center + Eigen::Vector3f(diff,-diff,diff);
					V_hex.row(current_filled_v+4) = center + Eigen::Vector3f(-diff,diff,-diff);
					V_hex.row(current_filled_v+5) = center + Eigen::Vector3f(diff,diff,-diff);
					V_hex.row(current_filled_v+6) = center + Eigen::Vector3f(-diff,diff,diff);
					V_hex.row(current_filled_v+7) = center + Eigen::Vector3f(diff,diff,diff);

					// fill F matrix
					// bottom (0,1,2,3)
					F_hex.row(current_filled_f+0) = Eigen::Vector3i(
														current_filled_v+0,
														current_filled_v+1,
														current_filled_v+2);
					F_hex.row(current_filled_f+1) = Eigen::Vector3i(
														current_filled_v+2,
														current_filled_v+1,
														current_filled_v+3);
					// left side (0,2,4,6)								
					F_hex.row(current_filled_f+2) = Eigen::Vector3i(
														current_filled_v+2,
														current_filled_v+4,
														current_filled_v+0);
					F_hex.row(current_filled_f+3) = Eigen::Vector3i(
														current_filled_v+2,
														current_filled_v+6,
														current_filled_v+4);
					// front side (2,3,6,7)								
					F_hex.row(current_filled_f+4) = Eigen::Vector3i(
														current_filled_v+2,
														current_filled_v+7,
														current_filled_v+6);
					F_hex.row(current_filled_f+5) = Eigen::Vector3i(
														current_filled_v+7,
														current_filled_v+2,
														current_filled_v+3);
					// right side (1,3,5,7)						
					F_hex.row(current_filled_f+6) = Eigen::Vector3i(
														current_filled_v+3,
														current_filled_v+5,
														current_filled_v+7);
					F_hex.row(current_filled_f+7) = Eigen::Vector3i(
														current_filled_v+3,
														current_filled_v+1,
														current_filled_v+5);
					// back side (0,1,4,5)							
					F_hex.row(current_filled_f+8) = Eigen::Vector3i(
														current_filled_v+1,
														current_filled_v+0,
														current_filled_v+5);
					F_hex.row(current_filled_f+9) = Eigen::Vector3i(
														current_filled_v+0,
														current_filled_v+4,
														current_filled_v+5);
					// top side (4,5,6,7)
					F_hex.row(current_filled_f+10) = Eigen::Vector3i(
														current_filled_v+6,
														current_filled_v+7,
														current_filled_v+5);
					F_hex.row(current_filled_f+11) = Eigen::Vector3i(
														current_filled_v+6,
														current_filled_v+5,
														current_filled_v+4);
					
														
					current_filled_v+=8;
					current_filled_f+=12;
					
					return 0;
				};
				// if this voxel is filled, put a cube there
				if(S_iso.row(index).value() == 1)
				{
					make_cube();

				}
	}
}

#endif