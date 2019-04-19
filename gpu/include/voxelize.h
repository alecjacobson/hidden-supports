#ifndef VOXELIZE_H
#define VOXELIZE_H

#include <Eigen/Core>
#include <igl/voxel_grid.h>
#include <igl/writeDMAT.h>
#include <igl/remove_unreferenced.h>

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
	const Eigen::VectorXi S, 
	const Eigen::MatrixXf &GV,
	const Eigen::Vector3i &side, 
	const double isovalue, 
	Eigen::MatrixXf &V_voxels, 
	Eigen::MatrixXi &F_voxels);

inline void	make_hex_from_visibility(
	const Eigen::VectorXi S, 
	const Eigen::MatrixXf &GV,
	const Eigen::Vector3i &side, 
	const int isovalue, 
	Eigen::MatrixXf &V_hex, 
	Eigen::MatrixXi &F_hex);

inline void	extract_voxel_surface(
	const Eigen::MatrixXf &V_voxels, 
	const Eigen::MatrixXi &F_voxels,
	Eigen::MatrixXf &V_quad, 
	Eigen::MatrixXi &F_quad,
	Eigen::MatrixXi &I_quad);


void removeRow(Eigen::MatrixXi& matrix, int rowToRemove)
{
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.bottomRows(numRows-rowToRemove);

    matrix.conservativeResize(numRows,numCols);
}

inline void surround_scene_in_grid(
	const float length,
	const Eigen::MatrixXf &V,
	const Eigen::MatrixXi &F,
	Eigen::Vector3i &side,
	Eigen::MatrixXf &GV)
{
	// make voxel grid

	Eigen::AlignedBox<float, 3> box(V.colwise().minCoeff(), V.colwise().maxCoeff());

	igl::voxel_grid(box, length, 3, GV, side);

	// int lowest = std::min(std::min(side(0),side(1)),side(2));
	// igl::voxel_grid(box, length, 2*(length-lowest)-1, GV, side);
}


inline void	make_voxels_from_visibility(
	const Eigen::VectorXi S, 
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

	Eigen::VectorXi S_iso = (S.array() >= (int)isovalue).select((int)isovalue, S.array()-S.array());
	// Eigen::VectorXi S_iso = S_iso_d.cast <int> ();
	int filled_count = (S_iso.array() == (int)isovalue).count();

	igl::writeDMAT("visibility_isovalue_"+std::to_string(isovalue)+".dmat", S_iso, false);

	float voxel_dims = std::abs(GV.row(0)(0) - GV.row(1)(0));

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
				if(S_iso.row(index).value() == (int)isovalue)
				{
					make_cube();

				}
	}
}

inline void	make_hex_from_visibility(
	const Eigen::VectorXi S, 
	const Eigen::MatrixXf &GV,
	const Eigen::Vector3i &side, 
	const int isovalue, 
	Eigen::MatrixXf &V_hex, 
	Eigen::MatrixXi &F_hex)
{
	int count = 0;

	int width = side(0);
    int height = side(1);
    int depth = side(2);
	// std::cout << side(0) << " * " << side(1) << " * " << side(2) << " = " << side(0)*side(1)*side(2)<< std::endl;

	Eigen::VectorXi S_iso = (S.array() >= isovalue).select(isovalue, S.array()-S.array());
	// Eigen::VectorXi S_iso = S_iso_d.cast <int> ();
	int filled_count = (S_iso.array() == isovalue).count();

	igl::writeDMAT("visibility_isovalue_"+std::to_string(isovalue)+".dmat", S_iso, false);

	float voxel_dims = std::abs(GV.row(0)(0) - GV.row(1)(0));

	// V_voxels.resize(GV.rows(),3);
	V_hex.resize(filled_count*8,3); 
	F_hex.resize(filled_count*6,4);

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
			V_hex.row(current_filled_v+0) = center + Eigen::Vector3f(-diff,-diff,-diff);
			V_hex.row(current_filled_v+1) = center + Eigen::Vector3f(diff,-diff,-diff);
			V_hex.row(current_filled_v+2) = center + Eigen::Vector3f(-diff,-diff,diff);
			V_hex.row(current_filled_v+3) = center + Eigen::Vector3f(diff,-diff,diff);
			V_hex.row(current_filled_v+4) = center + Eigen::Vector3f(-diff,diff,-diff);
			V_hex.row(current_filled_v+5) = center + Eigen::Vector3f(diff,diff,-diff);
			V_hex.row(current_filled_v+6) = center + Eigen::Vector3f(-diff,diff,diff);
			V_hex.row(current_filled_v+7) = center + Eigen::Vector3f(diff,diff,diff);

			// fill F matrix
			// bottom side (0,1,2,3)
			F_hex.row(current_filled_f+0) = Eigen::Vector4i(
												current_filled_v+0,
												current_filled_v+1,
												current_filled_v+3,
												current_filled_v+2);
			// back side (0,1,4,5)
			F_hex.row(current_filled_f+1) = Eigen::Vector4i(
												current_filled_v+0,
												current_filled_v+4,
												current_filled_v+5,
												current_filled_v+1);
			// left side (0,2,4,6)								
			F_hex.row(current_filled_f+2) = Eigen::Vector4i(
												current_filled_v+2,
												current_filled_v+6,
												current_filled_v+4,
												current_filled_v+0);

			// right side (1,3,5,7)
			F_hex.row(current_filled_f+3) = Eigen::Vector4i(
												current_filled_v+1,
												current_filled_v+5,
												current_filled_v+7,
												current_filled_v+3);


			// front side (2,3,6,7)					
			F_hex.row(current_filled_f+4) = Eigen::Vector4i(
												current_filled_v+3,
												current_filled_v+7,
												current_filled_v+6,
												current_filled_v+2);

			// top side (4,5,6,7)
			F_hex.row(current_filled_f+5) = Eigen::Vector4i(
												current_filled_v+7,
												current_filled_v+5,
												current_filled_v+4,
												current_filled_v+6);
												
			current_filled_v+=8;
			current_filled_f+=6;
			
			return 0;
		};
		// if this voxel is filled, put a cube there
		if(S_iso.row(index).value() == isovalue)
		{
			make_cube();
		}
	}

	std::cout << "how many voxels 'should have been' made " <<  filled_count << std::endl;
	std::cout << "number of voxels made " << count << std::endl;
}

inline void	extract_voxel_surface(
	const Eigen::MatrixXf &V_voxels, 
	const Eigen::MatrixXi &F_voxels,
	Eigen::MatrixXf &V_quad, 
	Eigen::MatrixXi &F_quad,
	Eigen::MatrixXi &I_quad)
{
 	Eigen::MatrixXf NV;
	Eigen::MatrixXi NF;
	Eigen::MatrixXi IM;
	Eigen::MatrixXi MI;
    igl::remove_duplicate_vertices(V_voxels, F_voxels, 1e-4, NV, IM, MI, NF);

	std::vector<int> removable_vertex_indices(NV.rows(), 0);
	for(int i = 0; i < NF.rows(); i++)
	{
		Eigen::Vector4i row = NF.row(i);
		removable_vertex_indices[row(0)] += 1;
		removable_vertex_indices[row(1)] += 1;
		removable_vertex_indices[row(2)] += 1;
		removable_vertex_indices[row(3)] += 1;
	}

	Eigen::MatrixXi F_removed;
	F_removed = NF;
	// for(int v = 0; v < removable_vertex_indices.size(); v++)
	// {
	// 	int removable_index = removable_vertex_indices[v];
	// 	std::cout << removable_index << std::endl;
		
	// }

	for(int i = 0; i < NF.rows(); i++)
	{
		Eigen::Vector4i row = NF.row(i);
		Eigen::Vector4i row_counts = Eigen::Vector4i(removable_vertex_indices[row(0)],
													removable_vertex_indices[row(1)],
													removable_vertex_indices[row(2)],
													removable_vertex_indices[row(3)]);
		if((row_counts.array() >= 24).any())
		{
			removeRow(F_removed, i);
			std::cout << "remove row" << std::endl;
		}
	}
	std::cout << F_removed.rows() << std::endl;

	std::cout << F_voxels.size() << std::endl;
	igl::remove_unreferenced(NV, F_removed, V_quad, F_quad, I_quad);
	
}

#endif
