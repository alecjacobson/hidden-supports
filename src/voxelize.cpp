#include "voxelize.h"

void voxelize_scene(const Eigen::MatrixXf &V,
	const Eigen::MatrixXi &F,
	Eigen::Vector3i &side,
	Eigen::MatrixXf &GV)
{
	// make voxel grid

	Eigen::AlignedBox<float, 3> box(V.colwise().minCoeff(), V.colwise().maxCoeff());

	//std::cout << V.colwise().minCoeff() << std::endl << V.colwise().maxCoeff() << std::endl;

	igl::voxel_grid(box, 64, 1, GV, side);
}


void make_voxels_from_visibility(
	const Eigen::VectorXf &S, 
	const Eigen::MatrixXf &GV,
	const Eigen::Vector3i &side, 
	const int isovalue, 
	Eigen::MatrixXf &V_voxels, 
	Eigen::MatrixXi &F_voxels)
{
	int width = side(0);
    int height = side(1);
    int depth = side(2);

	Eigen::VectorXf S_iso = (S.array() <= isovalue).select(-1, S);
	int filled_count = (S_iso.array() == -1).count();
	// std::cout << "number of filled voxels " << filled_count << std::endl;

	float voxel_dims = std::abs(GV.row(0)(0) - GV.row(1)(0));
	//std::cout << voxel_dims << std::endl;

	V_voxels.resize(filled_count*8,3);
	F_voxels.resize(filled_count*12,3);

	int current_filled_v = 0;
	int current_filled_f = 0;
	for(int z = 0; z < depth; z++)    
	{
        for(int y = 0; y < height; y++)
        {
            for(int x = 0; x < width; x++)
            {
				int index = y*width + z*width*height + x;
				Eigen::VectorXf center = GV.row(index);
				auto make_cube = [&](){
					float diff = voxel_dims/2;

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
					
														
					
					return 0;
				};
				// if this voxel is filled, put a cube there
				if(S_iso.row(index).value() != -1)
				{
					make_cube();

				}
				// increment, if possible
				if (current_filled_v+8 < index && current_filled_f < index)
				{
					current_filled_v+=8;
					current_filled_f+=12;
				}
			}
		}
	}
}

