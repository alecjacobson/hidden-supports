#include "compute_shadow_volume.h"

void compute_shadow_volume(const Eigen::MatrixXf &V,
	const Eigen::MatrixXi &F,
	const Eigen::MatrixXf &GV,
	const Eigen::Vector3f &view,
	Eigen::VectorXi &S)
{
	// intersect ray with boxed mesh
	bool was_hit;
	igl::Hit hit;

	igl::embree::EmbreeIntersector intersector;
	intersector.init(V, F);

	S.resize(GV.rows());

	for (int i = 0; i < GV.rows(); i++)
	{
		was_hit = intersector.intersectRay(view, GV.row(i)-view.transpose(), hit , 0.0, 1.0);
		if (was_hit)
		{
			// fill voxel grid cell
			S(i) = -1;
			//std::cout << "fill grid cell" << std::endl;

		}
		else
		{
			S(i) = 1;
		}

	}
}

void compute_shadow_volume(
	const Eigen::MatrixXf &V,
	const Eigen::MatrixXi &F,
	const Eigen::MatrixXf &GV,
	const Eigen::MatrixXf &views,
	Eigen::VectorXf &S)
{
	igl::embree::EmbreeIntersector intersector;
	intersector.init(V, F);

	S.resize(GV.rows());
	S = Eigen::VectorXf::Ones(GV.rows());

	for (int v = 0; v < views.rows(); v++)
	{
		const auto &view = views.row(v);
		// intersect ray with boxed mesh
		bool was_hit;
		igl::Hit hit;

		//Eigen::VectorXf S_i;
		//S_i.resize(GV.rows());

		for (int i = 0; i < GV.rows(); i++)
		{
			was_hit = intersector.intersectRay(view, GV.row(i) - view, hit, 0.0, 1.0);
			if (was_hit)
			{
				// fill voxel grid cell
				S(i) += 1;

			}
			else
			{
				S(i) = 0;
			}

		}
		//S = S.array() * S_i.array();
	}

	S = S / views.rows();
	S = (S.array() == 0).select(S.array() + 1, -S); //.select(1, -1);
}