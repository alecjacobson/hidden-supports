#ifndef GENERATE_VIEWS_H
#define GENERATE_VIEWS_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include <map>
#include <random>
#include <igl/sortrows.h>

void generate_views_on_plane(
	const Eigen::Vector3f &bottom_left,
	const Eigen::Vector3f &top_right,
	const float num_views,
	const std::vector<float> z_vals,
	Eigen::MatrixXf &views);

void generate_views_on_plane(
	const Eigen::Vector3f &bottom_left,
	const Eigen::Vector3f &top_right,
	const float num_views,
	const std::vector<float> z_vals,
	Eigen::MatrixXf &views)
{
	double mean_x = (bottom_left(0) + top_right(0)) / 2;
	double mean_y = (bottom_left(1) + top_right(1)) / 2;

	double std_dev_x = std::abs(bottom_left(0) - top_right(0)) / 6;
	double std_dev_y = std::abs(bottom_left(1) - top_right(1)) / 6;

	std::default_random_engine generator;
	std::normal_distribution<double> distribution_x(mean_x, std_dev_x);
	std::normal_distribution<double> distribution_y(mean_y, std_dev_y);

	views.resize(num_views,3);

	int total_generated = 0;

	while (total_generated < (int)(num_views))// / z_vals.size()))
	{
		int num_generated = 0;

		while(num_generated < z_vals.size())
		{
			double x = distribution_x(generator);
			double y = distribution_y(generator);
			if (x >= bottom_left(0) && x < top_right(0))
			{
				if (y >= bottom_left(1) && y < top_right(1))
				{
					views.row(total_generated) = Eigen::Vector3f(x, y, z_vals[num_generated]);
					total_generated++;
					num_generated++;
				}
			}
		
		}
	}
	Eigen::MatrixXf tmp = views;
	Eigen::MatrixXi index;
	igl::sortrows(tmp, true, views, index);
}

#endif
