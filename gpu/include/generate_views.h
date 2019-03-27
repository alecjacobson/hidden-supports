#ifndef GENERATE_VIEWS_H
#define GENERATE_VIEWS_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include <map>
#include <random>

void generate_views(
	const Eigen::Vector3f &bottom_left,
	const Eigen::Vector3f &top_right,
	const float num_views,
	Eigen::MatrixXf &views);

void generate_views(
	const Eigen::Vector3f &bottom_left,
	const Eigen::Vector3f &top_right,
	const float num_views,
	Eigen::MatrixXf &views)
{
	double mean_x = (bottom_left(0) + top_right(0)) / 2;
	double mean_y = (bottom_left(1) + top_right(1)) / 2;

	double std_dev_x = std::abs(bottom_left(0) - top_right(0)) / 6;
	double std_dev_y = std::abs(bottom_left(1) - top_right(1)) / 6;

	std::default_random_engine generator;
	std::normal_distribution<double> distribution_x(mean_x, std_dev_x);
	std::normal_distribution<double> distribution_y(mean_y, std_dev_y);

	views.resize(num_views, 3);
	int num_generated = 0;
	while (num_generated < num_views)
	{
		double x = distribution_x(generator);
		double y = distribution_y(generator);
		if (x >= bottom_left(0) && x < top_right(0))
		{
			if (y >= bottom_left(1) && y < top_right(1))
			{
				views.row(num_generated) = Eigen::Vector3f(x, y, -1);
				num_generated++;
			}
		}
	
	}

}

#endif