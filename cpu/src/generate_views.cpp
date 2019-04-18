#include "generate_views.h"
#include <iostream>

void generate_views(
	const Eigen::Vector3f &top_left,
	const Eigen::Vector3f &bottom_right,
	const float num_views,
	const std::vector<float> z_vals,
	Eigen::MatrixXf &views)
{
	double mean_x = (top_left(0) + bottom_right(0)) / 2;
	double mean_y = (top_left(1) + bottom_right(1)) / 2;


	double std_dev_x = std::abs(top_left(0) - bottom_right(0)) / 2;
	double std_dev_y = std::abs(top_left(1) - bottom_right(1)) / 2;

	std::default_random_engine generator;
	std::normal_distribution<double> distribution_x(mean_x, std_dev_x);
	std::normal_distribution<double> distribution_y(mean_y, std_dev_y);

	views.resize(num_views,3);

	int total_generated = 0;

	while (total_generated < (int)(num_views))// / z_vals.size()))
	{
		std::cout << total_generated << " of " << num_views << std::endl;
		int num_generated = 0;

		for(int i = 0; i < z_vals.size(); i++)
		{
			double x = distribution_x(generator);
			double y = distribution_y(generator);
			std::cout << x << ", " << y << std::endl;
			if (x >= top_left(0) && x < bottom_right(0))
			{
				if (y <= top_left(1) && y > bottom_right(1))
				{
					views.row(total_generated) = Eigen::Vector3f(x, y, z_vals[i]);
					
					total_generated++;
					num_generated++;
				}
			}
		
		}
	}


}
