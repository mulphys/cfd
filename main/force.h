namespace Force
{
	//const int	maxbonds=6;//max number of bonds for each particle
	extern const double
		rbond,
		large_mass,small_mass,
		average_mass,
		force_range_scale,
		frame_force_strength,
		attractive_force_strength,
		repulsive_force_radius,
		repulsive_force_strength,
		relaxation_time;
	double	FR(double r);
	double	FB(double r);
	double	F(double r, double	m);
}
