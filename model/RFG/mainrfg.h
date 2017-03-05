namespace ModRFG
{
	void rfg
	(	int ivel, //memory offset to velocity variables
		Domain *dom
	);
	void	genvel
	(	double	time,
		double	*x,
		double	time_scale,
		double	length_scale_x,
		double	length_scale_y,
		double	length_scale_z,
		double	uu,
		double	vv,
		double	ww,
		double	*vel
	);
};
