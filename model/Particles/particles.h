namespace ModParticles
{
///	void	InjectParticle(Domain *dom);
///	void	AdvanceParticles(double dt, Domain *dom);
	void	MoveParticle
	(
	 	Point	*origin,
		int	ip,//number of the particle to move
		double	dt,//time step
		double	pmass,//particle mass
		double	*fvel,//flow velocity at particle location
		double	*x,//particle coordinates
		double	*v //particle velocity
	);
	void	getBoundaryCells
	(	int	boundary,
		DCell	*cell_root,
		DCellList	*&boundary_cell_root
///		Domain	*dom
	);
//?	void	getBoundaryCellsStatic
//?	(	int	boundary,
//?		int	nc,
//?		Cell	*cell_root,
//?		CellList	*&boundary_cell_root
//?	);
};//END NAMESPACE 
