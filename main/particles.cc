
void	MoveParticle
(
 	Point	*origin,
	int	ip,//number of the particle to move
	double	dt,//time step
	double	pmass,//particle mass
	double	*fvel,//flow velocity at particle location
	double	*x,//particle coordinates
	double	*v //particle velocity
)
{
	static const double Cdrag=.1;
	int	imp=DIM*ip;
	Point	*p=origin+ip;
	for (int i=0; i<DIM; i++)
	{	double
			drag=Cdrag*(fvel[i]-v[i]),
			vnew=v[i]+drag/pmass*dt;
		x[i]+=0.5*(vnew+v[i])*dt;
		v[i]=vnew;
	}
}
