#include <stdlib.h>
#include <stdio.h>
//#include <iostream.h>
#include <string.h>
#include <math.h>
#include "main.h"
#include "io.h"
#include "vecalg.h"
#include "geom.h"
#include "func.h"
#include "var.h"
//#include "particles.h"
#include "tool.h"
#include "domain.h"
//#include "solver.h"
#include "force.h"

namespace Force
{
	//const int	maxbonds=6;//max number of bonds for each particle
	int	maxbranches=2, //10;//max no. of branches in the tree
		Nframenodes,Nparticles;
	const double	large_mass=10.,small_mass=1.0,
		rbond=0.5,
		average_mass=0.5*(large_mass+small_mass),
		force_range_scale=10.0,
		frame_force_strength=10.0,
		attractive_force_strength=20,
		repulsive_force_radius=1.0,
		repulsive_force_strength=5.0,
		relaxation_time=0.1;
	double	FR(double r)
	{//Short range repulsive force
		double	r0,x,force_strength,force;
		force=0.0;
		{	r0=force_range_scale*repulsive_force_radius;
			r*=force_range_scale;
			x=r-r0;
			force_strength=repulsive_force_strength;
		}
		if (x<0.0) force=1.0;
		else
			force=exp(-fabs(x));
//			force=-2.0*f0*x*exp(-x*x)+(x<-0.5?f0:x>0.5?-f0*exp(-x):0.0);///(r*r));
		return force_strength*force;
	}
	double	FB(double r)
	{//Bond force
		double	r0,x,force_strength,force;
		r0=force_range_scale*rbond;
		r*=force_range_scale;
		x=r-r0;
		force_strength=50;
		if (x<-0.5) force=1.0;
		else 
		if (x>0.5) force=-1.0;
		else
			force=-x;
		return force_strength*force;
	}
	double	F(double r, double	m)
	{
		double	r0,x,force_strength,ftail,force;
		force=0.0;
		if (m>average_mass)
		{	r0=10.0;
			r*=force_range_scale;
			x=r-r0;
			force_strength=frame_force_strength;ftail=1.0;///fabs(x);
		}
		else
		{	r0=force_range_scale*rbond;
			r*=force_range_scale;
			x=r-r0;
			force_strength=attractive_force_strength;
			ftail=exp(-fabs(x-1.0));
//		force=-f0*2.0*x*exp(-x*x)+(x<-.5?f0:0.0);
		}
		if (x<-.2) force=1.0;
		else 
		if (x>.2) force=-1.0*ftail;
		else
			force=-5*x;
//			force=-2.0*f0*x*exp(-x*x)+(x<-0.5?f0:x>0.5?-f0*exp(-x):0.0);///(r*r));
		return force_strength*force;
	}
}
