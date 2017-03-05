
/* 
 *  MODEL: RFG
 *  Template of a generic model
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "main.h"
#include "io.h"
#include "vecalg.h"
#include "geom.h"
#include "var.h"
#include "tool.h"
#include "domain.h"
#include "templates.cc"

extern "C"
{
	void	genspec_(int *N);
	void	genvel_(double *t, double *x, double *tt, double *tl, double *uu, double *vel);
};

namespace ModRFG
{
	const double 
		init_mass=1.0;
	int	iveln;//nodal velocity memory offset
	enum VariablesModRFG
	{
		velocity=0,
		maxVarModRFG
	};
void	genvel
//Simplified genvel: only normal UU components are needed
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
)
{	double	turb_leng_scale[DIM],
		UU[2*DIM];
	turb_leng_scale[0]=length_scale_x;
	turb_leng_scale[1]=length_scale_y;
	turb_leng_scale[2]=length_scale_z;
	for (int i=0; i<2*DIM; i++)UU[i]=0.0;
	UU[0]=uu;
	UU[2]=vv;
	UU[5]=ww;
	genvel_
	(
		&time,
		x,
		&time_scale,
		turb_leng_scale,
		UU,						
		vel
	);
}
void rfg(int iveln, Domain *dom)
{
	Variable	*variable=dom->variable;
	DNode	*node_root=dom->dnode_root,
		*node=node_root;
	DCell	*cell_root=dom->dcell_root,
		*cell=cell_root;
	if(node_root==NULL) return;
	do//Assign RFG field to the nodal velocities
	{	double
			t=runtime.current,
			turb_time_scale=.1,
			turb_leng_scale[DIM],
			UU[2*DIM],
			*x=node->x,
			*var=node->var,
			*vel=var+iveln;
		for (int i=0; i<DIM; i++)
		{	int	ii=DIM*i;
			vel[i]=0.0;
			turb_leng_scale[i]=1.0;
		}
		turb_leng_scale[0]=1.0;
		turb_leng_scale[1]=1.0;
		turb_leng_scale[2]=1.0;
		for (int i=0; i<2*DIM; i++)UU[i]=0.0;
		UU[0]=UU[2]=UU[5]=1.0;
		genvel_
		(
			&t,
			x,
			&turb_time_scale,
			turb_leng_scale,
			UU,						
			vel
		);
		node=node->next;
	}	while(node!=node_root);
}
};

using namespace ModRFG;

int	getNvarModRFG()
{
	return maxVarModRFG;
}
void	defVarModRFG
(
	int	ivar,
	char	*name,
	int	&type,
	int	&dim,
	Domain	*dom
)
{
	int idomain=dom-domain_root;
	*name='\0';
	switch(ivar)
	{
		case velocity:
			strcpy(name,"Velocity");
			type=nodes;
			dim=DIM;	
			break;
		default:
			fprintf
			(	stderr,
				"Specification of variable %d domain %d (%s) is incomplete\n",
				ivar+1,idomain+1,dom->name
			);
			return;
	}
}
void	initVarModRFG
(	int	ivar,
	char	*name,
	int	type,
	int	loc,
	int	dim,
	double	*val,
	Domain	*dom
)
{	iveln=dom->variable[velocity].loc;
	switch(type)
	{	case nodes:
			setZero(loc,dim,dom->dnode_root);
			break;
		case cells:
			setZero(loc,dim,dom->dcell_root);
			break;
		case boundary_faces:
			setZero(loc,dim,dom->dcell_root);
			break;
		case points:
		{	int	mp=dom->getMaxNoPoints();
			for (int i=0; i<mp*dim; i++)
				val[i]=0.0;
		}
		break;
		default:
			fprintf(stderr,"Can't initialize elements of type %d\n",type);
			exit(1);
	}
}
void	initModRFG(Domain *dom)
{
	int	ns=200;//RFG spectrum size
	printf("Initializing RFG with %d harmonics\n",ns);
	genspec_(&ns);//Initilizing RFG spectrum
	rfg(iveln,dom);
	//See mod_softbody.cc
	//Initialize points
	dom->initp();
	//Create particles
	if(dom->getMaxNoPoints()-dom->getNoPoints()<1) 
		ERROR("Not enough particles\n");
}
void	stepModRFG(double dt, Domain *dom)
{	const	double	dtool=1.0;
	Tool	*tool=dom->tool;
	static int i=0;
	if(option.verbose)printf("Iteration no.=%d\n",i++);
	rfg(iveln,dom);
	runtime.current+=dt;
}
