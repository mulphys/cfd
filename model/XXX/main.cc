
/* 
 *  MODEL: XXX
 *  Template of a generic model
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "../main.h"
#include "../io.h"
#include "../vecalg.h"
#include "../geom.h"
#include "../var.h"
#include "../tool.h"
#include "../domain.h"
#include "../gui.h"
#include "../templates.cc"

namespace ModXXX
{
	const REAL 
		init_mass=1.0,
		init_pres=1.0;
	enum VariablesModXXX
	{
		nvolume=0,
		cvolume,
		mass,
		addmass,
		velocity,
		velgradc,
		dmomentum,
		pressure,
		dpressure,
		pmass,//particle mass
		pvelocity,//particle velocity
		maxVarModXXX
	};
	int //Indexes to some variables
		ivoln,//Node-centered volume
		ivolc,//Cell-centered volume
		imasn,//Node-centered mass
		iadmn,//Mass source term
		iveln;//Node-centered velocity
};

using namespace ModXXX;


int	getNvarModXXX()
{
	return maxVarModXXX;
}
void	defVarModXXX
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
		case nvolume:
			strcpy(name,"NodeVolume");
			type=nodes;
			dim=1;
			break;
		case cvolume:
			strcpy(name,"CellVolume");
			type=cells;
			dim=1;
			break;
		case mass:
			strcpy(name,"Mass");
			type=nodes;
			dim=1;
			break;
		case addmass:
			strcpy(name,"AddedMass");
			type=nodes;
			dim=1;
			break;
		case velocity:
			strcpy(name,"Velocity");
			type=nodes;
			dim=DIM;	
			break;
		case velgradc:
			strcpy(name,"CellVelocityGradient");
			type=cells;
			dim=DIM*DIM;
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
void	initVarModXXX
(	int	ivar,
	char	*name,
	int	type,
	int	loc,
	int	dim,
	REAL	*val,
	Domain	*dom
)
{	switch(type)
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
void	initModXXX(Domain *dom)
{
	Tool	*&tool=dom->tool;
	Variable	*variable=dom->variable;
	ivoln=variable[nvolume].loc;
	ivolc=variable[cvolume].loc;
	imasn=variable[mass].loc;
	iadmn=variable[addmass].loc;
	iveln=variable[velocity].loc;
	DNode	*node_root=dom->dnode_root,
		*node=node_root;
	DCell	*cell_root=dom->dcell_root,
		*cell=cell_root;
	if(node_root==NULL) return;
	do//Initialize nodal variables
	{	REAL
			*var=node->var,
			vol=var[ivoln],voli=1./vol,
			*mas=var+imasn,
			*dmn=var+iadmn,
			*vel=var+iveln;
		*mas=init_mass;
		for (int i=0; i<DIM; i++)
		{	int	ii=DIM*i;
			vel[i]=0.0;
		}
		node=node->next;
	}	while(node!=node_root);
	//Initialize face variables
	//See mod_softbody.cc
	//Initialize points
	dom->initp();
	//Create particles
	if(dom->getMaxNoPoints()-dom->getNoPoints()<1) 
		ERROR("Not enough particles\n");
	if (tool!=NULL)
	{
		tool->mode.active=tool->mode.active?0:1;
		if (option.verbose)printf("Tool active = %s\n",tool->mode.active?"YES":"NO");
		tool->current_framenode=tool->first_framenode;
		if(tool->first_framenode==NULL)ERROR("No tool framenodes specified");
		tool->pos(tool->current_framenode->x);
		tool->setradius(tool->current_framenode->radius);
		tool->setrefinement(tool->current_framenode->refinement);
		tool->setVolInd(nvolume);
		//tool->setXoldInd(xold);
		tool->setype(tool->current_framenode->type);
		tool->setboundary();
	}
}
void	stepModXXX(REAL dt, Domain *dom)
{	const	REAL	dtool=1.0;
	Tool	*tool=dom->tool;

	time.current+=dt;
}
