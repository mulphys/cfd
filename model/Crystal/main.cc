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
#include "tool.h"
#include "domain.h"
#include "force.h"
#include "templates.cc"

/*
 * Generating new Tetra (T):
 * Pick up a triangular boundary face (ABC)
 * Compute the geometrical center point M
 * Determine the normal vector to ABC
 * Find point P lying on face-normal line starting at M
 *   at a distance PM=1/(2*sqrt(3))*Perimeter(ABC)
 * Find maxinum distance r0=MaxDistance(PA,PB,PC)
 * Consider a sphere S0=Sphere(P,r0)
 * Find all the boundary nodes B0, lying inside S0 and 
 *   on the same side of ABC as P
 * If B0=empty => Construct tetra T(ABCP); Stop;
 * Set h=PM
 * Compute d=MaxDistance(MA,MB,MC)
 * Select point Q on PM at a distance h1=0.5*(h-d^2/h) from M  (1)
 * Find max distance r1=MaxDistance(QA,QB,QC)                  (2)
 * Consider a sphere S1=Sphere(Q,r1)
 * Find subset B1 of B0, consisting of nodes inside S1
 * If B1!=empty
 * 	Call STMSV(B1)
 * 	Stop
 * Move the point P to avoid clustering in {B0,A,B,C,P}
 * Call STMSV({B1,P})
 * Stop.
 * 
 * Procedure STMSV(Nodes N): 
 * Select node P\in\,N such that Tetra PABC will have the biggest V/S ratio
 * among all other nodes in N, and such that not other node from N is 
 * inside PABC 
 *
 * NOTE: steps (1),(2) can be replaced by finding a 
 * circumsphere or ABCP
 */

/*

(H-h)^2=h^2+d^2 
H^2-2*H*h=d^2
h=0.5*(H-d^2/H)
   
*/


namespace ModCrystal
{
	enum Variables
	{
		volume=0,
		xold,
		maxvar
	};
	int	tool_boundary_nodes;
	void	Advance
	(	//Advance all variables by niter iterations 
		double	new_time,
		Domain	*dom
	);
}

using namespace ModCrystal;

void	ModCrystal::Advance
(	//Advance all variables by niter iterations 
	double	new_time,
	Domain	*dom
)
{	const	int	iz=1;
	DNode	*node;
	double dt=runtime.step,
		old_timestep=runtime.step0;
	if (option.debug)printf("Advancing to time %g\n",new_time);
}

int	getNvarModCrystal()
{
	return maxvar;
}
void	defVarModCrystal
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
		case volume:
			strcpy(name,"Volume");
			type=nodes;
			dim=1;
			break;
		case xold:
			strcpy(name,"Xold");
			type=nodes;
			dim=DIM;
			break;
		default:
			fprintf
			(	stderr,
				"Specification of variable %d domain %d (%s) is incomplete\n",
				ivar+1,idomain+1,name
			);
			return;
	}
}
void	initVarModCrystal
(	int	ivar,
	char	*name,
	int	type,
	int	loc,
	int	dim,
	double	*val, //val[0:size*DIM^rank]
	Domain	*dom
)
{	switch(type)
	{	case nodes:
			setZero(loc,dim,dom->dnode_root);
			break;
		case cells:
			setZero(loc,dim,dom->dcell_root);
			break;
		default:
			fprintf(stderr,"Can't initialize elements of type %d\n",type);
			exit(1);
	}
}
void	initModCrystal(Domain *dom)
//Called once before Step
{	//Create a binary tree
	Tool	*&tool=dom->tool;
	DCell	*&dcell_root=dom->dcell_root;
	BFaceList	*&bface_root=dom->bface_root;
	if (tool!=NULL)
	{	tool->mode.active=tool->mode.active?0:1;
		if (option.verbose)printf("Tool active = %s\n",tool->mode.active?"YES":"NO");
		tool->current_framenode=tool->first_framenode;
		if(tool->first_framenode==NULL)ERROR("No tool framenodes specified");
		tool->pos(tool->current_framenode->x);
		tool->setradius(tool->current_framenode->radius);
		tool->setrefinement(tool->current_framenode->refinement);
		tool->setVolInd(volume);
		tool->setXoldInd(xold);
		tool_boundary_nodes=tool->createLists(dom);
#ifndef GLOBAL_LIST_UPDATE
		tool_boundary_nodes=tool->updateLists(dom);
#endif
		tool->setype(tool->current_framenode->type);
		tool->setboundary();
	}
}

// ONE ITERATION STEP

void	stepModCrystal(double dt, Domain *dom)
//Advance all variables by one timestep 
{	Tool	*tool=dom->tool;
	DNode	*&dnode_root=dom->dnode_root;
	//REFINE THE GRID
	//Cycle the near boundary cells and split all cells 
	// where boundary or near-boundary nodes are too far apart
	//do//TOOL ADJUSTMENT CYLCE
//	while(tool_boundary_nodes>0)
	if (tool_boundary_nodes>0)
		tool_boundary_nodes=tool->growCrystal(dom);
	if (tool->first_framenode!=NULL&&tool_boundary_nodes==0)
	{//MOVE TOOL TO THE NEXT POSITION ALONG THE FRAME
		Frame
			*&first_framenode=tool->first_framenode,
			*&last_framenode=tool->last_framenode,
			*&current_framenode=tool->current_framenode;
		tool->relax(dom,0.3);
		if (current_framenode==last_framenode) current_framenode=first_framenode;
		else current_framenode=current_framenode->next;
		tool->setype(current_framenode->type);
		tool->pos(current_framenode->x);
		tool->setradius(current_framenode->radius);
		tool->setrefinement(current_framenode->refinement);
		tool->cleanSurface();
#ifdef GLOBAL_LIST_UPDATE
		tool_boundary_nodes=tool->createLists(dom);
#else
		tool_boundary_nodes=tool->updateLists(dom);
#endif
		tool->setboundary();
	}
	runtime.current+=runtime.step;
}

