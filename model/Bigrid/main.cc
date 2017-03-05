#include <stdlib.h>
#include <stdio.h>
//#include <iostream.h>
#include <string.h>
#include <math.h>
#include "../main.h"
#include "../io.h"
#include "../vecalg.h"
#include "../geom.h"
#include "../func.h"
#include "../var.h"
//#include "../particles.h"
#include "../tool.h"
#include "../domain.h"
#include "../force.h"
#include "../templates.cc"

#define	BRANCHING
//#define	GLOBAL_LIST_UPDATE

#define	MAXPRESITER	8	//number of pressure sub-iterations

#ifdef BRANCHING
	void	branch
	(//Symetric branching
		int	level,
		REAL	size,
		REAL	*x0,
		REAL	*e0,
		Tool	*tool
	);
	void	brancha
	(//Asymetric linear branching
		int	level,
		REAL	size,
		REAL	*x0,
		REAL	*e0,
		Tool	*tool
	);
	void	branch
	(//Asymmetric parabolic branching
		int	level,
		REAL	size,
		REAL	r,  //position on the parabolic branch
		REAL	*x0,
		REAL	*e0,
		Tool	*tool
	);
#endif

void	Advance
(	//Advance all variables by niter iterations 
	REAL	new_time,
	Domain	*dom
);

namespace ModBiGrid
{
	enum Variables
	{
		volume=0,
		xold,
		maxvar
	};
#ifdef BRANCHING
	int	maxbranches=2; //10;//max no. of branches in the tree
	const REAL	branch_size=3.0;
#endif
	int	tool_boundary_nodes;
}//END NAMESPACE MOD0

using namespace ModBiGrid;

int	getNvarModBiGrid()
{
	return maxvar;
}
void	defVarModBiGrid
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
void	initVarModBiGrid
(	int	ivar,
	char	*name,
	int	type,
	int	loc,
	int	dim,
	REAL	*val, //val[0:size*DIM^rank]
	Domain	*dom
)
{	switch(type)
	{
		case nodes:
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
//-{	int
//-		nrank=(int)pow(DIM,rank);
//-	if (option.debug)
//-	{	printf
//-		(	"Initializing variable %d: %s, rank=%d, totsize=%d x %d\n",
//-			ivar,name,rank,size,nrank
//-		);
//-	}
//-	switch(ivar)
//-	{
//-		case volume:
//-		case xold:
//-		if (dom->dnode_root!=NULL) 
//-		{	int
//-//-				nscl=dom->getNoScalars(),
//-//-				nvec=dom->getNoVectors(),
//-//-				nvar=nscl+nvec,
//-//-				mvar=nvec*DIM+nscl;
//-				mvar=dom->getVarBufSize(nodes);
//-			DNode
//-				*dnode_root=dom->dnode_root,
//-				*node=dnode_root;
//-			if (option.verbose)
//-			{	printf
//-				(	"Initializing %d*%d bytes for variable no. %d at each node ... ",
//-					mvar,sizeof(REAL),ivar
//-				);FLUSH;
//-			}
//-			do
//-			{	REAL	*var=node->var;
//-				for (int i=0; i<mvar; i++)
//-					var[i]=0.0;
//-				node=node->next;
//-			}	while(node!=dnode_root);
//-			if (option.verbose)printf("done\n");
//-		}
//-		break;
//-		default:
//-			for (int i=0; i<size*nrank; i++)
//-				val[i]=0.0;
//-			break;
//-	}
//-	if (option.debug)
//-	{	printf
//-		(	"Variable initialized\n"
//-		);FLUSH;
//-	}
//-}
///void	Tool::frameForce
///(	Point	*a,
///	REAL *xa,//Coordinates of a point 
///	REAL *X, //Global coordinates array
///	REAL *M, //Global mass-array
///	REAL *fframe //force of the frame
///)
///{//Copute the force of the frame
///	int	ia=a-origin; //index of point a
///	REAL
///		ma=M[ia],
///		dframe=LARGE;//distance to the frame
///	Point	*b=first_framenode;//particle 'b'
///	for (int i=0; i<DIM; i++)fframe[i]=0.0;
///	for (int ib=0; ib<Nframenodes; ib++,b=b->next)
///	{	if (a!=b)
///		{	int	jb=b-origin,jmb=DIM*jb;
///			REAL	d[DIM],dd,//distance between 'a' and 'b'
///				*xb=X+jmb,
///				mb=M[jb],//mass of particle b
///				force;
///			dd=0.0;
///			for (int k=0; k<DIM; k++)
///			{	REAL	r=xa[k]-xb[k];
///				dd+=r*r;
///				d[k]=r;
///			}
///			dd=sqrt(dd);
///			if (dd<SMALL)dd=SMALL;
///			if (dd<dframe)
///			{	dframe=dd;
///				force=ma*mb*F(dd,mb)/dd;
///				for (int i=0; i<DIM; i++)
///					fframe[i]=force*d[i];
///			}
///		}
///	}
///}
void	initModBiGrid(Domain *dom)
//Called once before Step
{	//Create a binary tree
	Tool	*&tool=dom->tool;
	DCell	*&dcell_root=dom->dcell_root;
	BFaceList	*&bface_root=dom->bface_root;
#ifdef BRANCHING
	int	nbranches=maxbranches;
	REAL	x[DIM],e[DIM];//initial position and direction
	//Select the starting position
	x[0]=x[1]=x[2]=0.0;
	//Select initial direction
	e[0]=e[1]=0.0;e[2]=1.0;
	if (tool!=NULL)
///		branch(nbranches,branch_size,x,e,tool);
		brancha(nbranches,branch_size,x,e,tool);
#endif
	if (tool!=NULL)
	{
		tool->mode.active=tool->mode.active?0:1;
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
///	}
///	dom->deleteList(bface_root);
///	dom->createCellFaceList(bface_root);
///	if (tool!=NULL)
///	{
		tool->setype(tool->current_framenode->type);
		tool->setboundary();
	}
}

// ONE ITERATION STEP

void	stepModBiGrid(REAL dt, Domain *dom)
//Advance all variables by one timestep 
{	Tool	*tool=dom->tool;
	DNode	*&dnode_root=dom->dnode_root;
		//REFINE THE GRID
		//Cycle the near boundary cells and split all cells 
		// where boundary or near-boundary nodes are too far apart
		//do//TOOL ADJUSTMENT CYLCE
//	while(tool_boundary_nodes>0)
		if (tool_boundary_nodes>0)
			tool_boundary_nodes=tool->step(dom);
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
	time.current+=time.step;
}
void	Advance
(	//Advance all variables by niter iterations 
	REAL	new_time,
	Domain	*dom
)
{	const	int	iz=1;
	DNode	*node;
	REAL dt=time.step,
		old_timestep=time.step0;
	if (option.debug)printf("Advancing to time %g\n",new_time);
}
#ifdef	BRANCHING
//SYMMETRIC LINEAR BRANCHING 
void	branch
(	int	level,
	REAL	stem_radius,
	REAL	*x0,
	REAL	*e0,
	Tool	*tool
)
{	int
		nstem,//number of points along the branch stem
		ntail; //points along the branch tail
	const	REAL
		refinement_factor=0.7,
		branching_factor=0.4,
		sqrt2=sqrt(2.),
		stem_length=6.0;
	REAL	r,r2,
		tail_length,
		branching,//to be determined
		segment_length=0.2*stem_radius,
		next_segment_length,next_radius,
		x[DIM],x1[DIM],e[DIM],e1[DIM];
	if (level--==0)return;
	next_radius=stem_radius/sqrt2;
	next_segment_length=segment_length/sqrt2;
	nstem=(int)ceil(stem_length/segment_length)+1;
	for (int i=0; i<DIM; i++)x[i]=x0[i];
	//Do nstem steps 
	for (int i=0; i<nstem; i++)
	{ //create particle 
		tool->addnode(x,(int)boundary,stem_radius,refinement_factor);
		for (int k=0; k<DIM; k++)
			x[k]+=segment_length*e0[k];
	}
	//Store the end of the stem
	for (int i=0; i<DIM; i++)x1[i]=x[i];
	//Select new direction
	do
	{	randvec(e);//random unit vector
		//Determine the new direction
		r2=0.0;
		for (int i=0; i<DIM; i++)
		{	REAL	a=(1.-branching_factor)*e0[i]+branching_factor*e[i];
			r2+=a*a;
			e[i]=a;
		}
		r=sqrt(r2);
		for (int i=0; i<DIM; i++) e[i]/=r;		
		//Determine the tail length
		branching=0.0;
		for (int i=0; i<DIM; i++) 
		{	REAL	b=e[i]-e0[i];
			branching+=b*b;
		}
		branching=sqrt(branching);
	} while(branching<0.05);
	tail_length=next_radius/(sqrt2*branching);
	ntail=(int)ceil(tail_length/next_segment_length)+1;
	//Do ntail steps 
	segment_length=tail_length/ntail;
	for (int i=0; i<ntail; i++)
	{	r=stem_radius+(REAL)(i+1)/(REAL)ntail*(next_radius-stem_radius);
		for (int k=0; k<DIM; k++)
			x[k]+=segment_length*e0[k];
		//create particle 
		tool->addnode(x,(int)boundary,r,refinement_factor);
	}
	if (level==0)
	{//open tip
		for (int k=0; k<DIM; k++)
			x[k]+=segment_length*e0[k];
 		tool->addnode(x,(int)outlet,r,refinement_factor);
	}
	//Create the first branch
	branch(level,next_radius,x1,e,tool);
	//Mirror-reflect e relative to e0:
	r=SCLP(e,e0);
	for (int i=0; i<DIM; i++)
		e[i]=2.*r*e0[i]-e[i];
	//Create second branch
	branch(level,next_radius,x1,e,tool);
}
//ASYMMETRIC LINEAR BRANCHING 
void	brancha
(	int	level,
	REAL	stem_radius,
	REAL	*x0,
	REAL	*e0,
	Tool	*tool
)
{	int
		nstem,//number of points along the branch stem
		ntail; //points along the branch tail
	const	REAL
		refinement_factor=0.7,///0.6,
		branching_factor=0.4,
		r0=0.55,
		r1=0.85,
		stem_length=9.0;
	REAL	r,r2,
		tail_length,
		branching,//to be determined
		segment_length=0.2*stem_radius,
		next_segment_length0,next_radius0,
		next_segment_length1,next_radius1,
		x[DIM],x1[DIM],e[DIM],e1[DIM];
	if (level--==0)return;
	next_radius0=r0*stem_radius;
	next_segment_length0=0.5*(1.0+r0)*segment_length;
	next_radius1=r1*stem_radius;
	next_segment_length1=0.5*(1.0+r1)*segment_length;
//next_segment_length0=next_segment_length1;//=segment_length;///DDD
	nstem=(int)ceil(stem_length/segment_length)+1;
	for (int i=0; i<DIM; i++)x[i]=x0[i];
	//Do nstem steps 
	for (int i=0; i<nstem; i++)
	{ //create particle 
		tool->addnode(x,(int)boundary,stem_radius,refinement_factor);
		for (int k=0; k<DIM; k++)
			x[k]+=segment_length*e0[k];
	}
	//Store the end of the stem
	for (int i=0; i<DIM; i++)x1[i]=x[i];
	//Select new direction
	do
	{	randvec(e);//random unit vector
		//Determine the new direction
		r2=0.0;
		for (int i=0; i<DIM; i++)
		{	REAL	a=(1.-branching_factor)*e0[i]+branching_factor*e[i];
			r2+=a*a;
			e[i]=a;
		}
		r=sqrt(r2);
		for (int i=0; i<DIM; i++) e[i]/=r;		
		//Determine the tail length
		branching=0.0;
		for (int i=0; i<DIM; i++) 
		{	REAL	b=e[i]-e0[i];
			branching+=b*b;
		}
		branching=sqrt(branching);
	} while(branching<0.4||branching>0.6);//limit the branching angle
	tail_length=r1*next_radius1/branching;
	ntail=(int)ceil(tail_length/next_segment_length1)+1;
	//Do ntail steps 
	segment_length=tail_length/ntail;
	for (int i=0; i<ntail; i++)
	{	r=stem_radius+(REAL)(i+1)/(REAL)ntail*(next_radius1-stem_radius);
		for (int k=0; k<DIM; k++)
			x[k]+=segment_length*e0[k];
		//create particle 
		tool->addnode(x,(int)boundary,r,refinement_factor);
	}
	if (level==0)
	{//open tip
		for (int k=0; k<DIM; k++)
			x[k]+=segment_length*e0[k];
 		tool->addnode(x,(int)presoutlet,r,refinement_factor);
	}
	//Create the first branch
	brancha(level,next_radius0,x1,e,tool);
	//Reduce the angle by half
	for (int k=0; k<DIM; k++)
		e[k]=0.5*(e[k]+e0[k]);
	//Mirror-reflect e relative to e0:
	r=SCLP(e,e0);
	for (int i=0; i<DIM; i++)
		e[i]=2.*r*e0[i]-e[i];
	//Create second branch
	brancha(level,next_radius1,x1,e,tool);
}
//ASYMMETRIC PARABOLIC BRANCHING 
void	branch
(	int	level,
	REAL	stem_radius,
	REAL	r,  //position on the parabolic branch
	REAL	*x0,//start coordinate
	REAL	*e0, //start direction of the small branch
	Tool	*tool
)
{	int
		nstem,//number of points along the branch stem
		ntail; //points along the branch tail
	static const	REAL
		rmin=0.1,rmax=2,
		branching_factor=0.4,
		sqrt2=sqrt(2.),
		stem_length=12.0;///12.0
	REAL
		tail_length,
		segment_length=0.2*stem_radius,
		next_segment_length,mass,next_radius,
		ee,e[DIM],//current unit direction vector
		e1[DIM],//second unit normal base-vector of the plane:
		        // SCLP(e0,e1)==0
		x[DIM],x1[DIM];
	if (level--==0)return;
	if (r<rmin)r=rmin;
	else
	if (r>rmax)r=rmax;
	//TODO: rewrite the rest below

//-	next_radius=stem_radius/sqrt2;
//-	next_segment_length=segment_length/sqrt2;
//-	nstem=(int)ceil(stem_length/segment_length)+1;
//-	for (int i=0; i<DIM; i++)x[i]=x0[i];
//-	//Do nstem steps 
//-	for (int i=0; i<nstem; i++)
//-	{	int	ip;
//-		for (int k=0; k<DIM; k++)
//-			x[k]+=segment_length*e0[k];
//-		//create particle 
//-		if ((ip=putp(x))<0) continue;
//-		setsclp(ip,tool_radius_index,stem_radius);
//-	}
//-	//Store the end of the stem
//-	for (int i=0; i<DIM; i++)x1[i]=x[i];
//-	//Select new direction
//-	randvec(e);//random unit vector
//-	//Determine the new direction
//-	ee=0.0;
//-	for (int i=0; i<DIM; i++)
//-	{	REAL	a=(1.-branching_factor)*e0[i]+branching_factor*e[i];
//-		ee+=a*a;
//-		e[i]=a;
//-	}
//-	ee=sqrt(ee);
//-	for (int i=0; i<DIM; i++) e[i]/=ee;		
//-	{	//Determine the tail length
//-		REAL	a=0.0;
//-		for (int i=0; i<DIM; i++) 
//-		{	REAL	b=e[i]-e0[i];
//-			a+=b*b;
//-		}
//-		a=sqrt(a);
//-		if (a<SMALL) return;//Braching can not be done
//-		tail_length=next_radius/(sqrt2*a);
//-		ntail=(int)ceil(tail_length/next_segment_length)+1;
//-	}
//-	//Do ntail steps 
//-	segment_length=tail_length/ntail;
//-	for (int i=0; i<ntail; i++)
//-	{	int	ip;
//-		REAL	size=stem_radius+(REAL)(i+1)/(REAL)ntail*(next_radius-stem_radius);
//-		for (int k=0; k<DIM; k++)
//-			x[k]+=segment_length*e0[k];
//-		//create particle 
//-		if ((ip=putp(x))<0) continue;
//-		setsclp(ip,tool_radius_index,size);
//-	}
//-	//Create the first branch
//-	branch(level,next_radius,x1,e);
//-	//Mirror-reflect e relative to e0:
//-	ee=SCLP(e,e0);
//-	for (int i=0; i<DIM; i++)
//-		e[i]=2.*ee*e0[i]-e[i];
//-	//Create second branch
//-	branch(level,next_radius,x1,e);
}
#endif
