/* 
 *  MODEL: Lattice model (diamond, foam)
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
#include "func.h"
#include "var.h"
#include "tool.h"
#include "domain.h"
#include "gui.h"
#include "templates.cc"

#define DIAMOND

extern void getvec(int loc, DCell *c, double *x, double *v);

namespace ModLattice
{
	const	int
		Nbonds=4,
		Nbonds_1=Nbonds-1;
	int
		ivolc,
		imasn,
		iveln,
		idmmn,
		infor,// node-forces (forces at the nodes)
		ixold;
	const	double
		sqrt3=sqrt(3.0),
		pi=4.0*atan(1),
		piover2=pi/2.0,
		node_mass=1.0,
		bond_length=5.0,
		small_length=0.3*bond_length,
		bond_angle=atan(sqrt(2.0)),
		alpha=2*bond_angle-piover2,
		sinalp=sin(alpha),
		cosalp=cos(alpha),
		bond0=bond_length*sinalp,
		bond1=bond_length*cosalp,
		Cstrain=1.0e0;//normal stress constant for boundary elasticity
	enum VariablesModLattice
	{
		cvolume,//cell-volume
		mass,
		velocity,
		dmomentum,
		nforce,
		xold,//old coordinates of vertexes
		maxVarModLattice
	};
	void	allocBonds
	(	Bond	*&root
	)
	{	//Allocate bonds
		Bond	*bond;
		root=new Bond;
		bond=root;
		bond->node=NULL;
		for (int ib=1; ib<Nbonds; ib++)
		{	bond->next=new Bond;
			bond=bond->next;
			bond->node=NULL;
///			bonds[ib]=bond;
		}
		bond->next=root;
	}
	void	createRootBonds
	(	double	x0[],//coordinates of 0-bond neighbor
		double	e1[],//vector in the plane of bond 1
		DNode	*node
	)
	{
		double	d,e[DIM],e0[DIM],*x=node->x;
		Bond	*&root=node->bond,*bond;///*bonds[Nbonds];
		allocBonds(root);
///		node->bond=bonds[0];
		//direction of 0-bond:
		for(int i=0;i<DIM;i++) {
			e0[i]=x0[i]-x[i];
		}
		d=LENGTH(e0);
#ifdef DEBUG
		if(d<SMALL)ERROR("Node-separation is too small to create a bond\n");
#endif
		for(int i=0;i<DIM;i++)e0[i]/=d;
		//Direction normal to e0:
		d=SCLP(e0,e1);
		for(int i=0;i<DIM;i++) {
			e1[i]-=d*e0[i];
		}
		d=LENGTH(e1);
#ifdef DEBUG
		if(d<SMALL)ERROR("Bond-angle is too small to create a bond\n");
#endif
		for(int i=0;i<DIM;i++)e1[i]/=d;
#ifdef	DIAMOND
		//Rotate e1 60deg around e0 
		VECP(e,e0,e1);//normal to e0,e1
		//normalize e:
		d=LENGTH(e);for(int i=0;i<DIM;i++)e[i]/=d;
		for(int i=0;i<DIM;i++) {
			e1[i]=0.5*(e1[i]+sqrt3*e[i]);
		}
#endif
		for(int i=0;i<DIM;i++)
			root->x[i]=x0[i];
		bond=root->next;
		while(1)
		{	double	d,*x=bond->x;
			for(int i=0;i<DIM;i++) {
				x[i]=-e0[i]*bond0+e1[i]*bond1;
			}
			bond=bond->next;
			if(bond==root)break;
			//Rotate e1 120deg around e0
			VECP(e,e0,e1);//normal to e0,e1
			//normalize e:
			d=LENGTH(e);for(int i=0;i<DIM;i++)e[i]/=d;
			for(int i=0;i<DIM;i++) {
				e1[i]=0.5*(-e1[i]+sqrt3*e[i]);
			}
		}
	}
	void	createBonds
	(
		double	e1[],//vector in the plane of bond 1
		DNode	*oldnode,
		DNode	*node
	)
	{
		double	d,e[DIM],e0[DIM],
			*x0=oldnode->x,//coordinates of 0-bond neighbor
			*x=node->x;
		Bond	*&root=node->bond,*bond;///*bonds[Nbonds];
		allocBonds(root);
///		node->bond=bonds[0];
		//direction of 0-bond:
		for(int i=0;i<DIM;i++)
			e0[i]=x0[i]-x[i];
		d=LENGTH(e0);
#ifdef DEBUG
		if(d<SMALL)ERROR("Node-separation is too small to create a bond\n");
#endif
		for(int i=0;i<DIM;i++)e0[i]/=d;
		//Direction normal to e0:
		d=SCLP(e0,e1);
		for(int i=0;i<DIM;i++) {
			e1[i]-=d*e0[i];
		}
		d=LENGTH(e1);
#ifdef DEBUG
		if(d<SMALL)ERROR("Bond-angle is too small to create a bond\n");
#endif
		for(int i=0;i<DIM;i++)e1[i]/=d;
#ifdef	DIAMOND
		//Rotate e1 60deg around e0
		VECP(e,e0,e1);//normal to e0,e1
		//normalize e:
		d=LENGTH(e);for(int i=0;i<DIM;i++)e[i]/=d;
		for(int i=0;i<DIM;i++) {
			e1[i]=0.5*(e1[i]+sqrt3*e[i]);
		}
#endif
		root->node=oldnode;
		bond=root->next;
		while(1)
		{	double	d,*x=bond->x;
			for(int i=0;i<DIM;i++)
				x[i]=-e0[i]*bond0+e1[i]*bond1;
			bond=bond->next;
			if(bond==root)break;
			//Rotate e1 120deg around e0
			VECP(e,e0,e1);//normal to e0,e1
			//normalize e:
			d=LENGTH(e);for(int i=0;i<DIM;i++)e[i]/=d;
			for(int i=0;i<DIM;i++) {
				e1[i]=0.5*(-e1[i]+sqrt3*e[i]);
			}
		}
	}
	void	relax
	(	int	ixold,
		int	ivolc,
		int	infor,
		double	relaxation
	)
	{

	}
};

using namespace ModLattice;

int	getNvarModLattice()
{
	return maxVarModLattice;
}
void	defVarModLattice
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
		case velocity:
			strcpy(name,"Velocity");
			type=nodes;
			dim=DIM;	
			break;
		case dmomentum:
			strcpy(name,"AddedMomentum");
			type=nodes;
			dim=DIM;
			break;
		case xold:
			strcpy(name,"OldNodeCoordinates");
			type=nodes;
			dim=DIM;
			break;
		case nforce:
			strcpy(name,"BoundaryForce");
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
void	initVarModLattice
(	int	ivar,
	char	*name,
	int	type,
	int	loc,
	int	dim,
	double	*val,
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
void	initModLattice(Domain *dom)
{	int
		mnvar=dom->getVarBufSize(nodes);
	Tool	*&tool=dom->tool;
	Variable	*variable=dom->variable;
	ivolc=variable[cvolume].loc;
	imasn=variable[mass].loc;
	iveln=variable[velocity].loc;
	idmmn=variable[dmomentum].loc;
	ixold=variable[xold].loc;
	infor=variable[nforce].loc;
	DNode	*node_root=dom->dnode_root,
		*node=node_root;
//	DCell	*cell_root=dom->dcell_root,
//		*cell=cell_root;
//	dom->indexNeibNodes();
//	dom->setBoundaryCellFlags();
//	dom->connectBoundaryCellsToFaces(dom->bface_root);
//	dom->createBoundaryNodeList(dom->bnode_root);
//	dom->createBoundaryCellList(dom->bcell_root);
//	dom->segm=new Segment;//needle
	if(node_root==NULL) return;
	do//Initialize nodal variables
	{	double
			x0[DIM],//coordinates of 0-bond
			e1[DIM],//vector in the plane of of bond 1
			*var=node->var,
			*mas=var+imasn,
			*vel=var+iveln,
			*dmm=var+idmmn,
			*bfr=var+infor;
		*mas=node_mass;
		for (int i=0; i<DIM; i++)
		{	int	ii=DIM*i;
			vel[i]=0.0;
			dmm[i]=0.0;
			bfr[i]=0.0;
			x0[i]=0.0;
		}
		x0[0]=-bond_length;
		//Vector in the plane of bond 1:
		e1[0]=0.0;e1[1]=1.0;e1[2]=0.0;
		createRootBonds(x0,e1,node);
		node=node->next;
	}	while(node!=node_root);
	{	Bond
			*oldbond=node->bond,
			*root=oldbond->next,
			*bond=root;
		double	*x=node->x;
		do
		{	DNode	*newnode=new DNode(mnvar,bond->x);
			double	*xold,e[DIM];
			dom->insert(node,newnode);
			newnode->state.boundary=1;	
			bond->node=newnode;
			//Create bonds
			if(oldbond->node==NULL)
				xold=oldbond->x;
			else
				xold=oldbond->node->x;
			for(int i=0;i<DIM;i++)e[i]=xold[i]-x[i];
			createBonds(e,node,newnode);
			oldbond=bond;
			bond=bond->next;
		}	while(bond!=root);
		node->state.boundary=0;
	}

//		if (tool!=NULL)
//		{
//			tool->mode.active=tool->mode.active?0:1;
//			if (option.verbose)printf("Tool active = %s\n",tool->mode.active?"YES":"NO");
//			tool->current_framenode=tool->first_framenode;
//			if(tool->first_framenode==NULL)ERROR("No tool framenodes specified");
//			tool->pos(tool->current_framenode->x);
//			tool->setradius(tool->current_framenode->radius);
//			tool->setrefinement(tool->current_framenode->refinement);
//			tool->setVolInd(nvolume);
//			tool->setXoldInd(xold);
//			Ntoolbnodes=tool->createLists(dom);
//			tool->setype(tool->current_framenode->type);
//			tool->setboundary();
//		}
}
void	stepModLattice(double dt, Domain *dom)
{	int
		mnvar=dom->getVarBufSize(nodes);
//		const	double	dtool=1.0;
//		Tool	*tool=dom->tool;
//		if 
//		(	tool!=NULL&&tool->first_framenode!=NULL
//			&&(int)((time.current-time.step)/dtool)!=(int)(time.current/dtool)
//			&&time.current>10.0
//		)
//		{//MOVE TOOL TO THE NEXT POSITION ALONG THE FRAME
//			Frame
//				*&first_framenode=tool->first_framenode,
//				*&last_framenode=tool->last_framenode,
//				*&current_framenode=tool->current_framenode;
//			tool->push();
//			if (current_framenode!=last_framenode)
//				current_framenode=current_framenode->next;
//			tool->setype(current_framenode->type);
//			tool->pos(current_framenode->x);
//			tool->setradius(current_framenode->radius);
//			tool->cleanSurface();
//			Ntoolbnodes=tool->createLists(dom);
//		}
	DNode	*node_root=dom->dnode_root,
		*node=node_root;
	if(node_root!=NULL)
	do
	{	double	*x=node->x;
		if(node->bond==NULL)
		{	Bond	*&bond_root=node->bond,*bond;
			allocBonds(bond_root);
			ERROR("Bonds not defined");///DDD
		}
		else
		{	if(node->state.boundary)
			{	int	inserted=0;
				Bond
					*oldbond=node->bond,
					*root=oldbond->next,
					*bond=root;
				do
				{	int	freebond=1;
					double	xbond[DIM],*z=bond->x;
					for(int i=0;i<DIM;i++)xbond[i]=x[i]+z[i];
					//check if the bond is close to another node
					DNode	*nd=node->next;
					do
					{	double	*y=nd->x,
							d=0.0;
#ifndef DEBUG
						if(nd->state.boundary)
						{
#endif
						for(int i=0;i<DIM;i++)
						{	double	r=y[i]-xbond[i];
							d+=r*r;
						}

						if(sqrt(d)<small_length)
						{
#ifdef DEBUG
							if(freebond==0){ERROR("Lattice conjestion");}
#endif
printf("Bond fusion\n");///DDD
							bond->node=nd;
							freebond=0;
#ifndef DEBUG
							break;
#endif
						}
#ifndef DEBUG
						}
#endif
						nd=nd->next;
					}	while(nd!=node);
					if(freebond)
					{//create new node
						DNode	*newnode;
						double	e[DIM],*xold;
						newnode=new DNode(mnvar,xbond);
						dom->insert(node,newnode);inserted++;
						newnode->state.boundary=1;	
						bond->node=newnode;
						//Create bonds
						if(oldbond->node==NULL)
							xold=oldbond->x;
						else
							xold=oldbond->node->x;
						for(int i=0;i<DIM;i++)e[i]=xold[i]-x[i];
						createBonds(e,node,newnode);
					}
					bond=bond->next;
				}	while(bond!=oldbond);
				node->state.boundary=0;
				for(int i=0;i<inserted;i++)node=node->next;
			}
		}
		node=node->next;
	}	while(node!=node_root);
	relax//Relax internal nodes
	(	ixold,
		ivolc,
		infor,
		0.5 //relaxation factor
	);
	runtime.current+=dt;
}
