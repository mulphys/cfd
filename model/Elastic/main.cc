/*
 *  MODEL: Viscous flow with particles and elastic walls
 *  Author: andrei.v.smirnov@gmail.com
 *  http://smirnov.mae.wvu.edu
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
#include "gasflow.h"
#include "templates.cc"

extern void getvec(int loc, DCell *c, double *x, double *v);

namespace ModElastic
{
	enum VariablesModElastic
	{
		nvolume=0,
		cvolume,
		xold,//old coordinates of vertexes
		eforce,//boundary forces
		teth,//cell heights
		trih,//boundary cell heights 
		trib,//boundary cell lateral positions
		//
		//               o--------------            
		//             _/|\_          ^
		//           _/  |  \_        |
		//         _/    |    \_      | height
		//       _/      |      \_    |
		//      /        |        \   v
		//     o---------+---------o----
		//     |<------->|
		//        base
		//
		maxVarModElastic
	};
	int
		Nparticles,
		Ntoolbnodes,
		ivoln,
		ivolc,
		iheig,
		itrib,
		ixold,
		ibfor;
	double
		Cshear,///shear stress constant for boundary elasticity
		Cstrain,//normal stress constant for boundary elasticity
		// Both Cshear, Cstrain are defined in model init-functions.
		Rgas=1.0,//redefined in Model.init functions
		Tamb=1.0,//Gas-constant of air [work/(mass*temp)
		Pamb=1.0e-0,//Ambient pressure [mass/(length*time^2)]
		Pin=1.2*Pamb,//Inlet pressure (unstable for bif1 at Pin>=1.2)
		Tin=Tamb,
		Cwall=0.5,//0.75;//Wall-friction coefficient [mass/(length^2*time)]
		visc=1.50e0;///-8//Viscosity, [length^2/time]
void	trihi
(
	double	*edge,//edges of the triangle [Nfv*DIM]
	double	*trib,//[Nfv]
	double	*heig,//[Nfv]
	double	*H //h[Nfv*DIM]
)
{	for (int i0=0; i0<Nfv; i0++)
	{	int
			j0=DIM*i0,
			i1=(i0+1)%Nfv,j1=DIM*i1;
		double	d,
			e[DIM],//unit vector in the direction of b
			*h=H+j0,
			*e0=edge+j0,
			*e1=edge+j1,
			ei=1./(LENGTH(e1));
		for(int i=0;i<DIM;i++)e[i]=ei*e1[i];
		d=-(SCLP(e0,e));
		trib[i0]=d*ei;
		for(int i=0;i<DIM;i++)
			h[i]=e0[i]+d*e[i];
		heig[i0]=LENGTH(h);
	}
}
void	init(Domain *dom)
{
	Tool	*&tool=dom->tool;
	Variable	*variable=dom->variable;
	ivoln=variable[nvolume].loc;
	ivolc=variable[cvolume].loc;
	iheig=variable[trih].loc;
	itrib=variable[trib].loc;
	ixold=variable[xold].loc;
	ibfor=variable[eforce].loc;
	DNode	*node_root=dom->dnode_root,
		*node=node_root;
	DCell	*cell_root=dom->dcell_root,
		*cell=cell_root;
	ModGasFlow::Rgas=ModElastic::Rgas;
	ModGasFlow::Tamb=ModElastic::Tamb;//Gas-constant of air [work/(mass*temp)
	ModGasFlow::visc=ModElastic::visc;///-8//Viscosity, [length^2/time]
	ModGasFlow::Pamb=ModElastic::Pamb;//Ambient pressure [mass/(length*time^2)]
	ModGasFlow::Pin =ModElastic::Pin;//Inlet pressure
	ModGasFlow::Tin =ModElastic::Tin;
	ModGasFlow::Cwall=ModElastic::Cwall;//Wall-friction coefficient [mass/(length^2*time)]
	dom->setCellVolumes(cvolume);
	dom->setNodeVolumes(nvolume,cvolume);
	dom->indexNeibNodes();
	dom->setBoundaryCellFlags();
	dom->setBoundaryVertexes();
	dom->connectBoundaryCellsToFaces(dom->bface_root);
	dom->createBoundaryNodeList(dom->bnode_root);
	dom->createBoundaryCellList(dom->bcell_root);
	//Relax the grid to put the nodes to the centroids: IMPORTANT!
	for(int i=0;i<100;i++)
	dom->relax//Relax internal nodes
	(	ixold,
	 	ivoln,
 		ivolc,
		0.0 //relaxation factor
	);
	//Initialize face variables
	if(dom->bface_root!=NULL)
	{	BFaceList
			*root=dom->bface_root,
			*face=root;
		do
		{	int	ivert=face->iface;
			double *y,
			*varf=face->var,
			*heig=varf+iheig,
			*base=varf+itrib,
			E[Nfv*DIM],///,d[Nfv],di[Nfv];//edge vectors lengths and inverses
			H[Nfv*DIM];
			DCell *cell=face->cell;
			DNode	**vert=cell->vert;
			y=vert[(ivert+Nv1)%Nv]->x;
			for(int iv=0; iv<Nfv; iv++)
			{	int	imv=DIM*iv;
				double
					*x=vert[(ivert+iv+1)%Nv]->x,
					*e=E+imv;
				for(int i=0; i<DIM; i++)
					e[i]=x[i]-y[i];
				y=x;
			}
			trihi(E,base,heig,H);
			face=face->next;
		}	while(face!=root);
	}
}
void	MoveBoundary(double	dt, Domain *dom)
{
	Variable	*variable=dom->variable;
	//Set boundary forces to zero
	if(dom->bnode_root!=NULL)
	{	BNodeList	*root=dom->bnode_root,
			*bnode=root;
		do
		{	DNode	*node=bnode->node;
			double
				*var=node->var,
				*bfor=var+ibfor;
			for (int i=0; i<DIM; i++)
				bfor[i]=0.0;
			bnode=bnode->next;
		}	while(bnode!=root);
	}
	//Assemble forces at the boundary
	if(dom->bface_root!=NULL)
	{	BFaceList
			*root=dom->bface_root,
			*face=root;
		do
		{	if(face->type>=boundary)
		{
			int
				ivert=face->iface,
				jvert=(ivert+Nv1)%Nv;
			DCell *cell=face->cell;
			DNode	**verts=cell->vert;
			double	*y,
				//denf=Damb,//density accross the face
				prsf=Pamb,//pressure accross the face
				area=face->area,///areaoverNfv=area/Nfv,
				*varn=verts[ivert]->var,//node-variables
				prsn=Pin,//pressure inside the cell
				dp=prsn-prsf,//pressure drop accross the boundary
				dpf=dp*area/(double)Nfv,
				*norm=face->norm,
				*varf=face->var,//face-varialbes
				*heig0=varf+iheig,
				*base0=varf+itrib,
				heig[Nfv],base[Nfv],
				E[Nfv*DIM],//edge vectors lengths
				H[Nfv*DIM];//height vectors
			//Compute edges and perimeter
			y=verts[jvert]->x;
			for(int iv=0; iv<Nfv; iv++)
			{	int	imv=DIM*iv;
				double
					*x=verts[(ivert+iv+1)%Nv]->x,
					*e=E+imv;
				for(int i=0; i<DIM; i++)
					e[i]=x[i]-y[i];;
				y=x;
			}
			trihi(E,base,heig,H);//Triangle strain
			for(int iv=0; iv<Nfv; iv++)
			{	//Compute forces at the boundary vertexes
				int	imv=DIM*iv,
					iv1=(iv+1)%Nfv,imv1=DIM*iv1,
					jvert1=(ivert+iv+1)%Nv;
				DNode	*vert=verts[jvert];
				double//Elasticity forces
					*e1=E+imv1,le1=LENGTH(e1),e1i=1./le1,
					*var=vert->var,
					*bfor=var+ibfor,//boundary force
					*hh=H+imv,
					h=heig[iv],hi,
					h0=heig0[iv],h0i=1./h0,
					b=base[iv],
					b0=base0[iv],
					dh=h0-h,
					db=b0-b,
					doublearea=le1*h,
					strain,shear;
#ifdef DEBUG
				if(h<SMALL)ERROR("Height too small in MoveGrid\n");
#endif
				hi=1./h;
				strain=Cstrain*dh*h0i*hi*doublearea;
				shear=Cshear*db*h0i*doublearea;
				//Assemble forces
				for (int i=0; i<DIM; i++)
					bfor[i]+=
						+shear*e1[i]
						-strain*hh[i]
						+dpf*norm[i];
				jvert=jvert1;
			}
		}
			face=face->next;
		}	while(face!=root);
	}
	//Move boundary
	if(dom->bnode_root!=NULL)
	{	const double	onethird=1./3.,
			flowdrag=0.90;//drag=0.5,cdrag=1./(1+drag).
		BNodeList
			*root=dom->bnode_root,
			*bnode=root;
		do//Update boundary node variables
		{	DNode	*node=bnode->node;
			double	*x=node->x,
				*var=node->var,
				vol=var[ivoln],
				cellsize=pow(vol,onethird),
				bigstep=0.1*cellsize,
				*bfor=var+ibfor,
				force=LENGTH(bfor);
			if(node->state.insidetool|node->state.toolsurface|node->state.fixed)goto loop;
			//Viscous motion
			if(force*dt<=bigstep)
			for(int i=0;i<DIM;i++)
				x[i]+=bfor[i]*dt;
			else
			{	double	factor=bigstep/force;
				for(int i=0;i<DIM;i++)
					x[i]+=factor*bfor[i];//*dt;
			}
			loop: bnode=bnode->next;
		}	while(bnode!=root);
	}
	dom->updateBoundaryFaceList(dom->bface_root);
}
};

using namespace ModElastic;

int	getNvarModElastic()
{
	return maxVarModElastic;
}
void	defVarModElastic
(
	int	ivar,
	char	*name,
	int	&type,
	int	&dimension,
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
			dimension=1;
			break;
		case cvolume:
			strcpy(name,"CellVolume");
			type=cells;
			dimension=1;
			break;
		case xold:
			strcpy(name,"OldNodeCoordinates");
			type=nodes;
			dimension=DIM;
			break;
		case eforce:
			strcpy(name,"ElasticForce");
			type=nodes;
			dimension=DIM;
			break;
		case teth:
			strcpy(name,"CellHeight");
			type=cells;
			dimension=4; //array of four heights
			break;
		case trih:
			strcpy(name,"BoundaryCellHeight");
			type=boundary_faces;
			dimension=DIM;
			break;
		case trib:
			strcpy(name,"Base");
			type=boundary_faces;
			dimension=DIM;
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
void	initVarModElastic
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
			setZero(loc,dim,dom->bface_root);
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
void	initElastic(Domain *dom)
{
	Tool	*&tool=dom->tool;
	DNode	*node_root=dom->dnode_root,
		*node=node_root;
	DCell	*cell_root=dom->dcell_root,
		*cell=cell_root;
	ModElastic::init(dom);
	//Initialize points
///	dom->initp();
	if(node_root==NULL) return;
	do//Initialize nodal variables
	{	double
			*var=node->var,
			vol=var[ivoln],voli=1./vol,
			*bfr=var+ibfor;
		for (int i=0; i<DIM; i++)
		{	
			bfr[i]=0.0;
		}
		node=node->next;
	}	while(node!=node_root);
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
		tool->setXoldInd(xold);
		Ntoolbnodes=tool->createLists(dom);
		tool->setype(tool->current_framenode->type);
		tool->setboundary();
	}
	Cshear=1.5;// normal stress
	Cstrain=1.5;// shear stress
}
void	stepElastic(double dt, Domain *dom)
{	const	double	dtool=.5;
	Tool	*tool=dom->tool;
	if 
	(	tool!=NULL&&tool->first_framenode!=NULL
		&&(int)((runtime.current-runtime.step)/dtool)!=(int)(runtime.current/dtool)
		&&runtime.current>2.0
	)
	{//MOVE TOOL TO THE NEXT POSITION ALONG THE FRAME
		Frame
			*&first_framenode=tool->first_framenode,
			*&last_framenode=tool->last_framenode,
			*&current_framenode=tool->current_framenode;
		tool->push();
		if (current_framenode!=last_framenode)
			current_framenode=current_framenode->next;
		tool->setype(current_framenode->type);
		tool->pos(current_framenode->x);
		tool->setradius(current_framenode->radius);
		tool->cleanSurface();
		Ntoolbnodes=tool->createLists(dom);
	}
	MoveBoundary(dt,dom);
	dom->relax//Relax internal nodes
	(	ixold,
		ivoln,
		ivolc,
		0.5 //relaxation factor
	);
	runtime.current+=dt;
}

