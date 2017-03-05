
/* 
 *  MODEL: Poisson
 *  Solves poisson equation 
 *	Using FEM
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
#include "gui.h"
#include "poisson.h"
#include "templates.cc"

namespace ModPoisson
{
	enum VariablesModPoisson
	{
		volumen=0,
		volumec,
		diagcoeffs,
		scalar,//variable to be solved
		force,//right-hand side of the poisson eq.
		source,//source term A (\ref{aisum})
		diffusion,//source term C (\ref{cisum})
		maxVarModPoisson
	};
void	setDiagCoeff
(	//Sets the inverted diagonal coefficients c1ii
	//for the stiffness matrix \ref{}
	int	ivolc,
	int	idiag,
	DNode	*node_root,
	DCell	*cell_root
)
{	const double	coeff=1./36.;
	DNode	*node=node_root;
	DCell	*cell=cell_root;
	if (node_root==NULL) return;
	do//Zero out the output terms
	{	double *var=node->var;
		var[idiag]=0.0;
		node=node->next;
	}	while(node!=node_root);
	if (cell_root==NULL)return;
	cell=cell_root;
	do //Compute the diagonal coefficients
	{	DNode	**verts=cell->vert;
		double 
			volc=cell->var[ivolc],volci,
			e[Nv][DIM];//tetrahedral edge-vectors					
#ifdef DEBUG
		if(volc<SMALL)ERROR("Volume too small in setDiagCoeff");
#endif
		volci=coeff/volc;
		EDGES(e,verts);
		for (int iv=0; iv<Nv; iv++)
		{	int
				iv0=(iv+1)%Nv,
				iv1=(iv+2)%Nv,
				iv2=(iv+3)%Nv;
			double
				*var=verts[iv]->var,
				area,a[DIM];//face area and area-vector
			//Face area vector
			VECP(a,e[iv1],e[iv2]);
			area=LENGTH(a);
			var[idiag]+=area*area*volci;
		}
		cell=cell->next;
	}	while(cell!=cell_root);
	node=node_root;
	do//Invert the diag coeffs
	{	double *var=node->var,
			*diag=var+idiag;
#ifdef DEBUG
		if(*diag<SMALL)
			ERROR("Diagonal coefficients too small");
#endif
		*diag=1./(*diag);//Store as inverted coeffs
		node=node->next;
	}	while(node!=node_root);
}
/* 
 * Assemble source terms for
 * FEM discretized Poissson equation
 *
 */
void	assemblePoissonSource
(
	int	iscr,// index to the source variable to be assembled
	int	ifor,// index to the force variable 
	int	ivolc,//index to the cell-volume variables
	DNode	*node_root,
	DCell	*cell_root
)
{	DNode	*node=node_root;
	DCell	*cell=cell_root;
	//Set the output variable to zero
	if(node_root==NULL)return;
	node=node_root;
	do
	{	
		node->var[iscr]=0.0;
		node=node->next;
	}	while(node!=node_root);
	//Assemble the output variable
	if(cell_root==NULL)return;
	cell=cell_root;
	do
	{	double
			volc=cell->var[ivolc],
			c0ii=0.10*volc,//\ref{}
			c0ij=0.05*volc;
		DNode	**verts=cell->vert;
		for(int iv=0;iv<Nv;iv++)
		{	DNode	*nodi=verts[iv];
			double	*vari=nodi->var,
				dfor=c0ii*vari[ifor];
			for(int jv=1;jv<Nv;jv++)
			{	DNode	*nodj=verts[(iv+jv)%Nv];
				double	*varj=nodj->var;
				dfor+=c0ij*varj[ifor];
			}
			vari[iscr]+=dfor;
		}
		cell=cell->next;
	}	while(cell!=cell_root);
}
void	assemblePoissonSourceGrad
(
	int	iscr,// index to the source variable to be assembled
	int	ifor,// index to the force vector variable 
	int	ivolc,//index to the cell-volume variables
	DNode	*node_root,
	DCell	*cell_root
)
{	DNode	*node=node_root;
	DCell	*cell=cell_root;
	//Set the output variable to zero
	if(node_root==NULL)return;
	node=node_root;
	do
	{
		node->var[iscr]=0.0;
		node=node->next;
	}	while(node!=node_root);
	//Assemble the output variable
	if(cell_root==NULL)return;
	cell=cell_root;
	do
	{	double
			volc=cell->var[ivolc],
			cij=0.5*volc,
			e[Nv][DIM];//tetrahedral edge-vectors
		DNode	**verts=cell->vert;
		EDGES(e,verts);//computing edge-vectors
		for(int iv=0;iv<Nv;iv++)
		{	DNode	*nodi=verts[iv];
			double	*vari=nodi->var,
				h[DIM],h0,h2,
				*e0=e[(iv+1)%Nv],
				*e1=e[(iv+2)%Nv],
				*e2=e[(iv+3)%Nv],
				dfor;
			//Compute height vector
			VECP(h,e1,e2);
			h0=LENGTH(h);
			for(int i=0;i<DIM;i++)h[i]/=h0;
			h0=SCLP(e0,h);
			for(int i=0;i<DIM;i++)
				h[i]*=h0;
			h2=h0*h0;
			dfor=0.0;
			for(int jv=0;jv<Nv;jv++)
			{	DNode	*nodj=verts[(iv+jv)%Nv];
				double	*varj=nodj->var,
					*f=varj+ifor;
				dfor+=SCLP(h,f);
			}
			vari[iscr]+=cij*dfor*volc/h2;
		}
		cell=cell->next;
	}	while(cell!=cell_root);
}
/* 
 * Assemble boundary source terms for
 * FEM discretized Poissson equation
 * with non-zero boundary gradient 
 */
void	assemblePoissonSourceNeuman
// Is called after assemblePoisssonSource
// if the non-zero gradinet boundary conditions
// are present
(
	int	isrc,// index to the boundary source variable to be assembled
	int	ifor,// index to the force variable 
	int	ibnd,// index to the boundary vector variable 
//-	int	idiag,// index to the diagonal coefficients 
	int	ivolc,//index to the cell-volume variables
	BFaceList	*bface_root,
	DNode	*node_root,
	DCell	*cell_root
)
{	const	double
		onesixth=1./6.,
		onesixth2=1./36.;
	DNode	*node=node_root;
	DCell	*cell=cell_root;
	BFaceList	*bface=bface_root;
	if(node_root==NULL)return;
	//Set to zero
	node=node_root;
	do
	{	double	*var=node->var;
		var[isrc]=0.0;
//-			var[idiag]=0.0;
		node=node->next;
	}	while(node!=node_root);
	//Assemble the source terms variable
	if(cell_root==NULL)return;
	cell=cell_root;
	do//Set the source terms
	{	double
			volc=cell->var[ivolc],
			c0ii=0.10*volc,//\ref{}
			c0ij=0.05*volc;
		DNode	**verts=cell->vert;
		for(int iv=0;iv<Nv;iv++)
		{	DNode	*nodi=verts[iv];
			double	*vari=nodi->var,
				dfor=c0ii*vari[ifor];
			for(int jv=0;jv<Nv1;jv++)
			{	DNode	*nodj=verts[(iv+jv+1)%Nv];
				double	*varj=nodj->var;
				dfor+=c0ij*varj[ifor];
			}
			vari[isrc]+=dfor;
		}
		cell=cell->next;
	}	while(cell!=cell_root);
//-		cell=cell_root;
//-		do //Assemble the diagonal coefficients
//-		{	DNode	**verts=cell->vert;
//-			double 
//-				volc=cell->var[ivolc],volci,
//-				e[Nv][DIM];//tetrahedral edge-vectors					
//-	#ifdef DEBUG
//-			if(volc<SMALL)ERROR("Volume too small in setDiagCoeff");
//-	#endif
//-			volci=onesixth2/volc;
//-			EDGES(e,verts);
//-			for (int iv=0; iv<Nv; iv++)
//-			{	int
//-					iv0=(iv+1)%Nv,
//-					iv1=(iv+2)%Nv,
//-					iv2=(iv+3)%Nv;
//-				double
//-					*var=verts[iv]->var,
//-					area,a[DIM];//face area and area-vector
//-				//Face area vector
//-				VECP(a,e[iv1],e[iv2]);
//-				area=LENGTH(a);
//-				var[idiag]+=area*area*volci;
//-			}
//-			cell=cell->next;
//-		}	while(cell!=cell_root);
	//Process the boundary conditions
	if(bface_root==NULL)return;
	bface=bface_root;
	do
	{	DCell	*cell=bface->cell;
		DNode	**verts=cell->vert;
		int	iface=bface->iface;
		double	area=bface->area,
			*norm=bface->norm,
			bii=onesixth*area,
			bij=0.25*area;
		for(int iv=0;iv<Nfv;iv++)
		{	int	ivert=(iface+iv+1)%Nv,
				i0=1;
			DNode	*nodi=verts[ivert];
			double	*vari=nodi->var,
				*bndi=vari+ibnd,
				src=bii*(SCLP(norm,bndi));//normal component of the boundary vector
			///vari[idiag]-=bii*(SCLP(norm,bndi));//normal component of the boundary vector
			///DDD: diagbnd should be of the type boundary_nodes with the connectivity to bounary-vertex array
			for (int jv=0; jv<Nfv1; jv++)
			{	int	jvert=(ivert+jv+i0)%Nv;
				if(jvert==iface)
				{	i0++;jv--;continue;
				}
				DNode	*nodj=verts[jvert];
				double	*varj=nodj->var,
					*bndj=varj+ibnd;
				src+=bij*(SCLP(norm,bndj));
			}
			vari[isrc]-=src;
		}
		bface=bface->next;
	}	while(bface!=bface_root);
}
/*
 * One step of Gauss-Seidel iteration of a Poisson FEM solver
 *
 */
void	stepPoisson
(
	int	iscl,//index to the variable to be updated
	int	isrc,//source term variable-index
	int	idif,//diffusion variable-index
	int	idiag,//diagonal matrix elements
	int	ivoln,
	int	ivolc,
	DNode	*node_root,
	DCell	*cell_root
)
{	const	int	max_iter=8; //maximum numbar of sub-cycles
	const	double	coeff=1./36.;//1/(4*9), because we use 
	                             //double areas
	if(node_root==NULL)return;
	if(cell_root==NULL)return;
	for(int iter=0;iter<max_iter;iter++)
	{	DCell	*cell;
		DNode	*node=node_root;
		do//Set the diffusion terms to zero
		{
			node->var[idif]=0.0;
			node=node->next;
		}	while(node!=node_root);
		cell=cell_root;
		do//Compute the diffusion term
		{	DNode	**verts=cell->vert;
			double	volc=cell->var[ivolc],volci=coeff/volc,
				e[Nv][DIM],//tetrahedral edge-vectors
				a[Nv][DIM];//areas 
			EDGES(e,verts);
			//Compute areas
			for(int iv=0;iv<Nv;iv++)
			{
				DNode	*nodi=verts[iv];
				double	*av=a[iv];//area vector
				AREA2(iv,e,av);
			}
			//Compute diffusion terms
			for(int iv=0;iv<Nv;iv++)
			{	DNode	*nodi=verts[iv];
				double	*vari=nodi->var,
					*ai=a[iv],
					*diffi;
					diffi=vari+idif;
				for(int jjv=0;jjv<Nv1;jjv++)
				{	int	jv=(jjv+iv+1)%Nv;
					DNode	*nodj=verts[jv];
					double	*varj=nodj->var,
						*aj=a[jv];
					*diffi+=(SCLP(ai,aj))*volci*varj[iscl];
				}
			}
			cell=cell->next;
		}	while(cell!=cell_root);
		node=node_root;
		do//Update the nodal variable
		{	double	*var=node->var;
			int	nodetype=node->type;
			if(nodetype!=presinlet&&nodetype!=presoutlet&&nodetype>maxElementStatus)
				var[iscl]=-var[idiag]*(var[isrc]+var[idif]);
			node=node->next;
		}	while(node!=node_root);
	}
}
};

using namespace ModPoisson;

int	getNvarModPoisson()
{
	return maxVarModPoisson;
}
void	defVarModPoisson
(
	int	ivar,
	char	*name,
	int	&type,
	int	&dim,
	Domain	*dom
)
{	int idomain=dom-domain_root;
	*name='\0';
	switch(ivar)
	{
		case volumec:
			strcpy(name,"CellVolume");
			type=cells;
			dim=1;
			break;
		case volumen:
			strcpy(name,"NodeVolume");
			type=nodes;
			dim=1;
			break;
		case diagcoeffs:
			strcpy(name,"DiagonalCoeffs");
			type=nodes;
			dim=1;
			break;
		case scalar:
			strcpy(name,"Temperature");
			type=nodes;
			dim=1;
			break;
		case force:
			strcpy(name,"HeatSource");
			type=nodes;
			dim=1;
			break;
		case source: //\ref{aisum}
			strcpy(name,"Source");
			type=nodes;
			dim=1;
			break;
		case diffusion: //\ref{cisum}
			strcpy(name,"Diffusion");
			type=nodes;
			dim=1;
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
void	initVarModPoisson
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
void	initModPoisson(Domain *dom)
{
	Tool	*&tool=dom->tool;
	Variable	*variable=dom->variable;
	int //Indexes to the variables
		ivoln=variable[volumen].loc,
		ivolc=variable[volumec].loc,
		idiag=variable[diagcoeffs].loc,
		iscl=variable[scalar].loc,
		ifor=variable[force].loc,
		isrc=variable[source].loc,
		idif=variable[diffusion].loc;
	DNode	*node_root=dom->dnode_root,
		*node=node_root;
	DCell	*cell_root=dom->dcell_root,
		*cell=cell_root;
	if(node_root==NULL) return;
	do//Initialize nodal variables
	{	double
			*var=node->var,
			*scl=var+iscl,
			*frc=var+ifor,
			*src=var+isrc,
			*dif=var+idif;
		*scl=0.0;
		*src=0.0;
		*frc=0.0;
		*dif=0.0;
		//	for (int i=0; i<DIM; i++)
		//	{	int	ii=DIM*i;
		//		vel[i]=0.0;
		//	}
		node=node->next;
	}	while(node!=node_root);
	//Initialize face variables
	//See mod_softbody.cc
	//Initialize points
	///dom->initp();
	///Create particles
	///if(dom->getMaxNoPoints()-dom->getNoPoints()<1) 
	///	ERROR("Not enough particles\n");
	if (tool!=NULL)
	{
		tool->mode.active=tool->mode.active?0:1;
		if (option.verbose)printf("Tool active = %s\n",tool->mode.active?"YES":"NO");
		tool->current_framenode=tool->first_framenode;
		if(tool->first_framenode==NULL)ERROR("No tool framenodes specified");
		tool->pos(tool->current_framenode->x);
		tool->setradius(tool->current_framenode->radius);
		tool->setrefinement(tool->current_framenode->refinement);
		tool->setVolInd(volumen);
		//tool->setXoldInd(xold);
		tool->setype(tool->current_framenode->type);
		tool->setboundary();
	}
	dom->setCellVolumes(volumec);
	dom->setNodeVolumes(volumen,volumec);
	//dom->setNodeVolumes(volumen);
	if(option.verbose)
	{	//Check volumes
		DNode
			*node_root=dom->dnode_root,	
			*node=node_root;
		double	totvoln=0.0,totvolc=0.0;
		if (node_root!=NULL)
		do//Zero out source terms
		{	double
				*var=node->var;
			totvoln+=var[ivoln];
			totvolc+=var[ivolc];
			node=node->next;
		}	while(node!=node_root);
		printf("Total volumes: nodal=%g, cell-based\n",0.25*totvoln,totvolc);FLUSH;
	}
	//Compute the c1ii coefficients \ref{}
	setDiagCoeff(ivolc,idiag,node_root,cell_root);
	dom->setBoundaryVertexes();
	do//Set initial and boundary conditions
	{	double
			*x=node->x,
			*var=node->var,
			*scl=var+iscl,
			*frc=var+ifor,
			*src=var+isrc,
			*dif=var+idif;
		*scl=0.0;
		*src=0.0;
		*frc=0.0;
		*dif=0.0;

		/*add some point force term */
		/*
		double source_radius=2.0;
		if(fabs(node->x[0])<source_radius
			&&fabs(node->x[1])<source_radius
			&&fabs(node->x[2]+3.0)<source_radius)	
		{
			//*frc=1.0; //sink
			*frc=-1.0; //source
		}
		*/
		/**********************/

		//	for (int i=0; i<DIM; i++)
		//	{	int	ii=DIM*i;
		//		vel[i]=0.0;
		//	}
		//Set the forces: delta function in the center
		//if(LENGTH(x)<0.1)*frc=1.0;
		//Set the boundary conditions
		if(node->type==(int)presinlet)*scl=1.0;
		if(node->type==(int)presoutlet)*scl=-1.0;
		node=node->next;
	}	while(node!=node_root);
}
void	stepModPoisson(double dt, Domain *dom)
{	const	double	dtool=1.0;
	Tool	*tool=dom->tool;
	Variable	*variable=dom->variable;
	int //Indexes to the variables
		ivoln=variable[volumen].loc,
		ivolc=variable[volumec].loc,
		idiag=variable[diagcoeffs].loc,
		iscl=variable[scalar].loc,
		ifor=variable[force].loc,
		isrc=variable[source].loc,
		idif=variable[diffusion].loc;
	//Assemble source terms
	assemblePoissonSource
	(
		isrc,// index to the source variable to be assembled 
		ifor,// index to the force variable
		ivolc,//index to the cell-volume variables
		dom->dnode_root,
		dom->dcell_root
	);
	//Update the Scalar 
	stepPoisson
	(
		iscl,//index to the variable to be updated
		isrc,//source term variable-index
		idif,//diffusion variable-index
		idiag,//diagonal matrix elements
		ivoln,
		ivolc,
		dom->dnode_root,
		dom->dcell_root
	);
	runtime.current+=dt;
}
