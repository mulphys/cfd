/*
 * Incompressible Viscous flow
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
#include "poisson.h"
#include "templates.cc"
#include "incomflow.h"

extern void getvec(int loc, DCell *c, double *x, double *v);
namespace ModIncomFlow
{
	enum VariablesModIncomFlow
	{
		nvolume=0,
		cvolume,
		diagcoeffs,
		velocity,
		delatation,
		velgradc,
		condif,
		pressure,
///		dpressure,//pressure correction
		prsource,//pressure source term (derived from prhs)
///		pmass,//particle mass
///		pvelocity,//particle velocity
		maxVarModIncomFlow
	};
	const	int		ndim2=DIM*DIM;
	double
		dti,//inverse time-step
		density=1.0,densi=1./density,
		visc=1.50e-0,///-8//Viscosity, [length^2/time]
		Pamb=1.00e-0,//Ambient pressure [mass/(length*time^2)]
		inlet_pressure=1.1*Pamb;//Inlet pressure
	DCellList	*inlet_cell_root=NULL,
		*inlet_cell;

/*************************** FLOW SOLVER  ******************/
void	ConvDiff
(
	int	ivesr,
	int	ivoln,
	int	ivolc,
	int	iveln,
	int	ivgrc,
	Domain *dom
)
{	Variable	*variable=dom->variable;
	DNode	*node_root=dom->dnode_root,
		*node=node_root;
	DCell	*cell_root=dom->dcell_root,
		*cell=cell_root;
	if (node_root==NULL)return;
	node=node_root;
	do//Set source terms to zero
	{	double
			*var=node->var,
			*dvs=var+ivesr;
		for (int i=0; i<DIM; i++)
		{
			dvs[i]=0.0;
		}
		node=node->next;
	}	while(node!=node_root);
	if (cell_root==NULL)return;
	cell=cell_root;
	do //Convection-Diffusion loop
	{	DNode	**verts=cell->vert;
		DCell	**neibs=(DCell**)cell->neib;
		ElementStatus	*facetypes=cell->facetype;
		double	*x0,*x1,
			*varc=cell->var,
			volc=varc[ivolc],
			volci=1.0/(volc),
			*vgrc=varc+ivgrc,
			newvgrc[ndim2],
			e[Nv][DIM];//edges
		EDGES(e,verts);
		//Initialize
		for (int i=0; i<ndim2; i++)
			newvgrc[i]=0.0;
		//Face Areas and Face Fluxes
		for (int ivert=0; ivert<Nv; ivert++)
		{	ElementStatus	facetype=facetypes[ivert];
			DNode	*vert=verts[ivert];
			DCell	*neib=neibs[ivert];
			int vertype=vert->type;
			double
				area,areai,//face area and it's inverse
				dvol,
				a[DIM],//face area vector
				velf[DIM],velcf[DIM],//face velocity
				dves[DIM],//momentum entering the cell
				xface[DIM],//face-center coordinates
				vgrf[ndim2],//face velocity gradient
				*varn=vert->var,//variables at the nodes
				voln=varn[ivoln],volni=1./voln,
				*veln=varn+iveln,//velocity
				*vesr=varn+ivesr;//momentum change
			if(facetype==composite)
			{	BFaceList	*face=(BFaceList*)neibs[ivert];
				facetype=face->type;
			}
			//Face area
			AREAVEC(ivert,e,a);
			area=LENGTH(a);
#ifdef DEBUG
			if (area<SMALL){ERROR("AdvanceFlow: Area too small");}
#endif
			areai=1.0/area;
			//Face-center averaged vel and velgrad 
			for (int i=0; i<DIM; i++)
			{	int	ii=i*DIM;
				double
					*vgrci=newvgrc+ii,
				  vcenter=0.0;
				for (int jv=1; jv<Nv; jv++)
				{	int	jvert=(ivert+jv)%Nv;
					double 
						*var=verts[jvert]->var;
					vcenter+=var[iveln+i];
				}
				vcenter/=(double)Nfv;
				//Cell-Center Velocity gradient
				for (int j=0; j<DIM; j++)
					vgrci[j]+=vcenter*a[j];
				velcf[i]=vcenter;//cell-center velocity
			}
			if (facetype==internal)
			{
				unsigned int	ivneib,nv;//neighbor vertex opposite to the current face
				DNode	**neibverts=neib->vert,*neibvert;
				DCell	**neibneibs=(DCell**)neib->neib;
				double
					*neibvarn,
					neibvoln,
					*neibveln,
					*neibvarc=neib->var,
					*neibvgrc=neibvarc+ivgrc,
					neibvolc=neibvarc[ivolc],
					totvolc=volc+neibvolc,
					totvolci=1./totvolc,
					velmid,velup;
				ivneib=(int)(mask[ivert]&cell->index)>>((uint)ivert<<1);
				neibvert=neibverts[ivneib];
				neibvarn=neibvert->var;
				neibvoln=neibvarn[ivoln];
				neibveln=neibvarn+iveln;
				//Face-velocity (central differencing)
			if(neibvert->type==outlet)
			for (int i=0; i<DIM; i++)
			{
				velf[i]=veln[i];
				///velf[i]=0.5*(veln[i]+neibveln[i]);
				for (int j=0; j<DIM; j++)
				{	int	k=i*DIM+j;
					//Face velocity gradient
					vgrf[k]=0.0;
				}
			}
			else
			{
				for (int i=0; i<DIM; i++)
				{
					velf[i]=0.5*(veln[i]+neibveln[i]);
					for (int j=0; j<DIM; j++)
					{	int	k=i*DIM+j;
						//Face velocity gradient
						vgrf[k]=totvolci*(vgrc[k]+neibvgrc[k]);
					}
				}
				//Hybrid face velocity 
				velmid=0.5*(LENGTH(velf));
				velup=SCLP(a,velf)*areai;
				if (velup>=velmid)
				{
					for (int i=0; i<DIM; i++)
						velf[i]=veln[i];
				}
				else
				if (velup<=-velmid)
				{
					for (int i=0; i<DIM; i++)
						velf[i]=neibveln[i];
				}
				}
			}
			else
			{//boundary face
				switch (facetype)
				{
				case presinlet:
					if(SCLP(a,velcf)>=0)
						for (int i=0; i<DIM; i++)
							velf[i]=0.0;
					else
						for (int i=0; i<DIM; i++)
							velf[i]=velcf[i];
					//Face velocity gradient
					for (int i=0; i<ndim2; i++)
						vgrf[i]=0.0;
				break;
				case outlet:
				case presoutlet:
					///if (SCLP(a,veln)>0.0)
					///{
					for (int i=0; i<DIM; i++)
							velf[i]=veln[i];/// velcf[i];
					///}
					///else
					///{	for (int i=0; i<DIM; i++)
					///		velf[i]=0.0;
					///}
					//Face velocity gradient
					for (int i=0; i<ndim2; i++)
						vgrf[i]=0.0;
				break;
				default:
					//Face velocity
					for (int i=0; i<DIM; i++)
						velf[i]=0.0;///Cwallfriction*velcf[i];
					//Face velocity gradient
					for (int i=0; i<ndim2; i++)
						vgrf[i]=volci*vgrc[i];///*Cwallfriction
				}
			}
			dvol=SCLP(a,velf);
			//Update momentum  term
			if(vertype==internal||vertype==boundary)
			{
				for (int i=0; i<DIM; i++)
				{	double
						*v=vgrf+i*DIM,
						ddif=visc*(SCLP(a,v)),
						dcon=dvol*velf[i],
						ddc=ddif-dcon;
					dves[i]=ddc;//momentum entering the face
					vesr[i]+=ddc;
				}
			}
			//Set boundary source terms
			if (facetype==boundary)
			{	//loop over opposite-face-nodes
				for (int jv=1; jv<Nv; jv++)
				{	int	jvert=(ivert+jv)%Nv;
					DNode	*vertj=verts[jvert];
					int	vertypej=vertj->type;
					double	*varj=vertj->var,
						*dvsj=varj+ivesr;
					if(vertypej==boundary)
					for (int i=0; i<DIM; i++)
					{
						dvsj[i]+=dves[i];
					}
				}
			}
		}//next ivert: end loop over cell-vertexes
		for (int i=0; i<ndim2; i++)
			vgrc[i]=newvgrc[i];
		cell=cell->next;
	}	while (cell!=cell_root);
	//Normalize 
	node=node_root;
	do
	{
		int nodetype=node->type;
		double
			*var=node->var,
			*dvs=var+ivesr,
			vol=var[ivoln],voli=1./vol;
		for (int i=0; i<DIM; i++)dvs[i]*=voli;
		node=node->next;
	}	while(node!=node_root);
}
void	AdvanceFlow
(
	double dt,
	Domain *dom
)
{	using namespace ModPoisson;
	Variable	*variable=dom->variable;
	int
		ivoln=variable[nvolume].loc,
		ivolc=variable[cvolume].loc,
		idiag=variable[diagcoeffs].loc,
		iveln=variable[velocity].loc,
		ivgrc=variable[velgradc].loc,
		idela=variable[delatation].loc,
		icond=variable[condif].loc,
		iprgr=icond,
		ipres=variable[pressure].loc,
		iprsr=variable[prsource].loc;
	double
		dti=1./dt;//inverse of the time step
	DNode	*node_root=dom->dnode_root,
		*node=node_root;
	DCell	*cell_root=dom->dcell_root,
		*cell=cell_root;

	//CONVECTION/DIFFUSION
	ConvDiff 
	(
		icond, // returns the convection/diffusion term
		ivoln, // nodal volumes
		ivolc, // cell volumes
		iveln, // nodal velocities
		ivgrc, // velocity gradient
		dom    // current domain
	);
	//Update velocity
	node=node_root;
	do
	{	int nodetype=node->type;
		double
			*var=node->var,
			*vel=var+iveln,
			*con=var+icond;
		if(nodetype==internal||nodetype==outlet)
		for (int i=0; i<DIM; i++)
		{
			vel[i]+=dt*con[i];
		}
		node=node->next;
	}	while(node!=node_root);

	//PRESSURE
	//Pressure source terms:
	dom->dVidXi(ivoln,icond,idela);//Delatation: Nabla.cond
	assemblePoissonSource //FEM-assembly of the pressure source terms
	(
		iprsr,// index to the pressure source variable to be assembled 
		idela,// index to the pressure source scalar 
		ivolc,// index to the cell-volume variables
		dom->dnode_root,
		dom->dcell_root
	);

///		assemblePoissonSourceGrad  ///<-- correcting poorly
///		(
///			iprsr,// index to the pressure source variable to be assembled 
///			iveln,// index to the pressure gradient source vector 
///			ivolc,// index to the cell-volume variables
///			dom->dnode_root,
///			dom->dcell_root
///		);

	//Poisson pressure solver:
	stepPoisson 
	(
		ipres,//index to the variable to be updated
		iprsr,//source term variable
		idela,//memory buffer for temporal variable
		      // at least scalar, but a vector is OK 
		idiag,//diagonal matrix elements
		ivoln,//index to nodal volumes
		ivolc,//index to cell volumes
		dom->dnode_root,
		dom->dcell_root
	);
	//Compute pressure gradient and update velocities	
	{	int iprgr=icond;//re-use velocity source array for 
		                //pressure-gradient storage 
		//Pressure gradient
		dom->dVdXi(ivoln,ipres,iprgr);
		//Update velocity for pressure
		node=node_root;
		do
		{	int	nodetype=node->type;
			double
				*var=node->var,
				*vel=var+iveln,
				*prg=var+iprgr;
			if(nodetype==internal)
			for (int i=0; i<DIM; i++)
			{
				vel[i]-=prg[i]*dt;
			}
			node=node->next;
		}	while(node!=node_root);
	}
}
/***************** ADVANCE PARTICLES ******************************/
///	#include "../particles.cc"
///	void	InjectParticle(Domain *dom)
///	{
///		int	ip,//index to the new particle
///			iveloc=dom->variable[velocity].loc;
///		double
///			xinj[DIM],dx,//injection point and dispersion
///			minj=1.0,//injection mass
///			vinj[DIM];//injection velocity
///			//Create injection coordinates array
///		Point	*origin=dom->origin,
///			*p;//current particle
///		DCell	*cell;
///		DNode	**verts;
///		if (inlet_cell==NULL) return;
///		cell=inlet_cell->cell;
///		verts=cell->vert;
///		for (int i=0; i<DIM; i++)
///		{	double	xcenter=0.0;
///			for (int iv=0; iv<Nv; iv++)
///				xcenter+=verts[iv]->x[i];
///			xcenter/=(double)Nv;
///			xinj[i]=xcenter;
///		}
///		dx=0.0;
///		for (int i=0; i<DIM; i++)
///		{	double	d=0.0;
///			for (int iv=0; iv<Nv; iv++)
///			{	d+=fabs(verts[iv]->x[i]-xinj[i]);
///			}
///			d/=(double)Nv;
///			dx+=d*d;
///		}
///		dx=0.5*sqrt(dx);
///		for (int i=0; i<DIM; i++)
///		{	double	r=2.0*RND-1.0;
///			xinj[i]+=r*dx;
///		}
///		inlet_cell=inlet_cell->next;
///		ip=dom->putp(xinj);
///		dom->setsclp(ip,pmass,minj);
///		getvec(iveloc,cell,xinj,vinj);
///		dom->setvecp(ip,pvelocity,vinj);
///		p=origin+ip;
///		p->cell=cell;
///	}
///	void	AdvanceParticles(double dt, Domain *dom)
///	{
///		int	iveloc=dom->variable[velocity].loc;
///		DCell	*root_cell=dom->dcell_root;
///		Point	*root=dom->origin,*&first=dom->first,*&last=dom->last,
///			*p=first;
///		double
///			*M=dom->variable[pmass].val,
///			*X=dom->coordinates[points].val,
///			*V=dom->variable[pvelocity].val;
///		if(dom->getNoPoints()==0)return;
///		do
///		{	int	ip=p-root,imp=DIM*ip;
///			ElementStatus	bc;//boundary conditions
///			double
///				*m=M+ip,
///				*x=X+imp,//particle position
///				*v=V+imp,//particle velocity
///				fvel[DIM];//flow velocity at particle location
///			DCell	*pcell=(DCell*)p->cell;//cell containing the particle
///			if (*m<0.0) goto next_loop;
///			ZERO3(fvel);
///			getvec(iveloc,pcell,x,fvel);
///			MoveParticle(root,ip,dt,*m,fvel,x,v);
///			if((bc=dom->findcell(x,pcell,(DCell*)p->cell))!=internal)
///			switch(bc)
///			{	case boundary:
///				*m=-1.0;
///				break;
///				default:
///				dom->delp(p);
///				break;
///			}
///			next_loop: p=p->next;
///		}	while(p!=last->next);
///	}
};//END NAMESPACE 

using namespace ModIncomFlow;

int	getNvarModIncomFlow()
{
	return maxVarModIncomFlow;
}
void	defVarModIncomFlow
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
		case diagcoeffs:
			strcpy(name,"DiagonalCoeffs");
			type=nodes;
			dim=1;
			break;
		case velocity:
			strcpy(name,"Velocity");
			type=nodes;
			dim=DIM;	
			break;
		case delatation:
			strcpy(name,"Delatation");
			type=nodes;
			dim=1;
			break;
		case velgradc:
			strcpy(name,"CellVelocityGradient");
			type=cells;
			dim=DIM*DIM;
			break;
		case condif:
			strcpy(name,"ConvectionDiffusion");
			type=nodes;
			dim=DIM;
			break;
		case pressure:
			strcpy(name,"Pressure");
			type=nodes;
			dim=1;
			break;
		case prsource:
			strcpy(name,"PressureSource");
			type=nodes;
			dim=1;
			break;
///			case dpressure:
///				strcpy(name,"PressureCorrection");
///				type=nodes;
///				dim=1;
///				break;
///			case pmass:
///				strcpy(name,"ParticleMass");
///				type=points;
///				dim=1;
///				break;
///			case pvelocity:
///				strcpy(name,"ParticleVelocity");
///				type=points;
///				dim=DIM;
///				break;
		default:
			fprintf
			(	stderr,
				"Specification of variable %d domain %d (%s) is incomplete\n",
				ivar+1,idomain+1,dom->name
			);
			return;
	}
}
void	initVarModIncomFlow
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
void	initModIncomFlow(Domain *dom)
{	using namespace ModPoisson;
	Variable	*variable=dom->variable;
	int
		ivoln=variable[nvolume].loc,
		ivolc=variable[cvolume].loc,
		idiag=variable[diagcoeffs].loc,
		iveln=variable[velocity].loc,
		icond=variable[condif].loc,
		iprsr=variable[prsource].loc,
		ipres=variable[pressure].loc;
	DNode	*node_root=dom->dnode_root,
		*node=node_root;
	DCell	*cell_root=dom->dcell_root,
		*cell=cell_root;
	dom->setCellVolumes(cvolume);
	dom->setNodeVolumes(nvolume,cvolume);
	setDiagCoeff(ivolc,idiag,node_root,cell_root);
	dom->setBoundaryVertexes();
	dom->indexNeibNodes();
	dom->setBoundaryCellFlags();
	if (node_root==NULL) return;
	do//Initialize
	{	int	type=node->type;
		double
			*var=node->var,
			vol=var[ivoln],voli=1./vol,
			*vel=var+iveln,
			*dvs=var+icond,
			*prs=var+ipres,
			*psr=var+iprsr;
		*prs=Pamb;
		for (int i=0; i<DIM; i++)
		{	int	ii=DIM*i;
			vel[i]=0.0;
			dvs[i]=0.0;
		}
		if(type==inlet)vel[2]=1.0;///DDD
		node=node->next;
	}	while(node!=node_root);
	//Initialize particles
	dom->initp();
	//Create the list of injection cells
	if (dom->dcell_root!=NULL)
	{	int	ninj=0;
		DCell
			*cell_root=dom->dcell_root,
			*cell=cell_root;
		if(inlet_cell_root==NULL)
		do
		{	ElementStatus	*facetype=cell->facetype;
			double	rinj;//coordinates of the injection point
			DNode	**verts=cell->vert;
			//Compute injection-cell coordinates
			rinj=0.0;
			for (int i=0; i<DIM; i++)
			{	double xc=0.0;//center coordinate
				if (i==2)continue;
				for (int iv=0; iv<Nv; iv++)
				{	double	*x=verts[iv]->x;
					xc+=x[i];
				}
				xc/=(double)Nv;
				rinj+=xc*xc;
			}
			rinj=sqrt(rinj);
			if (rinj>2.0)goto next_loop;//restrict injection cells
			for (int iface=0; iface<Nf; iface++)
			if (facetype[iface]==presinlet)
			{	ninj++;
				if (inlet_cell_root==NULL)
				{	inlet_cell_root=new DCellList;
					inlet_cell=inlet_cell_root;
				}
				else
				{	inlet_cell->next=new DCellList;
					inlet_cell=inlet_cell->next;
				}
				inlet_cell->cell=cell;
				inlet_cell->next=inlet_cell_root;//make it a loop
				break;
			}
			next_loop: cell=cell->next;
		}	while(cell!=cell_root);
		if(option.verbose)printf("Number of injection points: %d\n",ninj);
	}
	inlet_cell=inlet_cell_root;
#ifdef DEBUG
	dom->checkGrid(dom->dnode_root,dom->dcell_root);
#endif
}
void	stepModIncomFlow(double dt, Domain *dom)
{
//	static const	double	dtinj=1.0;//injection time interval
//	if ((int)(time.current-time.step)/dtinj!=(int)(time.current/dtinj))
//		InjectParticle(dom);
	AdvanceFlow
	(	
		dt,dom
	);
//	AdvanceParticles(dt,dom);
	runtime.current+=dt;
}
