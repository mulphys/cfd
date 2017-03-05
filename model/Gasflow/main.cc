/* 
 *  MODEL: Viscous flow with particles
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
#include "templates.cc"
#include "gasflow.h"

//#define NO_INLET_OUTFLOW

extern void getvec(int loc, DCell *c, double *x, double *v);
namespace ModGasFlow
{
	enum VariablesModGasFlow
	{
		nvolume=0,
		cvolume,
		mass,
		addmass,
		velocity,
#ifndef LOCALDIFFUSION
		velgradc,
#endif
		dmomentum,
		pressure,
		dpressure,
///		pmass,//particle mass
///		pvelocity,//particle velocity
		maxVarModGasFlow
	};
	const	int		ndim2=DIM*DIM;
	double
///DDD		Tamb=300.0,///300,//Ambient temperature [K]
		//Marginally stable with time.step=1.0e-2: diverged at time=220
		//time.step=1.0e-3: time=30 looked good
		//Rgas=2.87e+2,//Gas-constant of air [J/(kmol K)
		//visc=1.00e-5,//Air viscosity [m^2/s]
		//Pamb=1.00e+5,//Ambient pressure [Pa]
		//Pin=1.1*Pamb,//Inlet pressure
		//Centimeters, grams, time.step=1.0e-3: oscillatory at time=11
		//Rgas=2.86e-2,///Gas-constant of air [J/(kmol K)
		//Pamb=1.e4,///Ambient pressure [Pa]
		//Pin=1.01*Pamb,//Inlet pressure
		//Centimeters, grams, milliseconds: time.step=1.0e-2
		//Oscillatory at time=1.0, time.step=1.e-4:
		//	Rgas=10.0,///287,//Gas-constant of air [J/(kmol K)
		//	Pamb=100.0,///5,//Ambient pressure [Pa]
		//	Pin=1.001*Pamb,//Inlet pressure
		//milliseconds, meters, kilograms
		//Rgas=2.87e-4,//Gas-constant of air [work/(mass*temp)
		//visc=1.50e-0,///-8//Viscosity, [length^2/time]
		//Pamb=1.00e-1,//Ambient pressure [mass/(length*time^2)]
		Rgas=1.0,Tamb=1.0,//Gas-constant of air [work/(mass*temp)
		visc=1.50e-0,///-8//Viscosity, [length^2/time]
		Pamb=1.00e-0,//Ambient pressure [mass/(length*time^2)]
		Pin=1.2*Pamb,//Inlet pressure
		Tin=Tamb,
///		Cwall=0.5;//Wall-friction coefficient [mass/(length^2*time)]: ParticleFlow
		Cwall=10.;//Wall-friction coefficient [mass/(length^2*time)]
	DCellList	*inlet_cell_root=NULL,
		*inlet_cell;

double	InletDensity(double time)
{	double	den;
	//if(time<200.0)
	den=Pin/(Rgas*Tin);
	//else	den=Pamb/(Rgas*Tamb);
	return	den;
}
double	OutletDensity(double time)
{
	return	Pamb/(Rgas*Tamb);
}
//inline double	Pressure(double den)
/*delete "inline" because pgCC won't pass it */
double	Pressure(double den)
{
	return	den*Rgas*Tamb;
}
double	Pressure(double den, double temp)
{
	return	Rgas*den*temp;
}
/*************************** ADVANCE FLOW ******************/
//#define LIDDRIVENCAVITY
void	AdvanceFlow
(
	double inlet_density,
	double outle_density,
	double dt,
	Domain *dom
)
{
	Variable	*variable=dom->variable;
	int
		ivoln=variable[nvolume].loc,
		ivolc=variable[cvolume].loc,
		imasn=variable[mass].loc,
		iadmn=variable[addmass].loc,
		iveln=variable[velocity].loc,
#ifndef LOCALDIFFUSION
		ivgrc=variable[velgradc].loc,
#endif
		idmmn=variable[dmomentum].loc,
		iprsn=variable[pressure].loc,
		idprn=variable[dpressure].loc;
	DNode	*node_root=dom->dnode_root,
		*node=node_root;
	DCell	*cell_root=dom->dcell_root,
		*cell=cell_root;
	if (node_root==NULL)return;
	node=node_root;
	do//Set source terms to zero
	{	double
				*var=node->var,
				*damn=var+iadmn,
				*dmm=var+idmmn,
				*dpr=var+idprn;
		*damn=0.0;
		for (int i=0; i<DIM; i++)
		{
			dmm[i]=dpr[i]=0.0;
		}
		node=node->next;
	}	while(node!=node_root);
	if (cell_root==NULL)return;
	cell=cell_root;
	do //Set source terms: Loop over cells
	{	DNode	**verts=cell->vert;
		DCell	**neibs=(DCell**)cell->neib;
		ElementStatus	*facetypes=cell->facetype;
		double	*x0,*x1,
			*varc=cell->var,
			volc=varc[ivolc],
			volci=1.0/(volc),
#ifndef LOCALDIFFUSION
			*vgrc=varc+ivgrc,
			newvgrc[ndim2],
#endif
			e[Nv][DIM];//edges
		EDGES(e,verts);
#ifndef LOCALDIFFUSION
		//Initialize
		for (int i=0; i<ndim2; i++)
			newvgrc[i]=0.0;
#endif
		//Face Areas and Face Fluxes
		for (int ivert=0; ivert<Nv; ivert++)
		{	ElementStatus	facetype=facetypes[ivert];
			DNode	*vert=verts[ivert];
			DCell	*neib=neibs[ivert];
			double
				area,areai,//face area and it's inverse
				dm,//mass entering the cell
				dvol,
				denf,//face density
				a[DIM],//face area vector
				velf[DIM],velcf[DIM],//face velocity
				dmom[DIM],//momentum entering the cell
				xface[DIM],//face-centern coordinates
#ifdef	LOCALDIFFUSION
				dvdxface[DIM],//face-center derivative
				*neibcoords,//coordinates of the neighbor point
				*neibveldif,//velocity at the neighbor point
#else
				vgrf[ndim2],//face velocity gradient
#endif
				*varn=vert->var,//variables at the nodes
				voln=varn[ivoln],volni=1./voln,
				masn=varn[imasn],
				denn=masn*volni,
				prsn=varn[iprsn],

				prsf,//pressure at the face
				dpf,//pressure drop accross the face
				*admn=varn+iadmn,//added mass
				*veln=varn+iveln,//velocity
				*dmmn=varn+idmmn,//momentum change
				*dprn=varn+idprn;//pressure change
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
#ifndef LOCALDIFFUSION
					*vgrci=newvgrc+ii,
#endif
				  vcenter=0.0;
				for (int jv=1; jv<Nv; jv++)
				{	int	jvert=(ivert+jv)%Nv;
					double 
						*var=verts[jvert]->var;
					vcenter+=var[iveln+i];
				}
				vcenter/=(double)Nfv;
#ifndef LOCALDIFFUSION
				//Cell-Center Velocity gradient
				for (int j=0; j<DIM; j++)
					vgrci[j]+=vcenter*a[j];
#endif
				velcf[i]=vcenter;//cell-center velocity
			}
			if (facetype==internal)
///			if (neib!=NULL)///DDD
			{
				unsigned int	ivneib,nv;//neighbor vertex opposite to the current face
				DNode	**neibverts=neib->vert,*neibvert;
				DCell	**neibneibs=(DCell**)neib->neib;
				double
					*neibvarn,
					neibvoln,
					neibmasn,
					neibdenn,
					*neibveln,
					*neibvarc=neib->var,
#ifndef LOCALDIFFUSION
					*neibvgrc=neibvarc+ivgrc,
#endif
					neibvolc=neibvarc[ivolc],
					totvolc=volc+neibvolc,
					totvolci=1./totvolc,
					velmid,velup;
#ifdef COMPUTE_NEIB_VERTEX_INDEXES
				for (ivneib=0; ivneib<Nv; ivneib++)
					if (neibneibs[ivneib]==cell) break;
#ifdef DEBUG
				if (ivneib==Nv) BUG("GRID INCONSISTENCY IN CONVECTION\n");
#endif
#else
				ivneib=(int)(mask[ivert]&cell->index)>>((uint)ivert<<1);
#endif
				neibvert=neibverts[ivneib];
				neibvarn=neibvert->var;
				neibvoln=neibvarn[ivoln];
				neibmasn=neibvarn[imasn];
				neibveln=neibvarn+iveln;
				neibdenn=neibmasn/neibvoln;
				prsf=0.5*(prsn+neibvarn[iprsn]);
#ifdef LOCALDIFFUSION
				neibcoords=neibvert->x;
				neibveldif=neibveln;
#endif
				//Face-velocity (central differencing)
				for (int i=0; i<DIM; i++)
				{
					velf[i]=0.5*(veln[i]+neibveln[i]);
///					velf[i]=velcf[i];///DDD checkerboard
#ifndef LOCALDIFFUSION
					for (int j=0; j<DIM; j++)
					{	int	k=i*DIM+j;
						//Face velocity gradient
						vgrf[k]=totvolci*(vgrc[k]+neibvgrc[k]);
					}
#endif
				}
				//Hybrid face velocity and density
				velmid=0.5*(LENGTH(velf));
				velup=SCLP(a,velf)*areai;
				if (velup>=velmid)
				{	denf=denn;
					for (int i=0; i<DIM; i++)
						velf[i]=veln[i];
				}
				else
				if (velup<=-velmid)
				{	denf=neibdenn;
					for (int i=0; i<DIM; i++)
						velf[i]=neibveln[i];
				}
				else
					denf=0.5*(denn+neibdenn);
			}//facetype==internal
			else//boundary face
			{
				switch (facetype)
				{
				case inlet:
				case presinlet:
					denf=inlet_density;///InletDensity(time.current);
///					if(SCLP(a,velcf)>=0)
///						for (int i=0; i<DIM; i++)
///							velf[i]=0.0;
///					else
						for (int i=0; i<DIM; i++)
							velf[i]=velcf[i];
#ifndef LOCALDIFFUSION
					//Face velocity gradient
					for (int i=0; i<ndim2; i++)
						vgrf[i]=0.0;
#endif
				break;
				case outlet:
				case presoutlet:
					denf=OutletDensity(runtime.current);
///					if (SCLP(a,veln)>0.0)
///					{	
						for (int i=0; i<DIM; i++)
							velf[i]=veln[i];
///					}
///					else
///					{	for (int i=0; i<DIM; i++)
///							velf[i]=0.0;
///					}
#ifndef LOCALDIFFUSION
					//Face velocity gradient
					for (int i=0; i<ndim2; i++)
						vgrf[i]=0.0;
#endif
				break;
				case dove:
					denf=denn;
					//Face velocity gradient
					for (int i=0; i<ndim2; i++)
						vgrf[i]=varn[i];
					//Face velocity
					for (int i=0; i<DIM; i++)
						velf[i]=velcf[i];
				break;
				default:
					denf=denn;
#ifndef LOCALDIFFUSION
					//Face velocity gradient
					for (int i=0; i<ndim2; i++)
						vgrf[i]=0.0;///vgrn[i];
#endif
					//Face velocity
					for (int i=0; i<DIM; i++)
						velf[i]=0.0;///velcf[i];
				}
				prsf=Pressure(denf);//*Rgas*Tamb;///DDD
#ifdef LOCALDIFFUSION
				FACECENTER(ivert,verts,xface);
				neibcoords=xface;
				neibveldif=velf;
#endif
			}
			dvol=-SCLP(a,velf);
			dm=denf*dvol;
			*admn+=dm;//Added Mass
#ifdef LOCALDIFFUSION
			//Derivative at the face
			{	double	dx[DIM],dist2,disti,//distance to the neib point
					dv[DIM],//function difference
					*x=vert->x;//coordinates of this node
				dist2=0.0;
				for(int i=0;i<DIM;i++)
				{	double	d=neibcoords[i]-x[i];
					dv[i]=neibveldif[i]-veln[i];
					dist2+=d*d;
				}
				disti=1./sqrt(dist2);
				for(int i=0;i<DIM;i++)
					dvdxface[i]=dv[i]*disti;
			}
#endif
			//Pressure drop
			dpf=(prsf-prsn);///node_pressure-neib_pressure;
			//Convection/Diffusion/Pressure
			for (int i=0; i<DIM; i++)
			{	double
#ifdef LOCALDIFFUSION
					ddif=visc*area*dvdxface[i],
#else
					*v=vgrf+i*DIM,
					ddif=visc*(SCLP(a,v)),
#endif
					dcon=dvol*velf[i],
					dmm=ddif+dcon;
				dmom[i]=dmm;//momentum entering the face
				dmmn[i]+=dmm;
				dprn[i]+=a[i]*dpf;
			}
			//Set boundary source terms
			if 
			(	facetype!=internal&&facetype!=dove
			)
			{	int wall=facetype==boundary?1:0;
#ifdef LIDDRIVENCAVITY
				double	wall_velocity[]={1.0, 0.0, 0.0};
				if(facetype==inlet)	wall=1;
#endif
///DDD				double	DDD area=LENGTH(a),areai;
				for (int jv=1; jv<Nv; jv++)
				{	int	jvert=(ivert+jv)%Nv;
					double	*var=verts[jvert]->var,
						*velv=var+iveln,
						vnorm=areai*(SCLP(velv,a)),
						*admv=var+iadmn,
						*dmmv=var+idmmn,
						*dprv=var+idprn;
					*admv+=dm;//Added Mass
					for (int i=0; i<DIM; i++)
					{//Wall-tang component of velocity:
						double
							shear=(double)wall*(area*velv[i]-vnorm*a[i]);
#ifdef LIDDRIVENCAVITY
						if(facetype==inlet)
							shear-=wall_velocity[i]-velv[i];
#endif
						dmmv[i]+=dmom[i];//Added Momentum
						dprv[i]+=a[i]*dpf+Cwall*shear;
					}
				}
			}
		}//next ivert: looping over cell-vertexes
#ifndef LOCALDIFFUSION
		for (int i=0; i<ndim2; i++)
			vgrc[i]=newvgrc[i];
#endif
		cell=cell->next;
	}	while (cell!=cell_root);
	//Update nodal variables
	node=node_root;
	do
	{	double
			*var=node->var,
			vol=var[ivoln],voli=1./vol,
			den,
			*damn=var+iadmn,
			dm=*damn,
			*mas=var+imasn,
			*vel=var+iveln,
			*dmm=var+idmmn,
			*prs=var+iprsn,
			*dpr=var+idprn,
			dtvoli=dt*voli,
			dtmasi;
		if
		(	
			//node->type==(int)inlet
			//||node->type==(int)presinlet
			//||node->type==(int)outlet
			node->type==(int)presoutlet
			||node->type<=(int)maxElementStatus
		)	goto endloop;
		if(node->type!=(int)presinlet)
		{	*mas+=dm*dt;
			den=*mas*voli;
		}
		else
		{	den=inlet_density;
			*mas=den*vol;
		}
		*damn=0.0;
		dtmasi=dt/(*mas);
		*prs=Pressure(den);//*Rgas*Tamb;
		for (int i=0; i<DIM; i++)
		{
		 	vel[i]+=
		 		dmm[i]*dtvoli
				-dpr[i]*dtmasi;
		}
		endloop: node=node->next;
	}	while(node!=node_root);
#ifdef NO_INLET_OUTFLOW
	//Forbid momentum outflow from the inlet
	if (inlet_cell_root!=NULL)
	{
		DCellList	*c=inlet_cell_root;
		do
		{	DCell	*cell=c->cell;
			DNode	**verts=cell->vert;
			ElementStatus	*facetypes=cell->facetype;
			double	*x0,*x1,
				d[Nv][DIM];//edges
			//Edges
			x0=verts[0]->x;
			for (int ivert=0; ivert<Nv; ivert++)
			{	double
					*x1=verts[(ivert+1)%Nv]->x;
				for (int i=0; i<DIM; i++)
					d[ivert][i]=x1[i]-x0[i];
				x0=x1;
			}
			for (int ivert=0; ivert<Nv; ivert++)
			{
				DNode	*vert=verts[ivert];
				if (facetypes[ivert]==presinlet)
				{	double	
						a[DIM],area,areai,
						*edge0=d[ivert],
						*edge1=d[(ivert+1)%Nv],
						*edge2=d[(ivert+2)%Nv];
					//Face area
					VECP(a,edge1,edge2); 
					if (SCLP(a,edge0)<0.0)
						for (int i=0; i<DIM; i++) a[i]*=-0.5;
					else
						for (int i=0; i<DIM; i++) a[i]*= 0.5;
					area=LENGTH(a);
					areai=1./area;
					for (int jv=0; jv<Nfv; jv++)
					{	int	jvert=(ivert+jv+1)%Nv;
						double	*var=verts[jvert]->var,
							*veln=var+iveln,
							vnorm=areai*(SCLP(veln,a));
						if (vnorm>0.0)
						for (int i=0; i<DIM; i++)
						{//Subtract outflow of momentum
							veln[i]-=vnorm*areai*a[i];
						}
					}
				}
			}			
			c=c->next;
		}	while(c!=inlet_cell_root);
	}
#endif
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

using namespace ModGasFlow;

int	getNvarModGasFlow()
{
	return maxVarModGasFlow;
}
void	defVarModGasFlow
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
#ifndef LOCALDIFFUSION
		case velgradc:
			strcpy(name,"CellVelocityGradient");
			type=cells;
			dim=DIM*DIM;
			break;
#endif
		case dmomentum:
			strcpy(name,"AddedMomentum");
			type=nodes;
			dim=DIM;
			break;
		case pressure:
			strcpy(name,"Pressure");
			type=nodes;
			dim=1;
			break;
		case dpressure:
			strcpy(name,"PressureDrop");
			type=nodes;
			dim=DIM;
			break;
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
void	initVarModGasFlow
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
void	initModGasFlow(Domain *dom)
{
	Variable	*variable=dom->variable;
	int
		ivoln=variable[nvolume].loc,
		ivolc=variable[cvolume].loc,
		imasn=variable[mass].loc,
		iadmn=variable[addmass].loc,
		iveln=variable[velocity].loc,
		idmmn=variable[dmomentum].loc,
		idprn=variable[dpressure].loc,
		iprsn=variable[pressure].loc;
	DNode	*node_root=dom->dnode_root,
		*node=node_root;
	DCell	*cell_root=dom->dcell_root,
		*cell=cell_root;
	dom->setCellVolumes(cvolume);
	dom->setNodeVolumes(nvolume,cvolume);
	dom->setBoundaryVertexes();
	dom->indexNeibNodes();
	dom->setBoundaryCellFlags();
	if (node_root==NULL) return;
	do//Initialize
	{	double
			*var=node->var,
			vol=var[ivoln],voli=1./vol,
			*mas=var+imasn,
			*dmn=var+iadmn,
			*vel=var+iveln,
			*dmm=var+idmmn,
			*prs=var+iprsn,
			*dpr=var+idprn,
			den=OutletDensity(0.0);//density
		if(node->type==inlet||node->type==presinlet)
			den=InletDensity(0.0);
		*mas=den*vol;
		*prs=Pressure(den);
		*dmn=0.0;
		for (int i=0; i<DIM; i++)
		{	int	ii=DIM*i;
			vel[i]=0.0;
			dmm[i]=0.0;
			dpr[i]=0.0;
		}
//-			//LIDDRIVENCAVITY_TEST:
//-			if(node->type==inlet)
//-			{	double	*x=node->x;
//-				if(x[1]<0)node->type=(int)boundary;
//-				else
//-					for (int i=0; i<DIM; i++)
//-						vel[0]=1.0;
//-			}
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
void	stepModGasFlow(double dt, Domain *dom)
{
	static const	double	dtinj=1.0;//injection time interval
//	if ((int)(time.current-time.step)/dtinj!=(int)(time.current/dtinj))
//printf("time.current=%g, np=%d\n",time.current,dom->getNoPoints());FLUSH;///DDD
///	if (time.current>20.0)
///		InjectParticle(dom);
	AdvanceFlow
	(	InletDensity(runtime.current),
		OutletDensity(runtime.current),
		dt,dom
	);
///	AdvanceParticles(dt,dom);
	runtime.current+=dt;
}
