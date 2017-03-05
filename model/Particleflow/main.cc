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
#include "gasflow.h"
#include "templates.cc"

extern void getvec(int loc, DCell *c, double *x, double *v);

namespace ModParticleFlow
{
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
		///Rgas=0.5,//Gas-constant of air [work/(mass*temp)
		///Rgas=1.0,//Gas-constant of air [work/(mass*temp)
		Rgas=1.0,//Gas-constant of air [work/(mass*temp)
		Tamb=1.0,
		visc=1.50e-0,///-8//Viscosity, [length^2/time]
		Pamb=1.00e-0,//Ambient pressure [mass/(length*time^2)]
		Pin=2.0*Pamb,//Inlet pressure
		Tin=Tamb,
		///Cwall=0.1;//Wall-friction coefficient [mass/(length^2*time)]
		Cwall=0.1;//Wall-friction coefficient [mass/(length^2*time)]
	enum VariablesModParticleFlow
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
		maxVarModParticleFlow
	};
	DCellList	*inlet_cell_root=NULL,
		*inlet_cell;
void	InjectParticle(Domain *dom)
{
	int	ip,//index to the new particle
		iveloc=dom->variable[velocity].loc;
	double
		xinj[DIM],dx,//injection point and dispersion
		minj=1.0,//injection mass
		vinj[DIM];//injection velocity
		//Create injection coordinates array
	Point	*origin=dom->origin,
		*p;//current particle
	DCell	*cell;
	DNode	**verts;
	if (inlet_cell==NULL) return;
	//Randomize particle mass:
	{	double
			dm=0.9*minj, //deviation from the avarage
			r=2.0*RND-1.0;// random number between -1 and 1
		minj+=r*dm;
	}
	//Select cell-center:
	cell=inlet_cell->cell;
	verts=cell->vert;
	for (int i=0; i<DIM; i++)
	{	double	xcenter=0.0;
		for (int iv=0; iv<Nv; iv++)
			xcenter+=verts[iv]->x[i];
		xcenter/=(double)Nv;
		xinj[i]=xcenter;
	}
	//Randomize injection position:
	dx=0.0;
	for (int i=0; i<DIM; i++)
	{	double	d=0.0;
		for (int iv=0; iv<Nv; iv++)
		{	d+=fabs(verts[iv]->x[i]-xinj[i]);
		}
		d/=(double)Nv;
		dx+=d*d;
	}
	dx=0.5*sqrt(dx);
	for (int i=0; i<DIM; i++)
	{	double	r=2.0*RND-1.0;
		xinj[i]+=r*dx;
	}
	inlet_cell=inlet_cell->next;
	ip=dom->putp(xinj);
	dom->setsclp(ip,pmass,minj);
	getvec(iveloc,cell,xinj,vinj);
	dom->setvecp(ip,pvelocity,vinj);
	p=origin+ip;
	p->cell=cell;
}
void	MoveParticle
(
 	Point	*origin,
	int	ip,//number of the particle to move
	double	dt,//time step
	double	pmass,//particle mass
	double	*fvel,//flow velocity at particle location
	double	*x,//particle coordinates
	double	*v //particle velocity
)
{
	static const double Cdrag=.1;
	int	imp=DIM*ip;
	Point	*p=origin+ip;
	for (int i=0; i<DIM; i++)
	{	double
			drag=Cdrag*(fvel[i]-v[i]),
			vnew=v[i]+drag/pmass*dt;
		x[i]+=0.5*(vnew+v[i])*dt;
		v[i]=vnew;
	}
}
void	AdvanceParticles(double dt, Domain *dom)
{
	int	iveloc=dom->variable[velocity].loc;
	DCell	*root_cell=dom->dcell_root;
	Point	*root=dom->origin,*&first=dom->first,*&last=dom->last,
		*p=first;

	double
		*M=dom->variable[pmass].val,
		*X=dom->coordinates[points].val,
		*V=dom->variable[pvelocity].val;
	if(dom->getNoPoints()==0)return;
	do
	{	int	ip=p-root,imp=DIM*ip;
		ElementStatus	bc;//boundary conditions
		double
			*m=M+ip,
			*x=X+imp,//particle position
			*v=V+imp,//particle velocity
			fvel[DIM];//flow velocity at particle location
		DCell	*pcell=(DCell*)p->cell;//cell containing the particle
		if (*m<0.0) goto next_loop;
		ZERO3(fvel);
		getvec(iveloc,pcell,x,fvel);
		MoveParticle(root,ip,dt,*m,fvel,x,v);
#ifndef WITH_MPI
		if((bc=dom->findcell((double*)x,(DCell*)pcell,(DCell*&)p->cell))!=internal)
#else
		/*replace the above line with the following 2 lines 
			because pgCC compiler won't pass it*/
		DCell *p_cell= (DCell*)(p->cell); 
		if((bc=dom->findcell(x,pcell,p_cell))!=internal)
#endif
		switch(bc)
		{	case boundary:
			if(*m>0.0)*m=-*m;
			break;
			default:
			dom->delp(p);
			break;
		}
		next_loop: p=p->next;
	}	while(p!=last->next);
}
};//END NAMESPACE 

using namespace ModGasFlow;
using namespace ModParticleFlow;

int	getNvarModParticleFlow()
{
	return maxVarModParticleFlow;
}
void	defVarModParticleFlow
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
		case pmass:
			strcpy(name,"ParticleMass");
			type=points;
			dim=1;
			break;
		case pvelocity:
			strcpy(name,"ParticleVelocity");
			type=points;
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
void	initVarModParticleFlow
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
void	initModParticleFlow(Domain *dom)
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
	ModGasFlow::Rgas=ModParticleFlow::Rgas;
	ModGasFlow::Tamb=ModParticleFlow::Tamb;//Gas-constant of air [work/(mass*temp)
	ModGasFlow::visc=ModParticleFlow::visc;///-8//Viscosity, [length^2/time]
	ModGasFlow::Pamb=ModParticleFlow::Pamb;//Ambient pressure [mass/(length*time^2)]
	ModGasFlow::Pin =ModParticleFlow::Pin;//Inlet pressure
	ModGasFlow::Tin =ModParticleFlow::Tin;
	ModGasFlow::Cwall=ModParticleFlow::Cwall;//Wall-friction coefficient [mass/(length^2*time)]
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
void	stepModParticleFlow(double dt, Domain *dom)
{
	static const	double	dtinj=1.0;//injection time interval
//	if ((int)(time.current-time.step)/dtinj!=(int)(time.current/dtinj))
//printf("time.current=%g, np=%d\n",time.current,dom->getNoPoints());FLUSH;///DDD
	if (runtime.current>=0.0)
		InjectParticle(dom);
	AdvanceFlow
	(
		InletDensity(runtime.current),
		OutletDensity(runtime.current),
		dt,dom
	);
	AdvanceParticles(dt,dom);
	runtime.current+=dt;
}
