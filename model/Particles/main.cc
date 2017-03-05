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
#include "particles.h"

extern void getvec(int loc, DCell *c, double *x, double *v);

namespace ModParticles
{
	enum VariablesModParticles
	{
		pmass=0,//particle mass
		pvelocity,//particle velocity
		maxVarModParticles
	};
	DCellList	*inlet_cell_root=NULL,
		*inlet_cell;

void	InjectParticle(Domain *dom)
{
	int	ip;//index to the new particle
	double
		xc[DIM],//cell-center coordinates
		xinj[DIM],dx,//injection point and dispersion
		minj=.1,//injection mass
		vinj[DIM];//injection velocity
		//Create injection coordinates array
	Point	*origin=dom->origin,
		*p;//current particle
	DCell	*cell;
	DNode	**verts;
	if (inlet_cell==NULL) return;
	cell=inlet_cell->cell;
	verts=cell->vert;
	for (int i=0; i<DIM; i++)
	{	double	xcenter=0.0;
		for (int iv=0; iv<Nv; iv++)
			xcenter+=verts[iv]->x[i];
		xcenter/=(double)Nv;
		xc[i]=xcenter;
	}
	dx=0.0;
	for (int i=0; i<DIM; i++)
	{	double	d=0.0;
		for (int iv=0; iv<Nv; iv++)
		{	d+=fabs(verts[iv]->x[i]-xc[i]);
		}
		d/=(double)Nv;
		dx+=d*d;
	}
	dx=0.5*sqrt(dx);
	do
	{	for (int i=0; i<DIM; i++)
		{	double	r=2.0*RND-1.0;
			xinj[i]=xc[i]+r*dx;
		}
	}	while(dom->ihostcell(xinj,cell)>=0);
	inlet_cell=inlet_cell->next;
	ip=dom->putp(xinj);
	dom->setsclp(ip,pmass,minj);
	//Set injection velocity
	vinj[0]=1.0;vinj[1]=1.0;vinj[2]=1.0;
///	getvec(iveloc,cell,xinj,vinj);
	dom->setvecp(ip,pvelocity,vinj);
	p=origin+ip;
	p->cell=cell;
}
void	InjectParticle
(	Domain *dom, 
	double spread // usually set to 0.5
)
{
	int	ip;//index to the new particle
	double
		xinj[DIM],dx,//injection point and dispersion
		minj=.1,//injection mass
		vinj[DIM];//injection velocity
		//Create injection coordinates array
	Point	*origin=dom->origin,
		*p;//current particle
	DCell	*cell;
	DNode	**verts;
	if (inlet_cell==NULL) return;
	cell=inlet_cell->cell;
	verts=cell->vert;
	for (int i=0; i<DIM; i++)
	{	double	xcenter=0.0;
		for (int iv=0; iv<Nv; iv++)
			xcenter+=verts[iv]->x[i];
		xcenter/=(double)Nv;
		xinj[i]=xcenter;
	}
	dx=0.0;
	for (int i=0; i<DIM; i++)
	{	double	d=0.0;
		for (int iv=0; iv<Nv; iv++)
		{	d+=fabs(verts[iv]->x[i]-xinj[i]);
		}
		d/=(double)Nv;
		dx+=d*d;
	}
	dx=spread*sqrt(dx);
	for (int i=0; i<DIM; i++)
	{	double	r=2.0*RND-1.0;
		xinj[i]+=r*dx;
	}
	inlet_cell=inlet_cell->next;
	ip=dom->putp(xinj);
	dom->setsclp(ip,pmass,minj);
	//Set injection velocity
	vinj[0]=1.0;vinj[1]=1.0;vinj[2]=1.0;
//	dx=0.0;// directed from the center
//	for(int i=0;i<2;i++)
//		dx+=xinj[i]*xinj[i];
//	dx=sqrt(dx);
//	for(int i=0;i<2;i++)
//		vinj[i]=xinj[i]/dx;
//	vinj[2]=1.0;
///	getvec(iveloc,cell,xinj,vinj);
	dom->setvecp(ip,pvelocity,vinj);
	p=origin+ip;
	p->cell=cell;
}
void	InjectParticle
(	Domain *dom, 
	double *xinj,//injection coordinates 
	double *vinj //injection velocity
)
{
	int	ip;//index to the new particle
	double minj=1.0;///injection mass
		//Create injection coordinates array
	Point	*origin=dom->origin,
		*p;//current particle
	DCell	*cell;
	DNode	**verts;
	if((cell=dom->findcell(xinj))==NULL) return;
	verts=cell->vert;
	ip=dom->putp(xinj);
	dom->setsclp(ip,pmass,minj);
	//Set injection velocity
///	getvec(iveloc,cell,xinj,vinj);
	dom->setvecp(ip,pvelocity,vinj);
	p=origin+ip;
	p->cell=cell;
}
void	InjectParticle
(	Domain *dom, 
	double *xinj,//injection coordinates 
	double *vinj,//injection velocity
	double minj //particle mass
)
{
	int	ip;//index to the new particle
		//Create injection coordinates array
	Point	*origin=dom->origin,
		*p;//current particle
	DCell	*cell;
	DNode	**verts;
	if((cell=dom->findcell(xinj))==NULL) return;
	verts=cell->vert;
	ip=dom->putp(xinj);
	dom->setsclp(ip,pmass,minj);
	//Set injection velocity
///	getvec(iveloc,cell,xinj,vinj);
	dom->setvecp(ip,pvelocity,vinj);
	p=origin+ip;
	p->cell=cell;
	printf("");
	if(option.debug||option.verbose)
	{	printf
		(	"Particle injected: ip=%d; x=%g,%g,%g; v=%g,%g,%g; m=%g\n",
			ip,xinj[0],xinj[1],xinj[2],vinj[0],vinj[1],vinj[2],minj
		);FLUSH;
	}
}
void	MoveParticle
(
	double	dt,//time step
	double	pmass,//particle mass
	double	*fvel,//flow velocity at particle location
	double	*x,//particle coordinates
	double	*v //particle velocity
)
{
	static const double Cdrag=.1;
	for (int i=0; i<DIM; i++)
	{	double
			drag=Cdrag*(fvel[i]-v[i]),
			vnew=v[i]+drag/pmass*dt;
		x[i]+=0.5*(vnew+v[i])*dt;
		v[i]=vnew;
	}
}
void	Interaction
(
	Domain	*dom,
 	Point	*origin,
 	Point	*first,
 	Point	*last,
	int	ip,   //index to the current particle
	double	dt //time step
)
{	if(dom->getNoPoints()==0)return;
	const double 
	//	PI=4.0*atan(1.0),
		den=1.0;//particle density
	Point	
		*pi=origin+ip,// pointer to the current particle
		*pj=first;// pointer to the next particle 
	double
		*M=dom->variable[pmass].val,
		mi=M[ip],
		ri=den*pow(3.0/4.0/PI*mi/den,1.0/3.0),//particle radius
		*X=dom->coordinates[points].val,
		*V=dom->variable[pvelocity].val,
		*xi=X+ip*DIM,
		*vi=V+ip*DIM;
	do
	{	int	jp=pj-origin,imp=DIM*jp;
		if(ip!=jp)
		{
			double 
				*xj=X+jp*DIM,
				mj=M[jp],
				rj=den*pow(3.0/4.0/PI*mj/den,1.0/3.0),//particle radius
				d[DIM],d0,d2=0.0;
			for(int i=0;i<DIM;i++)
			{	double r=xj[i]-xi[i];
				d[i]=r;
				d2+=r*r;
			}
			d0=sqrt(d2);
			if(d0<ri+rj)
			{
				double 	m=mi+mj,
					*vj=V+jp*DIM,
					d1[DIM],//unit vector in direction d
					w[DIM],//CM velocity
					ui[DIM],uj[DIM],//particle velocities
							// in CM coordinates
					uni,unj;//lengths of ui and uj
				//Transfer to CM coordinates
				for(int k=0;k<DIM;k++)
				{	w[k]=(mi*vi[k]+mj*vj[k])/m;
					ui[k]=vi[k]-w[k];
					uj[k]=vj[k]-w[k];
					d1[k]=d[k]/d0;
				}
				uni=SCLP(ui,d1);
				unj=SCLP(uj,d1);
				if(uni-unj>0.0)
				for(int k=0;k<DIM;k++)
				{
					vi[k]=ui[k]-2.0*d1[k]*uni+w[k];
					vj[k]=uj[k]-2.0*d1[k]*unj+w[k];
				}
			}
		}
		pj=pj->next;
	}	while(pj!=last->next);
}
void	AdvanceParticles(double dt, Domain *dom)
{
	DCell	*root_cell=dom->dcell_root;
	
	Point	*origin=dom->origin,*&first=dom->first,*&last=dom->last,
		*p=first;
	
	double
		*M=dom->variable[pmass].val,
		*X=dom->coordinates[points].val,
		*V=dom->variable[pvelocity].val;
	if(dom->getNoPoints()==0)return;
	do
	{	int	ip=p-origin,imp=DIM*ip;
		ElementStatus	bc;//boundary conditions
		double
			*m=M+ip, //particle mass
			*x=X+imp,//particle position
			*v=V+imp,//particle velocity
			norm[DIM],//boundary-face outward normal vector
			xold[DIM],//old and middle particle position
			fvel[DIM];//flow velocity at particle location
		DCell	*pcell=(DCell*)p->cell;//cell containing the particle
		if (*m<0.0) goto next_loop;
		COPY3(fvel,v);//particles only
		COPY3(xold,x);//store the old particle position
//		getvec(iveloc,pcell,x,fvel);
		MoveParticle(dt,*m,fvel,x,v);
		Interaction(dom,origin,first,last,ip,dt);
#ifndef WITH_MPI	
 		if((bc=dom->findcell((double*)x,(DCell*)pcell,(DCell*&)p->cell,(double*)norm))!=internal)
#else
		/*replace the above line with the following 2 lines 
			because pgCC compiler won't pass it*/
		DCell *p_cell= (DCell*)(p->cell); 
		if((bc=dom->findcell(x,pcell,p_cell,norm))!=internal)
#endif
		/* ANALYSE BOUNDARY CONDITIONS */
		switch(bc)
		{	case noslip:
			*m=-1.0;//particle sticks to the wall
			break;
			case slip:
			case bounce:
			{	double
					reflect=bc==slip?1.0:2.0,// wall reflection parameter
					vel=-(SCLP(v,norm)); //wall-nomral velocity
				p->cell=pcell;
				// Reflect the particle trajectory
				for(int i=0;i<DIM;i++)
				{	x[i]=xold[i];
					v[i]+=reflect*vel*norm[i];
				}
			}
			break;
			default:
			dom->delp(p);
			break;
		}
		next_loop: p=p->next;
	}	while(p!=last->next);
}
void	getBoundaryCells
(	int	boundary,
///	Domain	*dom
	DCell	*cell_root,
	DCellList	*&boundary_cell_root
)
{	if (cell_root!=NULL)
	{	int	ninj=0;
		DCellList
			*boundary_cell=boundary_cell_root;
		DCell
///			*cell_root=dom->dcell_root,
			*cell=cell_root;
		if(boundary_cell_root==NULL)
		do
		{	ElementStatus	*facetype=cell->facetype;
			int iface;
			for (iface=0; iface<Nf; iface++)
				if(facetype[iface]!=internal&&facetype[iface]!=boundary)break;
			if(iface<Nf)goto next_loop;
			for (iface=0; iface<Nf; iface++)
			if (facetype[iface]==boundary)
			{	ninj++;
				if (boundary_cell_root==NULL)
				{	boundary_cell_root=new DCellList;
					boundary_cell=boundary_cell_root;
				}
				else
				{	boundary_cell->next=new DCellList;
					boundary_cell=boundary_cell->next;
				}
				boundary_cell->cell=cell;
				boundary_cell->next=boundary_cell_root;//make it a loop
				break;
			}
			next_loop: cell=cell->next;
		}	while(cell!=cell_root);
		if(option.verbose)printf("Number of injection points: %d\n",ninj);
	}
}
//?void	getBoundaryCellsStatic
//?(	int	boundary,
//?	int	nc,
//?	Cell	*cell_root,
//?	CellList	*&boundary_cell_root
//?)
//?{	if (nc>0)
//?	{	int	ninj=0;
//?		CellList	*boundary_cell=boundary_cell_root;
//?		if(boundary_cell_root==NULL)
//?		for(int ic=0;ic<nc;ic++)
//?		{///	ElementStatus	*facetype=cell->facetype;
//?			Cell	*cell=cell_root+ic;
//?			int *facetype=cell->cell,iface;
//?			for (iface=0; iface<Nf; iface++)
//?				if(facetype[iface]<0&&facetype[iface]!=boundary)break;
//?			if(iface<Nf)continue;
//?			for (iface=0; iface<Nf; iface++)
//?			if (facetype[iface]==boundary)
//?			{	ninj++;
//?				if (boundary_cell_root==NULL)
//?				{	boundary_cell_root=new CellList;
//?					boundary_cell=boundary_cell_root;
//?				}
//?				else
//?				{	boundary_cell->next=new CellList;
//?					boundary_cell=boundary_cell->next;
//?				}
//?				boundary_cell->cell=cell;
//?				boundary_cell->next=boundary_cell_root;//make it a loop
//?				break;
//?			}
//?		}
//?		if(option.verbose)printf("Number of injection points: %d\n",ninj);
//?	}
//?}
};//END NAMESPACE 

using namespace ModParticles;

int	getNvarModParticles()
{
	return maxVarModParticles;
}
void	defVarModParticles
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
void	initVarModParticles
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
void	initModParticles(Domain *dom)
{//Random particles injection
	DNode	*node_root=dom->dnode_root,
		*node=node_root;
	DCell	*cell_root=dom->dcell_root,
		*cell=cell_root;
///	dom->setCellVolumes(cvolume);
///	dom->setNodeVolumes(nvolume,cvolume);
	dom->setBoundaryVertexes();
	dom->indexNeibNodes();
	dom->setBoundaryCellFlags();
	if (node_root==NULL) return;
	dom->initp();//Initialize particles
	//Create the list of injection cells
	getBoundaryCells(presinlet,dom->dcell_root,inlet_cell_root);
	inlet_cell=inlet_cell_root;
#ifdef DEBUG
	dom->checkGrid(dom->dnode_root,dom->dcell_root);
#endif
	//Deterministic injection:
//	for (int i=0;i<dom->getMaxNoPoints();i++)
//		InjectParticle(dom,0.0);
}
void	stepModParticles(double dt, Domain *dom)
{//Random particles motion
//	static const	double	dtinj=1.0;//injection time interval
//	if ((int)(time.current-time.step)/dtinj!=(int)(time.current/dtinj))
	InjectParticle(dom);//stochastic injection every timestep
	AdvanceParticles(dt,dom);
	runtime.current+=dt;
}
void	initModParticles0(Domain *dom)
{//Two particle collision test
	DNode	*node_root=dom->dnode_root,
		*node=node_root;
	DCell	*cell_root=dom->dcell_root,
		*cell=cell_root;
///	dom->setCellVolumes(cvolume);
///	dom->setNodeVolumes(nvolume,cvolume);
	dom->setBoundaryVertexes();
	dom->indexNeibNodes();
	dom->setBoundaryCellFlags();
	if (node_root==NULL) return;
	dom->initp();//Initialize particles
	//Create the list of injection cells
	getBoundaryCells(presinlet,dom->dcell_root,inlet_cell_root);
	inlet_cell=inlet_cell_root;
#ifdef DEBUG
	dom->checkGrid(dom->dnode_root,dom->dcell_root);
#endif
	{	if (option.debug||option.verbose) printf("Injecting two particles\n");
		double x[DIM],v[DIM];//injection coordinates and velocicities
		x[0]=0.;x[1]=0.;x[2]=0.;
		v[0]=0.0;v[1]=0.0;v[2]=0.0;
		InjectParticle(dom,x,v);
		x[0]=-2.0;x[1]=0.0;x[2]=0.0;
		v[0]=2.0;v[1]=0.0;v[2]=0.0;
		InjectParticle(dom,x,v);
	}
}
void	stepModParticles0(double dt, Domain *dom)
{// two particle collision 
	AdvanceParticles(dt,dom);
	runtime.current+=dt;
}
