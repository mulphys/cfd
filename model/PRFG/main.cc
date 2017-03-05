
/* 
 *  MODEL: PRFG
 *  Particles in an RFG flow field
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
#include "mainrfg.h"
#include "particles.h"
#include "templates.cc"

//#define PRFG

extern "C"
{
	void	genspec_(int *N);
	void	genvel_(double *t, double *x, double *tt, double *tl, double *uu, double *vel);
};
extern void	getvec
(
 	int	vecloc,//offset of vector variable to be interpolated
	DCell	*cell,//cell where the variable is located
	double	*X,//point coordinates: X[0:DIM-1]
	//RETURNS:
	double	*V //vector at the point: V[0:DIM-1]
);

using namespace ModRFG;
using namespace ModParticles;

namespace ModPRFG
{
	const double 
		init_mass=1.0;
	enum VariablesModPRFG
	{
		pmass=0,
		velocity,
		pvelocity,//particle velocity
		maxVarModPRFG
	};
	int	iveln;
	DCellList	*inlet_cell_root=NULL,
		*inlet_cell;

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
	static const double Cdrag=1.0;
	int	imp=DIM*ip;
	Point	*p=origin+ip;
	for (int i=0; i<DIM; i++)
	{	double
//			drag=Cdrag*(fvel[i]-v[i]),
//			vnew=v[i]+drag/pmass*dt;
			vnew=fvel[i];
		x[i]+=0.5*(vnew+v[i])*dt;
		v[i]=vnew;
	}
}
void	InjectParticle(Domain *dom)
{
	int	ip;//index to the new particle
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
	dx=0.5*sqrt(dx);
	for (int i=0; i<DIM; i++)
	{	double	r=2.0*RND-1.0;
		xinj[i]+=r*dx;
	}
	inlet_cell=inlet_cell->next;
	ip=dom->putp(xinj);
	dom->setsclp(ip,pmass,minj);
	//Set injection velocity
	for(int i=0;i<DIM;i++)vinj[i]=0.0;//to zero
//	getvec(iveln,cell,xinj,vinj);//to flow velocity
	dom->setvecp(ip,pvelocity,vinj);
	p=origin+ip;
	p->cell=cell;
}
void	AdvanceParticles(double dt, Domain *dom)
{
	DCell	*root_cell=dom->dcell_root;
	Point	*root=dom->origin,*&first=dom->first,*&last=dom->last,
		*p=first;
	double
		*M=dom->variable[pmass].val,
		*X=dom->coordinates[points].val,
		*V=dom->variable[pvelocity].val;
	if(dom->getNoPoints()==0)return;
	do
	{	int	ip=p-root,imp=DIM*ip,
			ivert;//vertex opposite to the face crossed by the particle
		ElementStatus	bc;//boundary conditions
		double
			*m=M+ip,
			*x=X+imp,//particle position
			*v=V+imp,//particle velocity
			norm[DIM],//boundary-face outward normal vector
			xold[DIM],//old and middle particle position
			fvel[DIM];//flow velocity at particle location
		DCell	*pcell=(DCell*)p->cell;//cell containing the particle
		if (*m<0.0) goto next_loop;
		COPY3(xold,x);//store the old particle position
#ifdef PRFG
		genvel(runtime.current,x,1.0,1.0,1.0,1.0,1.0,1.0,1.0,fvel);
#else
		getvec(iveln,pcell,x,fvel);
#endif
		MoveParticle(root,ip,dt,*m,fvel,x,v);
#ifndef WITH_MPI
		if((bc=dom->findcell((double*)x,(DCell*)pcell,(DCell*&)p->cell,norm))!=internal)
#else
		/*replace the above line with the following 2 lines 
			because pgCC compiler won't pass it*/
		DCell *p_cell= (DCell*)(p->cell);
		if((bc=dom->findcell(x,pcell,p_cell,norm))!=internal)
#endif
		switch(bc)
		{	case noslip:
			*m=-1.0;//particle sticks to the wall
			break;
			case bounce:
			{	double	dist,vel;
				//Determine the face normal vector
				DNode	**verts=pcell->vert;
				p->cell=pcell;
				//Mirror-reflect the particle trajectory
				vel=-(SCLP(v,norm));
				for(int i=0;i<DIM;i++)
				{
					x[i]=xold[i];
					v[i]+=2.0*vel*norm[i];
				}
			}
			break;
///				{	double	d[3][DIM],norm[DIM],dmid[DIM],dist,vel;
///					//Determine the face normal vector
///					DNode	**verts=pcell->vert;
///					p->cell=pcell;
///					for(int j=0;j<3;j++)
///					for(int i=0;i<DIM;i++)
///						d[j][i]=verts[(ivert+j+1)%Nv]->x[i]-verts[(ivert+j)%Nv]->x[i];
///					VECP(norm,d[1],d[2]);
///					dist=LENGTH(norm);
///	#ifdef DEBUG
///					if(dist<SMALL)ERROR("AdvanceParticles: CELL-FACE UNDEFINED");
///	#endif
///					dist=1./dist;
///					if(SCLP(norm,d[0])>0)
///						MULVEC(norm,-dist);
///					else
///						MULVEC(norm, dist);
///					//now norm points inside the cell
///					//Mirror-reflect the particle trajectory
///	 				//Find the middle point
///					for(int i=0;i<DIM;i++)
///						dmid[i]=0.5*(x[i]-xold[i]);
///					dist=-(SCLP(dmid,norm));
///					vel=-(SCLP(v,norm));
///					for(int i=0;i<DIM;i++)
///					{	x[i]=xold[i];///+2.0*(dmid[i]+dist*norm[i]);
///						v[i]=v[i]+2.0*vel*norm[i];
///					}
///				}
			break;
			default:
			dom->delp(p);
			break;
		}
		next_loop: p=p->next;
	}	while(p!=last->next);
}
};

using namespace ModPRFG;

int	getNvarModPRFG()
{
	return maxVarModPRFG;
}
void	defVarModPRFG
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
		case velocity:
			strcpy(name,"Velocity");
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
void	initVarModPRFG
(	int	ivar,
	char	*name,
	int	type,
	int	loc,
	int	dim,
	double	*val,
	Domain	*dom
)
{	iveln=dom->variable[velocity].loc;
	switch(type)
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
void	initModPRFG(Domain *dom)
{
	int	ns=200;//RFG spectrum size
	if(option.verbose){printf("Initializing RFG with %d harmonics\n",ns);FLUSH;}
	genspec_(&ns);//Initilizing RFG spectrum
	if(option.verbose){printf("Initializing domain velocities with RFG\n");FLUSH;}
	rfg(iveln,dom);
	//See mod_softbody.cc
	//Initialize points
	if(option.verbose){printf("Initializing particles\n");FLUSH;}
	dom->initp();
	//Create particles
	if(dom->getMaxNoPoints()-dom->getNoPoints()<1) 
		ERROR("Not enough particles\n");
	ModParticles::getBoundaryCells(presinlet,dom->dcell_root,inlet_cell_root);
	inlet_cell=inlet_cell_root;
}
void	stepModPRFG(double dt, Domain *dom)
{	const	double	dtool=1.0;
	Tool	*tool=dom->tool;
	static int i=0;
	if(option.verbose)printf("Iteration no.=%d\n",i++);
	InjectParticle(dom);
	AdvanceParticles(dt,dom);
#ifndef PRFG
	rfg(iveln,dom);
#endif
	runtime.current+=dt;
}
