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

namespace ModMembrane
{
	enum VariablesModMembrane
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
		bforce,//boundary forces
		xold,//old coordinates of vertexes
		pmass,//particle mass
		pvelocity,//particle velocity
		height,//boundary cell heights 
		base,//boundary cell lateral positions
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
		maxVarModMembrane
	};
	const int	Nxfloor=16,Nyfloor=16,Nneedle=0,///2,
		ifloor=1;
	int
		Nparticles,
		Ntoolbnodes,
		ivoln,
		ivolc,
		imasn,
		iadmn,
		iveln,
		idmmn,
		idprn,
		iprsn,
		iheig,
		ibase,
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
///		Damb=Pamb/(Rgas*Tamb),//Ambient tensity
		Cwall=0.5,//0.75;//Wall-friction coefficient [mass/(length^2*time)]
		visc=1.50e0;///-8//Viscosity, [length^2/time]
	Point	*first_floor,*last_floor,
		*first_needle,*last_needle,
		*first_particle,*last_particle;
	double	dneedle[DIM];
double	InletDensity(double time)
{
///	static double
///		time_scale=50.0,
///		frequency=1./time_scale,
///		den0=Pamb/(Rgas*Tamb),
///		den1=Pin/(Rgas*Tin);
///	if(time<time_scale)return den1;
///	return den0;
///	return den0+dd*(1.0+sin(time*frequency));
///		if(time<200.0)den=Pin/(Rgas*Tin);
///		else	den=Pamb/(Rgas*Tamb);
///		return	den;
	return Pin/(Rgas*Tin);
}
double	OutletDensity(double time)
{
	return	Pamb/(Rgas*Tamb);
}
double	Pressure(double den)
{
	return den*Rgas*Tamb;
}
///	void	insideCircle
///	(	double *X,//Coordinates of the vertexes: X[0:Nfv*DIM-1]
///		double *E,//Edge-vectors: D[0:Nfv*DIM-1]
///		// RETURNS:
///		double *R,//Coordinates of the circle: R[DIM]
///		double *rr //Circle radius
///	)
///	{	double
///			B[Nfv*DIM];//Bisector unit vectors
///		for(int iv=0;iv<Nfv;iv++)
///		{	int	imv=DIM*iv,
///				iv1=(iv+1)%Nfv,
///				imv1=DIM*iv1;
///			double
///				*b=B+imv,
///				*e=E+imv,
///				*e1=E+imv1,
///				di,//inverse of d
///				d=0.0;
///			for(int i=0;i<DIM;i++)
///			{	double	r=e1[i]-e[i];
///				b[i]=r;
///				d+=r*r;
///			}
///			di=1./sqrt(d);
///			for(int i=0;i<DIM;i++)b[i]*=di;
///		}
///		{//Compute distance from node iv=0 to the center of the circle
///			int	iv=0,
///				iv1=iv+1, //%Nfv
///				iv2=iv+2,//%Nfv
///				imv=DIM*iv,
///				imv1=DIM*iv1,
///				imv2=DIM*iv2;
///			double
///				*b=B+imv,
///				*b1=B+imv1,
///				*b2=B+imv2,
///				*e=E+imv,
///				*e1=E+imv1,
///				*e2=E+imv2,
///				bb1=SCLP(b,b1),
///				bb2=SCLP(b,b2),
///				b1b2=SCLP(b1,b2),
///				eb2=SCLP(e,b2),
///				e1b=SCLP(e1,b),
///				e2b=SCLP(e2,b),
///				bb=(e1b+e2b*bb1+eb2*b1b2*bb1)/(1.0-bb1*bb2*b1b2);
///			for(int i=0;i<DIM;i++)
///				R[i]=X[imv+i]+bb*b[i];
///			*rr=bb*SCLP(b,e1)/LENGTH(e1);
///		}
///	}
void	trihi
(
	double	*edge,//edges of the triangle [Nfv*DIM]
	double	*base,//[Nfv]
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
		base[i0]=d*ei;
		for(int i=0;i<DIM;i++)
			h[i]=e0[i]+d*e[i];
		heig[i0]=LENGTH(h);
	}
}
void	trihi0
(
	double	edge[],//edges of the triangle
	double	heig[],
	double	base[]
)
{
	for (int i0=0; i0<Nfv; i0++)
	{	int
			j0=DIM*i0,
			i1=(i0+1)%Nfv,j1=DIM*i1;
		double	d,
			e[DIM],//unit vector in the direction of b
			*b=base+j0,
			*h=heig+j0,
			*e0=edge+j0,
			*e1=edge+j1,
			ei=1./(LENGTH(e1));
		for(int i=0;i<DIM;i++)e[i]=ei*e1[i];
		d=-(SCLP(e0,e));
		for(int i=0;i<DIM;i++)
		{	double	r=d*e[i];
			b[i]=r;
			h[i]=e0[i]+r;
		}
	}
/* 
b0+b1==b
a^2-b0^2==c^2-b1^2
a^2-b0^2==c^2-(b-b0)^2=c^2-b^2+2*b*b0-b0^2
a^2==c^2-b^2+2*b*b0
base=b0=(a^2-c^2+b^2)/(2*b)
height=sqrt(a^2-base^2)
*/
///		for (int i=0; i<Nfv; i++)
///		{	int
///				i1=(i+1)%Nfv,
///				i2=(i+2)%Nfv;
///			double
///				a=d[i],a2=a*a,
///				b=d[i1],
///				c=d[i2],
///				bb=(a2+b*b-c*c)/(2.*b);
///			heig[i]=sqrt(a2-bb*bb);
///			base[i]=bb;
///		}
}
void	init(Domain *dom)
{
	Tool	*&tool=dom->tool;
	Variable	*variable=dom->variable;
	ivoln=variable[nvolume].loc;
	ivolc=variable[cvolume].loc;
	imasn=variable[mass].loc;
	iadmn=variable[addmass].loc;
	iveln=variable[velocity].loc;
	idmmn=variable[dmomentum].loc;
	idprn=variable[dpressure].loc;
	iprsn=variable[pressure].loc;
	iheig=variable[height].loc;
	ibase=variable[base].loc;
	ixold=variable[xold].loc;
	ibfor=variable[bforce].loc;
	DNode	*node_root=dom->dnode_root,
		*node=node_root;
	DCell	*cell_root=dom->dcell_root,
		*cell=cell_root;
	ModGasFlow::Rgas=ModMembrane::Rgas;
	ModGasFlow::Tamb=ModMembrane::Tamb;//Gas-constant of air [work/(mass*temp)
	ModGasFlow::visc=ModMembrane::visc;///-8//Viscosity, [length^2/time]
	ModGasFlow::Pamb=ModMembrane::Pamb;//Ambient pressure [mass/(length*time^2)]
	ModGasFlow::Pin =ModMembrane::Pin;//Inlet pressure
	ModGasFlow::Tin =ModMembrane::Tin;
	ModGasFlow::Cwall=ModMembrane::Cwall;//Wall-friction coefficient [mass/(length^2*time)]
	dom->setCellVolumes(cvolume);
	dom->setNodeVolumes(nvolume,cvolume);
	if(option.debug){printf("indexNeibNodes\n");fflush(stdout);}///DDD
	dom->indexNeibNodes();
	if(option.debug){printf("setBoundaryCellFlags\n");fflush(stdout);}///DDD
	dom->setBoundaryCellFlags();
	if(option.debug){printf("setBoundaryVertexes\n");fflush(stdout);}///DDD
	dom->setBoundaryVertexes();
	if(option.debug){printf("connectBoundaryCellsToFaces\n");fflush(stdout);}///DDD
	dom->connectBoundaryCellsToFaces(dom->bface_root);
	if(option.debug){printf("createBoundaryNodeList\n");fflush(stdout);}///DDD
	dom->createBoundaryNodeList(dom->bnode_root);
	if(option.debug){printf("createBoundaryCellList\n");fflush(stdout);}///DDD
	dom->createBoundaryCellList(dom->bcell_root);
	if(option.debug){printf("Relax the grid to put the nodes to the centroids\n");fflush(stdout);}///DDD
	//Relax the grid to put the nodes to the centroids: IMPORTANT!
	for(int i=0;i<100;i++)
	dom->relax//Relax internal nodes
	(	ixold,
	 	ivoln,
 		ivolc,
		0.0 //relaxation factor
	);
	//Initialize face variables
	if(option.debug){printf("Initialize face variables\n");fflush(stdout);}///DDD
	if(dom->bface_root!=NULL)
	{	BFaceList
			*root=dom->bface_root,
			*face=root;
		do
		{	int	ivert=face->iface;
			double *y,
			*varf=face->var,
			*heig=varf+iheig,
			*base=varf+ibase,
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
#include "particles.cc"
void	AdvanceParticles
( Point	*first,
	Point	*last,
	double dt, 
	Domain *dom
)
{ int	iveloc=dom->variable[velocity].loc;
	DCell	*root_cell=dom->dcell_root;
	Point	*root=dom->origin,
		*p=first,*lastnext=last->next;
	double
		*M=dom->variable[pmass].val,
		*X=dom->coordinates[points].val,
		*V=dom->variable[pvelocity].val;
	do
	{	int	ip=p-root,imp=DIM*ip;
		ElementStatus	bc;//boundary conditions
		double
			*m=M+ip,//particle mass
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
			*m=-1.0;
			break;
			default:
			dom->delp(p);
			break;
		}
		next_loop: p=p->next;
	}	while(p!=lastnext);
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
				prsn=varn[iprsn],//pressure inside the cell
				dp=prsn-prsf,//pressure drop accross the boundary
				dpf=dp*area/(double)Nfv,
				//masn=varn[imasn],//density inside the cell
				//voln=varn[ivoln],volni=1./voln,
				//denn=masn*volni,//density inside the cell
				*norm=face->norm,
				*varf=face->var,//face-varialbes
				*heig0=varf+iheig,
				*base0=varf+ibase,
				heig[Nfv],base[Nfv],
///				r,//radius of the inside circle
///				R[DIM],//center of mass of the triangle
///				X[Nfv*DIM],//edge coordinates
///				e[Nfv*DIM],//edge unit direction vectors
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
//			insideCircle(X,E,R,&r);
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
///				masf=var[imasn],//mass of fluid at the boundary
///				masidt=dt/masf,
///				velfold[DIM],
///				*velf=var+iveln,//velocity of the boundary node
///				bmas=var[ibmas],
///				*bvel=var+ibvel,
				*bfor=var+ibfor,
				force=LENGTH(bfor);
			if(node->state.insidetool|node->state.toolsurface|node->state.fixed)goto loop;
//				for(int i=0;i<DIM;i++)velfold[i]=velf[i];
//				for(int i=0;i<DIM;i++)
//					velf[i]=(velf[i]+masidt*bfor[i])*flowdrag;
//				if(force*dt<=bigstep)
//				for(int i=0;i<DIM;i++)
//					x[i]+=Cwall*0.5*(velfold[i]+velf[i])*dt;
//				else
//				{	double	factor=bigstep/force;
//					for(int i=0;i<DIM;i++)
//						x[i]+=factor*bfor[i];//*dt;
//				}
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

using namespace ModMembrane;
using namespace ModGasFlow;

int	getNvarModMembrane()
{
	return maxVarModMembrane;
}
void	defVarModMembrane
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
		case xold:
			strcpy(name,"OldNodeCoordinates");
			type=nodes;
			dim=DIM;
			break;
		case bforce:
			strcpy(name,"BoundaryForce");
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
		case height:
			strcpy(name,"Height");
			type=boundary_faces;
			dim=DIM;
			break;
		case base:
			strcpy(name,"Base");
			type=boundary_faces;
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
void	initVarModMembrane
(	int	ivar,
	char	*name,
	int	type,
	int	loc,
	int	dimension,
	double	*val,
	Domain	*dom
)
{	switch(type)
	{	case nodes:
			setZero(loc,dimension,dom->dnode_root);
			break;
		case cells:
			setZero(loc,dimension,dom->dcell_root);
			break;
		case boundary_faces:
			setZero(loc,dimension,dom->bface_root);
			break;
		case points:
		{	int	mp=dom->getMaxNoPoints();
			for (int i=0; i<mp*dimension; i++)
				val[i]=0.0;
		}
		break;
		default:
			fprintf(stderr,"Can't initialize elements of type %d\n",type);
			exit(1);
	}
}
void	initMembrane(Domain *dom)
{
	Tool	*&tool=dom->tool;
	DNode	*node_root=dom->dnode_root,
		*node=node_root;
	DCell	*cell_root=dom->dcell_root,
		*cell=cell_root;
	ModMembrane::init(dom);
	//Initialize points
	dom->initp();
	if(node_root==NULL) return;
	do//Initialize nodal variables
	{	double
			*var=node->var,
			vol=var[ivoln],voli=1./vol,
			*mas=var+imasn,
			*dmn=var+iadmn,
			*vel=var+iveln,
			*dmm=var+idmmn,
			*prs=var+iprsn,
			*dpr=var+idprn,
			*bfr=var+ibfor,
			den=ModMembrane::InletDensity(0.0);
		*mas=den*vol;
		*prs=ModGasFlow::Pressure(den,ModMembrane::Tin);
		*dmn=0.0;
		for (int i=0; i<DIM; i++)
		{	int	ii=DIM*i;
			vel[i]=0.0;
			dmm[i]=0.0;
			dpr[i]=0.0;
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
	ModGasFlow::visc=1.5e1;
	Cshear=1.5;// normal stress
	Cstrain=1.5;// shear stress
}
void	stepMembrane(double dt, Domain *dom)
{	const	double	dtool=.5;
	double
		inlet_density=ModGasFlow::InletDensity(runtime.current),
		outlet_density=ModGasFlow::OutletDensity(runtime.current);
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
///	AdvanceFlow(inlet_density,outlet_density,dt,dom);
	MoveBoundary(dt,dom);
	dom->relax//Relax internal nodes
	(	ixold,
		ivoln,
		ivolc,
		0.5 //relaxation factor
	);
	runtime.current+=dt;
}
#define FLOOR_LEVEL	-3.0
void	initMicrobe(Domain *dom)
{
	Tool	*&tool=dom->tool;
	DNode	*node_root=dom->dnode_root,
		*node=node_root;
	DCell	*cell_root=dom->dcell_root,
		*cell=cell_root;
	ModMembrane::init(dom);
	//Initialize points
	dom->initp();
	if(node_root==NULL) return;
	do//Initialize nodal variables
	{	double
			*var=node->var,
			vol=var[ivoln],voli=1./vol,
			*mas=var+imasn,
			*dmn=var+iadmn,
			*vel=var+iveln,
			*dmm=var+idmmn,
			*prs=var+iprsn,
			*dpr=var+idprn,
			*bfr=var+ibfor,
			den=ModMembrane::InletDensity(0.0);
		*mas=den*vol;
		*prs=ModGasFlow::Pressure(den,ModMembrane::Tin);
		*dmn=0.0;
		for (int i=0; i<DIM; i++)
		{	int	ii=DIM*i;
			vel[i]=0.0;
			dmm[i]=0.0;
			dpr[i]=0.0;
			bfr[i]=0.0;
		}
		node=node->next;
	}	while(node!=node_root);
	//Create the floor
	if(dom->getMaxNoPoints()-dom->getNoPoints()<Nxfloor*Nyfloor) 
		ERROR("Not enough points to build the floor\n");
	first_floor=dom->origin;
	last_floor=first_floor;
	{	int	Nfl[2];
		double	fmin[2],fmax[2],df[2],
					xinj[DIM],dx,//injection point and dispersion
					minj=1.0,//injection mass
					vinj[DIM];//injection velocity
		Point	*origin=dom->origin;
		for (int i=0; i<2; i++)
		{	Nfl[i]=16;
			fmin[i]=-5.0;
			fmax[i]= 5.0;
			df[i]=(fmax[i]-fmin[i])/(double)(Nfl[i]-1);
		}
		minj=5.0;
		xinj[ifloor]=FLOOR_LEVEL;
		for (int i=0; i<DIM; i++) vinj[i]=0.0;
		for (int ix=0; ix<Nxfloor; ix++)
		{	xinj[(ifloor+DIM-1)%DIM]=fmin[0]+(double)ix*df[0];
			for (int iy=0; iy<Nyfloor; iy++)
			{	int	ip;
				xinj[(ifloor+1)%DIM]=fmin[1]+(double)iy*df[1];
				ip=dom->putp(xinj);
				dom->setsclp(ip,pmass,minj);
				dom->setvecp(ip,pvelocity,vinj);
				last_floor=origin+ip;
				last_floor->cell=cell;
			}
		}
	}
	//Create particles
	if(dom->getMaxNoPoints()-dom->getNoPoints()<1) 
		ERROR("Not enough particles\n");
	first_particle=last_floor->next;
	last_particle=first_particle;
	//Create the needle
	if(dom->getMaxNoPoints()-dom->getNoPoints()<Nneedle) 
	{	fprintf
		(	stderr,
			"No. of particles %d is not enough to create the needle\n",
			last_floor-dom->origin
		);exit(1);
	}
	first_needle=last_particle->next;
	last_needle=first_needle;
	if(Nneedle>0)
	{
		double	x0[DIM],x1[DIM],
					xinj[DIM],//injection point and dispersion
					minj=1.0,//injection mass
					vinj[DIM];//injection velocity
		Point	*origin=dom->origin;
		x0[0]=4.0;x0[1]=4.0;x0[2]=4.0;
		x1[0]=4.05;x1[1]=4.05;x1[2]=4.05;
		for(int i=0; i<DIM; i++)
		{	dneedle[i]=(x0[i]-x1[i])/(double)(Nneedle-1);
			vinj[i]=0.0;
		}
		minj=25.0;
		for (int in=0; in<Nneedle; in++)
		{	int	ip;
			for(int i=0; i<DIM; i++)
				xinj[i]=x0[i]+(double)in*dneedle[i];
			ip=dom->putp(xinj);
			dom->setsclp(ip,pmass,minj);
			dom->setvecp(ip,pvelocity,vinj);
			last_needle=origin+ip;
			last_needle->cell=cell;
		}
	}
	if (tool!=NULL)
	{//Position the tool
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
	ModGasFlow::visc=1.5e1;
	Cshear=1.5;// normal stress
	Cstrain=1.5;// shear stress
}
void	stepMicrobe(double dt, Domain *dom)
{	const	double	dtool=.5;
	double
		inlet_density=ModGasFlow::InletDensity(runtime.current),
		outlet_density=ModGasFlow::OutletDensity(runtime.current);
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
///		tool->relax(dom,0.3);
		tool->push();
		if (current_framenode!=last_framenode)
			current_framenode=current_framenode->next;
		tool->setype(current_framenode->type);
		tool->pos(current_framenode->x);
		tool->setradius(current_framenode->radius);
///		tool->setrefinement(current_framenode->refinement);
		tool->cleanSurface();
		Ntoolbnodes=tool->createLists(dom);
///			Ntoolbnodes=tool->updateLists(dom);
//		tool->setboundary();
	}
	AdvanceFlow(inlet_density,outlet_density,dt,dom);
	if (dom->dnode_root!=NULL)
	{	DNode	*node_root=dom->dnode_root,*node=node_root;
	do//Fix the floor
	{	if(node->x[ifloor]<=FLOOR_LEVEL)node->state.fixed=1;
		else node->state.fixed=0;
		node=node->next;
	}	while(node!=node_root);
	}
	MoveBoundary(dt,dom);
	dom->relax//Relax internal nodes
	(	ixold,
		ivoln,
		ivolc,
		0.5 //relaxation factor
	);
//	AdvanceParticles(first_particle,last_particle,dt,dom);
	runtime.current+=dt;
}
#undef FLOOR_LEVEL
//Case Breath:
void	initBreath(Domain *dom)
{
	Tool	*&tool=dom->tool;
	DNode	*node_root=dom->dnode_root,
		*node=node_root;
	DCell	*cell_root=dom->dcell_root,
		*cell=cell_root;
	ModMembrane::init(dom);
	if (node_root==NULL) return;
	do//Initialize Nodal variables
	{	double
			*var=node->var,
			vol=var[ivoln],voli=1./vol,
			*mas=var+imasn,
			*dmn=var+iadmn,
			*vel=var+iveln,
			*dmm=var+idmmn,
			*prs=var+iprsn,
			*dpr=var+idprn,
			den=ModMembrane::OutletDensity(0.0);//density
		*mas=den*vol;
		*prs=ModMembrane::Pressure(den);
		*dmn=0.0;
		for (int i=0; i<DIM; i++)
		{	int	ii=DIM*i;
			vel[i]=0.0;
			dmm[i]=0.0;
			dpr[i]=0.0;
		}
		node=node->next;
	}	while(node!=node_root);
	//Initialize points
	dom->initp();
	Cshear=10.0;
	Cstrain=10.0;
	ModGasFlow::Cwall=1.;
}
void	stepBreath(double dt, Domain *dom)
{
	ModGasFlow::AdvanceFlow
	(	ModMembrane::InletDensity(runtime.current),
		ModMembrane::OutletDensity(runtime.current),
		dt,dom
	);
	MoveBoundary(dt,dom);
	dom->relax//Relax internal nodes
	(	ixold,
		ivoln,
		ivolc,
		0.5 //relaxation factor
	);
 	runtime.current+=dt;
}
