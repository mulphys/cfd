#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "main.h"
#include "io.h"
#include "vecalg.h"
#include "geom.h"
#include "func.h"

/* SOLVER FUNCTIONS */
//?void getBrickBC(char *p, char *s, char *filename){}
//?	void	exterVectorCellFaceCentral
//?	//Centeral-interpolates cell-centered vector to the faces
//?	(
//?		int	nbf,//number of boundary faces
//?		struct Cell	*cell,
//?		struct Face	*face,
//?		double	*U, //cell-center vector
//?		double	*W, //boundary-face vector
//?		double	*V  //extrapolated face vector
//?	)
//?	{	//Loop over boundary faces
//?		for (int iface=0; iface<nbf; iface++)
//?		{	int
//?				icell,//gloval neighbor cell index
//?				mface=DIM*iface;
//?			double	*u,*v,*w;//local variables
//?			struct Face	*f=face+iface;
//?			v=V+mface;
//?			w=W+mface;//boundary-vector
//?			icell=f->cell[0];
//?			u=U+DIM*icell;
//?			switch(f->cell[1])
//?			{
//?			case inlet:
//?				//Copy boundary values
//?				for (int i=0; i<DIM; i++)
//?					v[i]=w[i];
//?				break;
//?			case outlet:
//?				//Zero normal derivative
//?				for (int i=0; i<DIM; i++)
//?					v[i]=u[i];
//?				break;
//?			}
//?		}
//?	}
//?	void	interVectorCellFaceCentral
//?	//Centeral-interpolates cell-centered vector to the faces
//?	(
//?		int	nf, //number of faces
//?		int	nbf,//number of boundary faces
//?		struct Cell	*cell,
//?		struct Face	*face,
//?		double	*U, //cell-center vector
//?		double	*V  //face-interpolated variable
//?	)
//?	{	//Loop over internal faces
//?		for (int iface=nbf; iface<nf; iface++)
//?		{	struct Face	*f=face+iface;
//?			double	a,
//?				*v=V+DIM*iface,//local face interpolated variable
//?				*u0,*u1;//two neighbor cell-center vectors 
//?			int	*ncell=f->cell;//neighbor cell indexes
//?			u0=U+DIM*ncell[0];
//?			u1=U+DIM*ncell[1];
//?			for (int i=0;i<DIM;i++)
//?				v[i]=0.5*(u0[i]+u1[i]);
//?		}
//?	}
//?void	mapVecCell2Face
//?//Centeral-interpolates cell-centered vector to the faces
//?(
//?	int	nf, //number of faces
//?	int	nbf,//number of boundary faces
//?	struct Cell	*cell,
//?	struct Face	*face,
//?	double	*U, //cell-center vector
//?	double	*W, //boundary-face vector
//?	double	*V  //face-interpolated variable
//?)
//?{	//Loop over boundary faces
//?	for (int iface=0; iface<nbf; iface++)
//?	{	int
//?			icell,//gloval neighbor cell index
//?			mface=DIM*iface;
//?		double	*u,*v;//local variables
//?		struct Face	*f=face+iface;
//?		v=V+mface;
//?		icell=f->cell[0];
//?		u=U+DIM*icell;
//?		//	if (W==NULL)
//?		//	{
//?		//		for (int i=0; i<DIM; i++)
//?		//			v[i]=u[i];//Zero normal derivative
//?		//	}
//?		//	else
//?		switch(f->cell[1])
//?		{
//?		case inlet:
//?			//Copy boundary values
//?			if (W!=NULL)
//?				for (int i=0; i<DIM; i++)
//?					v[i]=W[mface+i];
//?			else
//?				for (int i=0; i<DIM; i++)
//?					v[i]=0.0;
//?					//v[i]=u[i];//Zero normal derivative
//?			break;
//?		case outlet:
//?			//Extrapolate the internal value to the boundary
//?			for (int i=0; i<DIM; i++)
//?				v[i]=u[i];
//?			break;
//?		default:
//?			if (W==NULL)
//?			for (int i=0; i<DIM; i++)
//?				v[i]=0.0;//Zero at all other boundaries
//?				//v[i]=u[i];//Extrapolate the internal value to the boundary
//?			else
//?			for (int i=0; i<DIM; i++)
//?				v[i]=W[mface+i];//Boundary value at all other boundaries
//?				//v[i]=u[i];//Zero normal derivative
//?			//fprintf(stderr,"\nnbf=%d, iface=%d\n",nbf,iface);
//?			//BUG("mapVectorCellFace: Non-boundary face in the wrong part of the array");
//?		}
//?	}
//?	//Loop over internal faces
//?	for (int iface=nbf; iface<nf; iface++)
//?	{	struct Face	*f=face+iface;
//?		double
//?			*v=V+DIM*iface,//local face interpolated variable
//?			*u0,*u1;//two neighbor cell-center vectors 
//?		int	*ncell=f->cell;//neighbor cell indexes
//?		u0=U+DIM*ncell[0];
//?		u1=U+DIM*ncell[1];
//?		for (int i=0;i<DIM;i++)
//?			v[i]=0.5*(u0[i]+u1[i]);
//?	}
//?}
//?void	mapVecCell2FaceUpwind
//?//Centeral-interpolates cell-centered vector to the faces
//?(
//?	int	nf, //number of faces
//?	int	nbf,//number of boundary faces
//?	struct Cell	*cell,
//?	struct Face	*face,
//?	double	*U, //cell-center vector
//?	double	*W, //boundary-face vector
//?	double	*V  //face-interpolated variable
//?)
//?{	//Loop over boundary faces
//?	for (int iface=0; iface<nbf; iface++)
//?	{	int
//?			icell,//gloval neighbor cell index
//?			mface=DIM*iface;
//?		double	*u,*v;//local variables
//?		struct Face	*f=face+iface;
//?		v=V+mface;
//?		icell=f->cell[0];
//?		u=U+DIM*icell;
//?		//	if (W==NULL)
//?		//	{
//?		//		for (int i=0; i<DIM; i++)
//?		//			v[i]=u[i];//Zero normal derivative
//?		//	}
//?		//	else
//?		switch(f->cell[1])
//?		{
//?		case inlet:
//?			//Copy boundary values
//?			if (W!=NULL)
//?				for (int i=0; i<DIM; i++)
//?					v[i]=W[mface+i];
//?			else
//?				for (int i=0; i<DIM; i++)
//?					v[i]=0.0;
//?					//v[i]=u[i];//Zero normal derivative
//?			break;
//?		case outlet:
//?			//Zero normal derivative
//?			for (int i=0; i<DIM; i++)
//?				v[i]=u[i];
//?			break;
//?		default:
//?			for (int i=0; i<DIM; i++)
//?				v[i]=0.0;//Zero value at all other boundaries
//?				//v[i]=u[i];//Zero normal derivative
//?			//fprintf(stderr,"\nnbf=%d, iface=%d\n",nbf,iface);
//?			//BUG("mapVectorCellFace: Non-boundary face in the wrong part of the array");
//?		}
//?	}
//?	//Loop over internal faces
//?	for (int iface=nbf; iface<nf; iface++)
//?	{	struct Face	*f=face+iface;
//?		double	a,
//?			*norm=f->norm,
//?			*v=V+DIM*iface,//local face interpolated variable
//?			*u0,*u1;//two neighbor cell-center vectors 
//?		int	*ncell=f->cell;//neighbor cell indexes
//?		u0=U+DIM*ncell[0];
//?		u1=U+DIM*ncell[1];
//?		for (int i=0;i<DIM;i++)
//?			v[i]=0.5*(u0[i]+u1[i]);
//?		a=SCLP(norm,v);
//?		if(a> SMALL) for (int i=0;i<DIM;i++)v[i]=u0[i];
//?		else
//?		if(a<-SMALL) for (int i=0;i<DIM;i++)v[i]=u1[i];
//?	}
//?}
//?void	conVecFace2Cell
//?(	//Convection of face variables to the cell-center
//?	int	nc,
//?	struct Cell	*cell,
//?	struct Face	*face,
//?	double	*Vol,//cell-volumes 
//?	double	*V,//convected face-variable
//?	double	*U,//face-velocity
//?	double	*C //cell-centered convection term
//?)
//?{	//static double	onethird=1./3.;
//?	for (int ic=0; ic<nc; ic++)
//?	{	struct Cell	*c=cell+ic;
//?		int
//?			//*cn=c->vert,//cn[i]=index of vertex i=0:Nv in the node array
//?			*cf=c->face,//cf[i]=index of face i in the global face array
//?			*cc=c->cell;//cc[i]=index of neighbor cell i in the global cell array
//?		double//	d=0.0;
//?			voli=1./Vol[ic],
//?			*con=C+DIM*ic;
//?		for (int k=0; k<DIM; k++)
//?			con[k]=0.0;
//?		for (int iface=0; iface<Nf; iface++)
//?		{	int
//?				jf=cf[iface],mf=DIM*jf,//global face index
//?				jc=cc[iface],//global neighbor-cell index
//?				direction=1-(((c->mirror.inside>>iface)%2)<<1);
//?			double d,
//?				area,//face areas
//?				*vel,//face-velocity
//?				*var,//face-variable
//?				*norm;//face normal vector (direction determined by direction)
//?			struct Face	*f;
//?			//if(jf>=0)direction=1;
//?			//else{direction=-1;jf=-jf-1;}
//?			f=face+jf;
//?			norm=f->norm;
//?			area=f->area;
//?			vel=U+mf;
//?			var=V+mf;
//?			d=(double)direction*voli*area*(SCLP(vel,norm)); //DEBUG: *onethird;
//?			for (int k=0; k<DIM; k++)
//?				con[k]+=d*var[k];
//?		}
//?	}
//?}
//?void	conCellVecUp
//?(	//Convection of cell variables to the cell-center 
//?	//with upwind interpolation 
//?	int	nc,
//?	struct Cell	*cell,
//?	struct Face	*face,
//?	double	*Vol,//cell-volumes 
//?	double	*V,//cell-centered variable
//?	double	*U,//cell-centered velocity
//?	double	*W,//boundary velocity
//?	double	*C //cell-centered convection term
//?)
//?{	static double	 zerovec[]={0.0,0.0,0.0};
//?	//static const double	onethird=1./3.;
//?	int mc=DIM*nc;
//?	for (int i=0; i<mc; i++) C[i]=0.0;
//?	for (int ic=0; ic<nc; ic++)
//?	{	struct Cell	*c=cell+ic;
//?		int
//?			*cn=c->vert,//cn[i]=index of vertex i=0:Nv in the node array
//?			*cf=c->face,//cf[i]=index of face i in the global face array
//?			*cc=c->cell;//cc[i]=index of neighbor cell i in the global cell array
//?		double//	d=0.0;
//?			voli=1./Vol[ic],
//?			*con=C+DIM*ic;
//?		for (int iface=0; iface<Nf; iface++)
//?		{	int
//?				jf=cf[iface],//global face index
//?				jc=cc[iface],//global neighbor-cell index
//?				direction=1-(((c->mirror.inside>>iface)%2)<<1);
//?			double a,d,
//?				area,//face areas
//?				vav[DIM],//face-average velocity
//?				*vel,//face-velocity
//?				*var,//face-variable
//?				*norm;//face normal vector (direction determined by direction)
//?			struct Face	*f;
//?			f=face+jf;
//?			norm=f->norm;
//?			area=f->area;
//?			vel=U+DIM*ic;
//?			var=V+DIM*ic;
//?			if (jc>=0)
//?			{	//do upwinding
//?				double
//?					*vi=U+DIM*ic,
//?					*vj=U+DIM*jc;
//?				for (int k=0; k<DIM; k++)vav[k]=0.5*(vi[k]+vj[k]);
//?				a=(double)direction*(SCLP(norm,vav));
//?				if (a<0.0)
//?				{
//?					vel=vj;
//?					var=V+DIM*jc;
//?				}
//?				a=(double)direction*(SCLP(norm,vel));
//?			}
//?			else //jc<0 => boundry face
//?			{
//?				switch(jc)
//?				{
//?					case inlet:
//?						//Get inlet conditions
//?						vel=W+DIM*jf;
//?						var=W+DIM*jf;//DEBUG: should be a variable, not velcity
//?						a=(double)direction*(SCLP(norm,vel));
//?						break;
//?					case outlet:
//?					//Zero normal derivative
//?					if ((a=direction*SCLP(norm,vel))>=0.0)
//?					{
//?						for (int i=0; i<DIM; i++) W[DIM*jf+i]=vel[i];
//?					}
//?					else a=0.0;		
//?					break;
//?					default:
//?						a=0.0;			
//?					break;
//?				}
//?			}
//?			d=a*area; //*onethird;//DEBUG: why onethird?
//?			for (int k=0; k<DIM; k++)
//?				con[k]+=d*var[k];
//?		}
//?		for (int k=0; k<DIM; k++) con[k]*=voli;
//?	}
//?}
//?void	difCellVec
//?//Computes diffusion of a cell-centered vector
//?(
//?	int	nc, //number of cells
//?	int	nbc,//number of boundary cells
//?	struct Cell	*cell,
//?	struct Face	*face,
//?	double	*Vol,//cell-volume
//?	double	*Y,//face-center coordinates
//?	double	*Z,//cell-center coordinates
//?	double	*V,//cell-centered variable
//?	double	*W,//face-vector at the boundary
//?	double	*D //cell-centered diffusion term
//?)
//?{	//Loop over boundary cells
//?//	static double	onethird=1./3.;
//?	int	mc=nc*DIM;
//?	for (int i=0; i<mc; i++)D[i]=0.0;
//?	for (int icell=0; icell<nbc; icell++)
//?	{	struct Cell	*c=cell+icell;
//?		int
//?			imcell=DIM*icell,
//?			*cf=c->face,//cell-face indexes
//?			*cc=c->cell;//neighbor-cell indexes
//?		double
//?			vol=Vol[icell],voli,//cell-volume abd its inverse
//?			*z=Z+imcell,//center coordinates of this cell
//?			*v=V+imcell,//vector variable at this cell
//?			d,   //distance between this cell and its neighbor
//?			di, //and its inverse
//?			*dif=D+imcell;
//?#ifdef DEBUG
//?		if (vol<SMALL)BUG("Cell volume too small in difVectorCell");
//?#endif
//?		voli=1./vol;
//?		for (int lface=0; lface<Nf; lface++)
//?		{	//Face normal derivative
//?			int
//?				jcell=cc[lface],//global index to the neighb cell
//?				//direction,
//?				iface=cf[lface];
//?			double
//?				area,//face area
//?				*z1,//center coordinates of the neighbor cell
//?				*v1;//vector variable at the neighbor cell
//?			struct Face	*f;
//?			//if (iface<0)
//?			//{	//direction=-1;
//?			//	iface=-iface-1;
//?			//}
//?			//else
//?			//	direction=1;
//?			f=face+iface;
//?			area=f->area;
//?			if (jcell>=0)
//?			{	int	jmcell=DIM*jcell;
//?				z1=Z+jmcell;//center coordinates of the neighbor cell
//?				v1=V+jmcell;//vector variable at the neighbor cell
//?			}
//?			else//boundary
//?			{
//?				//	if (W==NULL)
//?				//	{
//?				//		z1=Y+DIM*iface;//coordinates of the neighbor face
//?				//		v1=nullvector;
//?				//	}
//?				//	else
//?				switch(jcell)
//?				{
//?				//	case inlet:
//?				//		//use boundary values
//?				//		//Copy boundary values
//?				//		z1=Y+DIM*iface;
//?				//		v1=W+DIM*iface;
//?				//		break;
//?				case outlet:
//?					//Boundary var = inside var => Neuman
//?					z1=Y+DIM*iface;//boundary face coordinates
//?					v1=V+imcell;//variable at the neighbor cell equals to that of this cell
//?					break;
//?				default:
//?					//use boundary values
//?					//Copy boundary values
//?					z1=Y+DIM*iface;
//?					v1=W+DIM*iface;
//?					break;
//?				}
//?			}
//?			d=0.0;
//?			for (int k=0; k<DIM; k++)
//?			{	double	dk=z1[k]-z[k];
//?				d+=dk*dk;
//?			}
//?			d=sqrt(d);
//?#ifdef DEBUG
//?			if (d<SMALL)BUG("Cell separation too small in difVectorCell");
//?#endif
//?			di=area/d;//*onethird;//DEBUG
//?			for (int i=0; i<DIM; i++)
//?				dif[i]+=(v1[i]-v[i])*di;//face-normal derivative
//?		}
//?		for (int i=0; i<DIM; i++) dif[i]*=voli;
//?	}
//?	for (int icell=nbc; icell<nc; icell++)
//?	{	struct Cell	*c=cell+icell;
//?		int
//?			imcell=DIM*icell,
//?			*cf=c->face,//cell-face indexes
//?			*cc=c->cell;//neighbor-cell indexes
//?		double
//?			vol=Vol[icell],voli,//cell-volume abd its inverse
//?			*z=Z+imcell,
//?			*v=V+imcell,//vector variable at this cell
//?			*dif=D+imcell;
//?#ifdef DEBUG
//?		if (vol<SMALL)BUG("Cell volume too small in difVectorCell");
//?#endif
//?		voli=1./vol;
//?		for (int i=0; i<DIM; i++) dif[i]=0.0;
//?		for (int lface=0; lface<Nf; lface++)
//?		{	//Face normal derivative
//?			int
//?				jcell=cc[lface],//global index to the neighb cell
//?				jmcell=DIM*jcell,
//?				//direction,
//?				iface=cf[lface];
//?			double
//?				d,   //distance between this cell and its neighbor
//?				di, //and its inverse
//?				area,//face area
//?				*z1=Z+jmcell,
//?				*v1=V+jmcell;//vector variable at the neighbor cell
//?			struct Face	*f;
//?			//if (iface<0)
//?			//{	//direction=-1;
//?			//	iface=-iface-1;
//?			//}
//?			//else
//?			//	direction=1;
//?			f=face+iface;
//?			area=f->area;
//?			d=0.0;
//?			for (int k=0; k<DIM; k++)
//?			{	double	dk=z1[k]-z[k];
//?				d+=dk*dk;
//?			}
//?			d=sqrt(d);
//?#ifdef DEBUG
//?			if (d<SMALL)BUG("Cell separation too small in difVectorCell");
//?#endif
//?			di=area/d;//*onethird;//DEBUG
//?			for (int i=0; i<DIM; i++)
//?				dif[i]+=(v1[i]-v[i])*di;//face-normal derivative
//?		}
//?		for (int i=0; i<DIM; i++) dif[i]*=voli;
//?	}
//?}
//?void	delatCellVec2NodeScl0
//?//Computes delatation of a face-vector at the internal nodes
//?// delatation=(\nabla\,,variable)
//?// and convection at the boundary nodes
//?(	
//?	int	nc,//number of cells
//?	int	nn,//number of nodes
//?	int	nbc,//number of boundary cells//DEBUG
//?	struct Cell	*cell,
//?	struct Face	*face,
//?	struct Node	*node,	
//?	struct BoundaryVertex	*Bv,
//?	double	*Vol,//node-based volume
//?	double	*U,//cell-velocity for boundary convection
//?	double	*V,//cell-vector variable which delatation is to be computed
//?	double	*D //scalar node-variable representing the delatation of V
//?)
//?{
//?	static const double TwoThirds3i=27./8.;
//?	fluxCellVec2NodeScl
//?	//Computes delatation of a face-vector at the nodes
//?	// delatation=(\nabla\,,variable)
//?	(	nc,//number of cells
//?		nn,//number of nodes
//?		cell,face,node,Bv,
//?		U,//cell-velocity for boundary convection
//?		V,//cell-vector variable 
//?		D //scalar node-variable representing the delatation of V
//?	);
//?	for (int in=0; in<nn; in++)
//?		D[in]*=TwoThirds3i/Vol[in];//TODO: replace with mul( Voli)
//?}
//?void	delatCellVec2CellScl
//?//Computes delatation of a cell-vector at the cell centers
//?// delatation=(\nabla\,,variable)
//?(	
//?	int	nc,//number of cells
//?	int	nn,//number of nodes
//?	int	nbc,//number of boundary cells
//?	struct Cell	*cell,
//?	struct Face	*face,
//?	struct Node	*node,	
//?	double	*Y,//face-center coordinates
//?	double	*Z,//cell-center coordinates
//?	double	*Volc,//cell-based volume
//?	double	*U,//cell-velocity for boundary convection
//?	double	*V,//cell-vector variable which delatation is to be computed
//?	double	*W,//face-velocity at the boundary
//?	double	*D //scalar cell-variable representing the delatation of V
//?)
//?{	//Set array to zero
//?	//Loop over boundary cells
//?	for (int icell=0; icell<nbc; icell++)
//?	{	struct Cell	*c=cell+icell;
//?		int
//?			imcell=DIM*icell,
//?			*cv=c->vert,
//?			*cf=c->face,
//?			*cc=c->cell;
//?		double
//?			*v=V+imcell,
//?			voli=1./Volc[icell],
//?			delat=0.0;//delatation
//?		for (int jf=0; jf<Nf; jf++)
//?		{	int
//?				//jnode=cv[jf],//global node-index of the current vertex
//?				jface=cf[jf],//face opposite to the vertex
//?				jcell=cc[jf],//neighbor cell across the face
//?				direction=1-(((c->mirror.inside>>jf)%2)<<1);
//?				//defines direction of norm blow
//?			struct Face	*f=face+jface;
//?			double	d,
//?				area=f->area,
//?				dirarea=direction*area,
//?				*norm=f->norm;
//?			if (jcell>=0) //internal face
//?			{	int jmcell=DIM*jcell;
//?				double
//?					vf[DIM],//face interpolated variable
//?					*vv=V+jmcell;
//?				//Center interpolate at the face
//?				for (int i=0; i<DIM; i++)
//?					vf[i]=0.5*(v[i]+vv[i]);
//?				d=(double)dirarea*(SCLP(vf,norm));
//?				//Upwind
//?				//	if(d>=0.0)
//?				//		d=(double)direction*SCLP(v,norm)*area;
//?				//	else
//?				//		d=(double)direction*SCLP(vv,norm)*area;
//?			}
//?			else //boundary face
//?			{	int jmface=DIM*jface;
//?				double 
//?					dx,du[DIM],
//?					*w=W+jmface,
//?					*y=Y+jmface,
//?					*u=U+imcell,
//?					*z=Z+imcell;
//?				switch(jcell)
//?				{
//?					case inlet:
//?					case outlet:
//?						d=(double)dirarea*(SCLP(v,norm));
//?					break;
//?					default:
//?						//implicity approximation of
//?						//derivatives tried below 
//?						//did not work (exploded)
//?						//Distance to the boundary
//?						//dx=0.0;
//?						//for (int i=0; i<DIM; i++)
//?						//{	double	x=y[i]-z[i];
//?						//	dx+=x*x;
//?						//	du[i]=W[i]-U[i];
//?						//	//du[i]=U[i];
//?						//}
//?						//dx=sqrt(dx);
//?						//if (dx>SMALL)
//?						//	d=(double)dirarea*(SCLP(du,norm))/dx;
//?						//else i
//?							d=0.0;
//?					break;
//?				}
//?			}
//?			delat+=d;
//?		}
//?		D[icell]=delat*voli;
//?	}
//?	//Internal cells
//?	for (int icell=nbc; icell<nc; icell++)
//?	{	struct Cell	*c=cell+icell;
//?		int
//?			imcell=DIM*icell,
//?			*cv=c->vert,
//?			*cf=c->face,
//?			*cc=c->cell;
//?			double
//?				*v=V+imcell,
//?				voli=1./Volc[icell],
//?				delat=0.0;
//?		for (int jf=0; jf<Nf; jf++)
//?		{	int
//?				//jnode=cv[jf],//global node-index current vertex
//?				jface=cf[jf],//face opposite to the vertex
//?				jcell=cc[jf],jmcell=DIM*jcell,//neighbor cell across the face
//?				direction=1-(((c->mirror.inside>>jf)%2)<<1);
//?			struct Face	*f=face+jface;
//?			double	d,
//?				area=f->area,
//?				*norm=f->norm,
//?				vf[DIM],//face interpolated variable
//?				*vv=V+jmcell;
//?			//Center interpolate at the face
//?			for (int i=0; i<DIM; i++)
//?				vf[i]=0.5*(v[i]+vv[i]);
//?			d=SCLP(vf,norm);
//?			//Upwind
//?			if(direction*d>=0.0)
//?				d=SCLP(v,norm);
//?			else
//?				d=SCLP(vv,norm);
//?			delat+=(double)direction*area*d;
//?		}
//?		D[icell]=delat*voli;
//?	}
//?}
//?void	delatCellVec2NodeScl
//?//Computes delatation of a cell-vector at the internal nodes
//?// delatation=(\nabla\,,variable)
//?// and convection at the boundary nodes
//?(	
//?	int	nc,//number of cells
//?	int	nn,//number of nodes
//?	int	nbc,//number of boundary cells
//?	struct Cell	*cell,
//?	struct Face	*face,
//?	struct Node	*node,
//?	struct BoundaryVertex	*Bv,
//?	double	*Vol,//node-based volume
//?	//double	*U,//cell-velocity for boundary convection
//?	double	*V,//cell-vector variable which delatation is to be computed
//?	double	*D //scalar node-variable representing the delatation of V
//?)
//?{	//Set array to zero
//?	for (int i=0; i<nn; i++) D[i]=0.0;
//?	//fluxCellVec2NodeSclUp
//?	//Loop over boundary cells
//?//	for (int icell=0; icell<nbc; icell++)
//?	for (int icell=0; icell<nc; icell++)//DEBUG
//?	{	struct Cell	*c=cell+icell;
//?		int
//?			imcell=DIM*icell,
//?			*cv=c->vert,
//?			*cf=c->face,
//?			*cc=c->cell;
//?		double
//?			*v=V+imcell;
//?		for (int jf=0; jf<Nf; jf++)
//?		{	int
//?				jnode=cv[jf],//global node-index of the current vertex
//?				jface=cf[jf],//face opposite to the vertex
//?				jcell=cc[jf],//neighbor cell across the face
//?				direction=1-(((c->mirror.inside>>jf)%2)<<1);
//?				//defines direction of norm blow
//?			struct Face	*f=face+jface;
//?			double	d,
//?				area=f->area,
//?				*norm=f->norm;
//?			if (jcell>=0) //internal face
//?			{	int jmcell=DIM*jcell;
//?				double
//?					vf[DIM],//face interpolated variable
//?					*vv=V+jmcell;
//?				//Center interpolate at the face
//?				for (int i=0; i<DIM; i++)
//?					vf[i]=0.5*(v[i]+vv[i]);
//?				d=(double)direction*area*(SCLP(vf,norm));
//?				//Upwind
//?				//if(d>=0.0)
//?				//	d=(double)direction*SCLP(v,norm)*area;
//?				//else
//?				//	d=(double)direction*SCLP(vv,norm)*area;
//?			}
//?			else //boundary face
//?			{//	double *u;
//?				//	switch(jcell)
//?				//	{
//?					//	case inlet:
//?					//	case outlet:
//?					//	{	//Update the element boundary vertices
//?					//		double	c=SCLP(v,norm);
//?					//		d=c<0.0?0.0:(double)direction*area*c;
//?					//		for (int j=0; j<Nfv; j++)
//?					//			D[cv[(jf+j+1)%Nv]]+=d;
//?					//			//D[cv[(jf+j+1)%Nv]]-=d;//DEBUG
//?					//	}
//?					//	break;
//?					//	default:
//?					{	//d=0.0;
//?						//Update the element boundary vertices
//?						double	c=SCLP(v,norm);
//?						//d=c<0.0?0.0:(double)direction*area*c;
//?						d=(double)direction*area*c;//*2. to account for boundary integral (?)
//?						for (int jv=0; jv<Nfv; jv++)
//?						{	int	jvert=cv[(jf+jv+1)%Nv];
//?							D[jvert]-=d;
//?						}
//?						//d*=Vol[jnode];//to account for boundary integral (?)
//?						//d=0.0;//to account for boundary integral (?)
//?					}
//?					//	break;
//?				//	}
//?			}
//?			//D[jnode]+=d;
//?			D[jnode]+=d;//DEBUG
//?		}
//?	}
//?	//Internal cells
//?	if(0)//DEBUG
//?	for (int icell=nbc; icell<nc; icell++)
//?	{	struct Cell	*c=cell+icell;
//?		int
//?			imcell=DIM*icell,
//?			*cv=c->vert,
//?			*cf=c->face,
//?			*cc=c->cell;
//?		double
//?			*v=V+imcell;
//?		for (int jf=0; jf<Nf; jf++)
//?		{	int
//?				jnode=cv[jf],//global node-index current vertex
//?				jface=cf[jf],//face opposite to the vertex
//?				jcell=cc[jf],jmcell=DIM*jcell,//neighbor cell across the face
//?				direction=1-(((c->mirror.inside>>jf)%2)<<1);
//?			struct Face	*f=face+jface;
//?			double	d,
//?				area=f->area,
//?				*norm=f->norm,
//?				vf[DIM],//face interpolated variable
//?				*vv=V+jmcell;
//?			//Center interpolate at the face
//?			for (int i=0; i<DIM; i++)
//?				vf[i]=0.5*(v[i]+vv[i]);
//?			d=SCLP(vf,norm);
//?			//Upwind
//?			if(direction*d>=0.0)
//?				d=SCLP(v,norm);
//?			else
//?				d=SCLP(vv,norm);
//?			D[jnode]+=(double)direction*area*d;
//?			//D[jnode]-=(double)direction*area*d;//DEBUG
//?		}
//?	}
//?	//Normalize
//?	for (int i=0; i<nn; i++) D[i]/=Vol[i];//to be replaced with '*'
//?}
//?void	mapCellVec2NodeSclBnd
//?//Map cell-center vector boundary node scalar
//?//(used to compute the boundary source term in NS)
//?(
//?	int	nbf,//number of boundary faces
//?	int	nbn,//number of boundary nodes
//?	struct Cell	*cell,
//?	struct Face	*face,
//?	struct Node	*node,	
//?	struct BoundaryVertex	*Bv,
//?	//double	*Volc,//cell-based volume
//?	double	*Barea,//Boundary node-based areas (temporal variable)
//?	double	*V,//cell-vector
//?	double	*B //boundary node scalar representing the normal projection of V
//?)
//?{	for (int i=0; i<nbn; i++)B[i]=Barea[i]=0.0;
//?	//Add contributions from the boundary nodes
//?	for (int iface=0; iface<nbf; iface++)
//?	{	struct Face	*f=face+iface;
//?		int	
//?			*fv=f->vert,
//?			icell=f->cell[0],
//?			ibc=f->cell[1];//boundary condition
//?		//struct Cell	*c=cell+icell;//internal cell
//?		double	d,barea;
//?		//	switch(ibc)
//?		//	{
//?			//	case outlet:
//?			//	case inlet:
//?			{	double
//?				*norm=f->norm,//face-normal vector
//?				//volc=Volc[icell],
//?				*vec=V+DIM*icell;
//?				barea=f->area;
//?				d=SCLP(vec,norm);
//?				//d=-SCLP(vec,norm);//DEBUG
//?				//if(d<0.0) d=0.0;//Upwinding
//?			}
//?			//		break;
//?			//	default:
//?			//		d=0.0;
//?		//	}
//?		for (int iv=0; iv<Nfv; iv++)
//?		{	int
//?				inode=fv[iv],
//?				ibn=node[inode].ind;
//?			B[ibn]+=d*barea;
//?			Barea[ibn]+=barea;
//?		}
//?	}
//?	//Normalize
//?	for (int i=0; i<nbn; i++)B[i]/=Barea[i];
//?}
//?void	delatFaceVec2NodeScl
//?//Computes delatation of a face-vector at the nodes
//?// delatation=(\nabla\,,variable)
//?(	
//?	int	nc,//number of cells
//?	int	nn,//number of nodes
//?	struct Cell	*cell,
//?	struct Face	*face,
//?	double	*Vol,//node-based volume
//?	double	*V,//face-vector variable which gradient is to be computed
//?	double	*D //scalar node-variable representing the delatation of V
//?)
//?{
//?	fluxFaceVec2NodeScl
//?	//Computes delatation of a face-vector at the nodes
//?	// delatation=(\nabla\,,variable)
//?	(	
//?		nc,//number of cells
//?		nn,//number of nodes
//?		cell,
//?		face,
//?		V,//face-vector variable which gradient is to be computed
//?		D //scalar node-variable representing the delatation of V
//?	);
//?	div(nn,D,Vol);//TODO: replace with mul( Voli)
//?}
//?void	fluxFaceVec2NodeScl
//?//Computes delatation of a face-vector at the nodes
//?// delatation=(\nabla\,,variable)
//?(	
//?	int	nc,//number of cells
//?	int	nn,//number of nodes
//?	struct Cell	*cell,
//?	struct Face	*face,
//?	double	*V,//face-vector variable which gradient is to be computed
//?	double	*D //scalar node-variable representing the delatation of V
//?)
//?{
//?	static double	onethird=1./3.;
//?	for (int in=0; in<nn; in++) D[in]=0.0;
//?	for (int ic=0; ic<nc; ic++)
//?	{	struct Cell	*c=cell+ic;
//?		int
//?			*cn=c->vert,//cn[i]=index of vertex i=0:Nv in the node array
//?			*cf=c->face,//cf[i]=index of face i in the global face array
//?			*cc=c->cell;//cc[i]=index of neighbor cell i in the global cell array
//?		for (int i=0; i<Nf; i++)
//?		{	int
//?				jn=cn[i],//global node index
//?				jf=cf[i],//global face index
//?				jc=cc[i],//global neighbor-cell index
//?				direction=1-(((c->mirror.inside>>i)%2)<<1);//=jf>=0?1:-1,
//?			double
//?				d,
//?				area,//face areas
//?				*norm,//face normal vector (direction determined by direction)
//?				*v;
//?			struct Face	*f;
//?			//if(jf>=0)direction=1;
//?			//else{direction=-1;jf=-jf-1;}
//?			f=face+jf;
//?			v=V+DIM*jf;//face-center vector variable
//?			norm=f->norm;
//?			area=f->area;
//?			d=(double)direction*area*(SCLP(v,norm))*onethird;
//?			D[jn]+=d;
//?			if (jc<0)
//?			{	//boundary cell
//?				for (int j=0; j<Nfv; j++)
//?					D[cn[(i+j+1)%Nv]]+=d;
//?			}
//?		}
//?	}
//?}
//?void	fluxCellVec2NodeScl
//?//Computes flux of a cell-vector at the nodes
//?// delatation=(\nabla\,,variable)
//?(	
//?	int	nc,//number of cells
//?	int nn,//number of nodes
//?	struct Cell	*cell,
//?	struct Face	*face,
//?	struct Node	*node,
//?	struct	BoundaryVertex	*Bv,
//?	double	*U,//cell-velocity to compute boundary convection terms
//?	double	*V,//cell-vector variable which delatation is to be computed
//?	//RETURNS:
//?	double	*D //scalar node-variable representing the delatation of V
//?)
//?{	static double	onethird=1./3.;
//?//	static double	twothirds2=4./9.;
//?	static double	twothirds=2./3.;
//?	static double	twonineths=2./9.;
//?	static double Two2oThree3=4./27.;
//?	for (int in=0; in<nn; in++) D[in]=0.0;
//?	for (int ic=0; ic<nc; ic++)
//?	{	struct Cell	*c=cell+ic;
//?		int
//?			*cn=c->vert,//cn[i]=index of vertex i=0:Nv in the node array
//?			*cf=c->face,//cf[i]=index of face i in the global face array
//?			*cc=c->cell;//cc[i]=index of neighbor cell i in the global cell array
//?		double
//?			*u=U+DIM*ic,//cell-center velocity 
//?			*v=V+DIM*ic;//cell-center vector variable
//?		for (int i=0; i<Nv; i++)
//?		{	int
//?				jn=cn[i],//global node index
//?				jf=cf[i],//global face index
//?				jc=cc[i],//global cell index
//?				direction=1-(((c->mirror.inside>>i)%2)<<1);//=jf>=0?1:-1,
//?			double
//?				d,
//?				area,dirarea,//face areas
//?				*norm;//face normal vector (direction determined by direction)
//?			struct Face	*f;
//?			f=face+jf;
//?			area=f->area;
//?			norm=f->norm;
//?			dirarea=Two2oThree3*(double)direction*area;
//?			if (node[jn].type>=0)
//?			{
//?				d=(SCLP(v,norm));
//?			}
//?			else //Boundary node
//?			{	if (node[jn].type!=outlet) 
//?				{	double	*bnorm=Bv[node[jn].ind].norm;//boundary node norm
//?					d=SCLP(u,norm)*SCLP(u,bnorm);
//?				}
//?				else	d=0.0;
//?			}
//?			//	if (jc==outlet)
//?			//	{	//boundary cell
//?			//		double	a=SCLP(u,norm),
//?			//			dd=a*a*twonineths*area;
//?			//		for (int j=0; j<Nfv; j++)
//?			//			D[cn[(i+j+1)%Nv]]+=dd;
//?			//	}
//?			D[jn]+=dirarea*d;
//?		}
//?	}
//?}
//?void gradFaceScl2CellVec
//?//Computes gradient of a face-scalaer at the cell-centers
//?// delatation=(\nabla\,,scalar)
//?(
//?	int	nc,//number of cells
//?	struct Cell	*cell,
//?	struct Face	*face,
//?	double	*Vol,//cell volumes
//?	double	*P,//face-scalar variable which gradient is to be computed
//?	double	*G //cell-center vector representing the gradient of P
//?)
//?{	static double	onethird=1./3.;
//?	for (int ic=0; ic<nc; ic++)
//?	{	struct Cell	*c=cell+ic;
//?		int
//?			*cn=c->vert,//cn[i]=index of vertex i=0:Nv in the node array
//?			*cf=c->face,//cf[i]=index of face i in the global face array
//?			*cc=c->cell;//cc[i]=index of neighbor cell i in the global cell array
//?		double	d=0.0,
//?			*g=G+DIM*ic,//cell-center gradient
//?			vol=Vol[ic],voli;//cell-volume and its inverse
//?#ifdef DEBUG
//?		if(vol<SMALL)BUG("Cell volume too small in gradScalarFaceCell");
//?#endif
//?		voli=1./vol;
//?		for (int i=0; i<DIM; i++)
//?		{	double d=0.0;
//?			for (int iface=0; iface<Nf; iface++)
//?			{	int
//?					jn=cn[iface],//global node index
//?					jf=cf[iface],//global face index
//?					jc=cc[iface],//global neighbor-cell index
//?					direction=1-(((c->mirror.inside>>iface)%2)<<1);
//?				double
//?					dirarea,//face areas
//?					p,//local face-scalar value
//?					*norm;//face normal vector (direction determined by direction)
//?				struct Face	*f;
//?				//if(jf>=0)direction=1;
//?				//else{direction=-1;jf=-jf-1;}
//?				f=face+jf;
//?				p=P[jf];//local face-center variable
//?				norm=f->norm;
//?				dirarea=(double)direction*f->area;
//?				d+=dirarea*norm[i]*p*onethird;
//?			}
//?			g[i]=d*voli;
//?		}
//?	}
//?}
//?void gradCellScl2CellVec
//?//Computes gradient of a face-scalaer at the cell-centers
//?// delatation=(\nabla\,,scalar)
//?(
//?	int	nc,//number of cells
//?	struct Cell	*cell,
//?	struct Face	*face,
//?	double	*Vol,//cell volumes
//?	double	*P,//cell-scalar variable which gradient is to be computed
//?	double	*G //cell-center vector representing the gradient of P
//?)
//?{	static double	onethird=1./3.;
//?	for (int ic=0; ic<nc; ic++)
//?	{	struct Cell	*c=cell+ic;
//?		int
//?			*cn=c->vert,//cn[i]=index of vertex i=0:Nv in the node array
//?			*cf=c->face,//cf[i]=index of face i in the global face array
//?			*cc=c->cell;//cc[i]=index of neighbor cell i in the global cell array
//?		double	d=0.0,
//?			*g=G+DIM*ic,//cell-center gradient
//?			p=P[ic],//cell-center value
//?			vol=Vol[ic],voli;//cell-volume and its inverse
//?#ifdef DEBUG
//?		if(vol<SMALL)BUG("Cell volume too small in gradScalarFaceCell");
//?#endif
//?		voli=1./vol;
//?		for (int i=0; i<DIM; i++)
//?		{	double d=0.0;
//?			for (int iface=0; iface<Nf; iface++)
//?			{	int
//?					jn=cn[iface],//global node index
//?					jf=cf[iface],//global face index
//?					jc=cc[iface],//global neighbor-cell index
//?					direction=1-(((c->mirror.inside>>iface)%2)<<1);
//?				double
//?					dirarea,//face areas
//?					pf,//face-interpolated scalar value
//?					*norm;//face normal vector (direction determined by direction)
//?				struct Face	*f;
//?				//if(jf>=0)direction=1;
//?				//else{direction=-1;jf=-jf-1;}
//?				f=face+jf;
//?				if (jc>=0)
//?				{//Neighbor exists
//?					pf=0.5*(P[jc]+p);//local face-center variable
//?				}
//?				else
//?				{//Boundary cell: no neighbor
//?					pf=p;
//?				}
//?				norm=f->norm;
//?				dirarea=(double)direction*f->area;
//?				d+=dirarea*norm[i]*pf*onethird;
//?			}
//?			g[i]=d*voli;
//?		}
//?	}
//?}
//?void delatFaceVec2CellScl
//?//Computes delatation of a face-vector at the cell-centers
//?// delatation=(\nabla\,,vector)
//?(
//?	int	nc,//number of cells
//?	struct Cell	*cell,
//?	struct Face	*face,
//?	double	*Vol,//cell volumes
//?	double	*V,//face-vector variable which gradient is to be computed
//?	double	*D //scalar cell-variable representing the delatation of V
//?)
//?{
//?	fluxFaceVec2CellScl
//?	//Computes delatation of a face-vector at the cell-centers
//?	// delatation=(\nabla\,,variable)
//?	(
//?		nc,//number of cells
//?		cell,
//?		face,
//?		V,//face-vector variable which gradient is to be computed
//?		D //scalar cell-variable representing the delatation of V
//?	);
//?	div(nc,D,Vol);
//?}
//?void fluxFaceVec2CellScl
//?//Computes delatation of a face-vector at the cell-centers
//?// delatation=(\nabla\,,variable)
//?(
//?	int	nc,//number of cells
//?	struct Cell	*cell,
//?	struct Face	*face,
//?	double	*V,//face-vector variable which gradient is to be computed
//?	double	*D //scalar cell-variable representing the delatation of V
//?)
//?{	static double	onethird=1./3.;
//?	for (int ic=0; ic<nc; ic++)
//?	{	struct Cell	*c=cell+ic;
//?		int
//?			*cn=c->vert,//cn[i]=index of vertex i=0:Nv in the node array
//?			*cf=c->face,//cf[i]=index of face i in the global face array
//?			*cc=c->cell;//cc[i]=index of neighbor cell i in the global cell array
//?		double	d=0.0;
//?		for (int i=0; i<Nf; i++)
//?		{	int
//?				jn=cn[i],//global node index
//?				jf=cf[i],//global face index
//?				jc=cc[i],//global neighbor-cell index
//?				direction=1-(((c->mirror.inside>>i)%2)<<1);
//?			double
//?				area,//face areas
//?				*v,
//?				*norm;//face normal vector (direction determined by direction)
//?			struct Face	*f;
//?			//if(jf>=0)direction=1;
//?			//else{direction=-1;jf=-jf-1;}
//?			f=face+jf;
//?			v=V+DIM*jf;//local face-center variable
//?			norm=f->norm;
//?			area=f->area;
//?			d+=(double)direction*area*(SCLP(v,norm))*onethird;
//?		}
//?		D[ic]=d;
//?	}
//?}
//?//Assemble the stiffness matrix coeffs from the source-terms
//?//../doc/fem.tex
//?void	assemblePoissonNode
//?(
//?	int	nc,//number of cells
//?	int	nn,//number of nodes
//?	int	nbf,//number of boundary faces 
//?	int	nbv,//number of boundary vertexes 
//?	struct Cell	*cell,
//?	struct Face	*face,
//?	struct Node	*node,
//?	struct BoundaryVertex	*Bv,
//?	double	*Vol,//cell-volumes array
//?	double	*S,//source terms (nodal)
//?//	double	*B,//boundary source terms (nodal)
//?	//RETURNS:
//?	double	*A //coefficients to be assembled (nodal)
//?)
//?{//	const double	twothirds=2./3.;
//?	const double	onesixth=1./6.;
//?	for (int in=0; in<nn; in++) A[in]=0.0;
//?	for (int ic=0; ic<nc; ic++)
//?	{	struct Cell	*c=cell+ic;
//?		int	*vert=c->vert;
//?		double
//?			vol=Vol[ic],
//?			c0ii=0.10*vol,
//?			c0ij=0.05*vol;
//?		//Loop through the vertexes
//?		for (int iv=0; iv<Nv; iv++)
//?		{	int	ivert=vert[iv];//this cell's vertex opposite to the current face
//?			double	a=c0ii*S[ivert];
//?			//Stiffness matrix diagonal elements:
//?//DEBUG			if (node[ivert].type<0)continue;
//?			for (int jv=0; jv<Nv1; jv++)
//?			{	int	jvert=vert[(iv+jv+1)%Nv];
//?				//conv[mncon*jvert]+=c0ij*source[jvert*nsrc];
//?				a+=c0ij*S[jvert];
//?			}
//?			A[ivert]+=a;//blew up 0814
//?			//A[ivert]+=a;//blew-up 1103
//?		}
//?	}
//?	//	for (int ibf=0; ibf<nbf; ibf++)
//?	//	{	//Boundary source terms
//?	//		struct Face	*f=face+ibf;
//?	//		int
//?	//			*fv=f->vert;
//?	//		double
//?	//			area=f->area,
//?	//			bii=onesixth*area,///twothirds
//?	//			bij=0.25*area;///0.5
//?	//		for (int iv=0; iv<Nfv; iv++)
//?	//		{	int
//?	//				ivert=fv[iv],//index to the global vertex array
//?	//				ibn=node[ivert].ind;//index to the boundary vertex array
//?	//			double	b=bii*B[ibn];
//?	//			for (int jv=0; jv<Nfv1; jv++)
//?	//			{	int
//?	//					jvert=fv[(iv+jv+1)%Nfv],
//?	//					jbn=node[jvert].ind;
//?	//				b+=bij*B[jbn];
//?	//			}
//?	//			//A[ivert]+=b; //bad direction
//?	//			A[ivert]-=b;//bad direction
//?	//		}
//?	//	}
//?}
//?void	stepPoissonNodeGSNeuman
//?//One Gauss-Seidel iteration of a Poisson equation in FEM discretization
//?(
//?	int	nc,//number of cells
//?	int	nn,//number of nodes
//?	int	nbc,//number of boundary cells
//?	struct Cell	*cell,
//?	struct Face	*face,
//?	struct Node	*node,
//?	double	*Volc,//cell-volumes array
//?	double	*Voln,//node-volumes array
//?	double	*Y,//face-centers
//?	double	*Z,//cell-centers
//?	double	*A,//Source terms (assembpled coeffs ../doc/fem.tex:\ref{poisiter})
//?	double	*B,//Boundary source terms: vector whose face-normal component
//?	         // determines a face-normal derivative of P
//?	double	*C,//Neighbor's contributions
//?	//RETURNS:
//?	double	*P //Node-scalar variable to be solved
//?)
//?{	static const double factor=1./(double)(1-Nfv/Nv);
//?	for (int in=0; in<nn; in++) C[in]=0.0;
//?	//Boundary Cells (Neuman condition)
//?	for (int ic=0; ic<nbc; ic++)
//?	{	struct Cell	*c=cell+ic;
//?		int
//?			mcell=DIM*ic,
//?			*cn=c->vert,
//?			*cc=c->cell,
//?			*cf=c->face;
//?		double
//?			volc=Volc[ic],
//?			*z=Z+mcell,
//?			*b=B+mcell;
//?		for (int jf=0; jf<Nf; jf++)
//?		{
//?			if (cc[jf]<0)
//?			{	//Boundary face
//?				int
//?					jface=cf[jf],
//?					direction=1-(((c->mirror.inside>>jf)%2)<<1);
//?				struct Face	*f=face+jface;
//?				double
//?					*y=Y+DIM*jface,
//?					*norm=f->norm,
//?					pint=P[cn[jf]],//value at internal node
//?					r[DIM],//distance between the face and cell-centers
//?					d2,
//?					snr,//=SCLP(norm,r)
//?					Pface,Bn;
//?				for (int k=0; k<DIM; k++) r[k]=y[k]-z[k];
//?				d2=SCLP(r,r);
//?				snr=direction*SCLP(norm,r);
//?#ifdef DEBUG
//?				if (snr<0.0) BUG("snr NEGATIVE IN stepPoissonNodeGSNeuman");
//?#endif
//?				/*
//?					(Pface-Pmid)/d*SCLP(norm,r)/d=Bn
//?					Bn=SCLP(B(cell),norm)
//?					Pmid=(Pint+Nfv*Pface)/Nv
//?					Pface=(Pint/Nv+Bn*d2/snr)/(1-Nfv/Nv)     [1]
//?					Pface=SUM(P(node),Nfv)/Nfv               [2]
//?				  [1,2]=>	Pnode+=Pface*vol(cell)/vol(node)
//?				 */
//?				Bn=direction*SCLP(norm,b);
//?				Pface=(pint/Nv+Bn*d2/snr)*factor;
//?				for (int iv=0; iv<Nfv; iv++)
//?				{	int
//?						ivert=(jf+iv+1)%Nv,
//?						inode=cn[ivert];
//?					P[inode]+=Pface*volc/Voln[inode];
//?				}
//?			}
//?		}
//?	}
//?	//Internal cells
//?//	for (int ic=nbc; ic<nc; ic++)
//?	for (int ic=0; ic<nc; ic++)
//?	{	struct Cell	*c=cell+ic;
//?		int
//?			*cv=c->vert,//cell-vertexes
//?			*cf=c->face,//cf[i]=index of face i in the global face array
//?			*cc=c->cell;
//?		double
//?			vol=Volc[ic],vol9i;
//?#ifdef DEBUG
//?		if(vol<SMALL)BUG("Volume too small in PoissonNodeGS");
//?#endif
//?		vol9i=1./(9*vol);
//?		for (int i=0; i<Nf; i++)
//?		{	int
//?				icv=cv[i],//global vertex index
//?				icf=cf[i],//global face index
//?				//icc=cc[i],//global neighbor-cell index
//?				diri=1-(((c->mirror.inside>>i)%2)<<1);
//?			double
//?				b,//../doc/fem.tex:\ref{poisiter}
//?				ai,//face area
//?				*ni;//face normal vector (direction determined by direction)
//?			struct Face	*fi;
//?			if (cc[i]<0) continue;//DEBUG
//?			//if(icf>=0)diri=1;
//?			//else{icf=-icf-1;diri=-1;}
//?			fi=face+icf;
//?			ni=fi->norm;
//?			ai=(double)diri*fi->area;
//?			b=0.0;
//?			for (int j=0; j<Nv1; j++)
//?			{	int
//?					ij=(i+j+1)%Nv,
//?					dirj=1-(((c->mirror.inside>>ij)%2)<<1),
//?					jcv=cv[ij],
//?					jcf=cf[ij];//global face index
//?				double
//?					c1ij,
//?					aj,*nj;
//?				struct Face	*fj;
//?				//if(jcf>=0)dirj=1;
//?				//else{jcf=-jcf-1;dirj=-1;}
//?				fj=face+jcf;
//?				nj=fj->norm;
//?				aj=(double)dirj*fj->area;
//?				c1ij=(SCLP(ni,nj))*ai*aj;
//?				b+=c1ij*P[jcv];
//?			}
//?			C[icv]+=b*vol9i;
//?		}
//?	}
//?	//upgrade all nodes \ref{poisiter}
//?	for (int i=0; i<nn; i++)
//?		//if (node[i].type>=0)//!=presinlet)
//?//{printf("%d: B=%g, c1ii=%g, dP=%g, P=",i,B[i],node[i].c1ii,B[i]/node[i].c1ii);
//?			P[i]=(A[i]-C[i])/node[i].c1ii;
//?//printf("%g\n",P[i]);
//?//}//DEBUG
//?}
//?void	stepPoissonNodeGS
//?//One Gauss-Seidel iteration of a Poisson equation in FEM discretization
//?(
//?	int	nc,//number of cells
//?	int	nn,//number of nodes
//?	struct Cell	*cell,
//?	struct Face	*face,
//?	struct Node	*node,
//?	double	*Volc,//cell-volumes array
//?	double	*A,//Source terms (assembpled coeffs ../doc/fem.tex:\ref{poisiter})
//?	//USES FOR TEMPORAL DATA STORAGE:
//?	double	*B,//Array of node-scalars used for diffusion terms
//?	//RETURNS:
//?	double	*P //Node-scalar variable to be solved
//?)
//?{	static const double factor=1./(double)(1-Nfv/Nv);
//?	//Assemble B[*]
//?	for (int in=0; in<nn; in++) B[in]=0.0;
//?	//Internal cells
//?	for (int ic=0; ic<nc; ic++)
//?	{	struct Cell	*c=cell+ic;
//?		int
//?			*cv=c->vert,//cell-vertexes
//?			*cf=c->face;//cf[i]=index of face i in the global face array
//?			//*cc=c->cell;
//?		double
//?			vol=Volc[ic],vol9i;
//?#ifdef DEBUG
//?		if(vol<SMALL)BUG("Volume too small in PoissonNodeGS");
//?#endif
//?		vol9i=1./(9*vol);
//?		for (int i=0; i<Nf; i++)
//?		{	int
//?				ivert=cv[i],//global vertex index
//?				iface=cf[i],//global face index
//?				//icc=cc[i],//global neighbor-cell index
//?				diri=1-(((c->mirror.inside>>i)%2)<<1);
//?			double
//?				b,//../doc/fem.tex:\ref{poisiter}
//?				ai,//face area
//?				*ni;//face normal vector (direction determined by direction)
//?			struct Face	*fi;
//?			fi=face+iface;
//?			ni=fi->norm;
//?			ai=(double)diri*fi->area;
//?			b=0.0;
//?			for (int j=0; j<Nv1; j++)
//?			{	int
//?					ij=(i+j+1)%Nv,
//?					dirj=1-(((c->mirror.inside>>ij)%2)<<1),
//?					jvert=cv[ij],
//?					jface=cf[ij];//global face index
//?				double
//?					c1ij,
//?					aj,*nj;
//?				struct Face	*fj;
//?				fj=face+jface;
//?				nj=fj->norm;
//?				aj=(double)dirj*fj->area;
//?				c1ij=(SCLP(ni,nj))*ai*aj;
//?				b+=c1ij*P[jvert];
//?			}
//?			B[ivert]+=b*vol9i;
//?		}
//?	}
//?	//upgrade all nodes \ref{poisiter}
//?	for (int i=0; i<nn; i++)
//?		if (node[i].type!=presinlet)
//?			P[i]=(A[i]-B[i])/node[i].c1ii;
//?}
//?void	grad
//?(	double	*d[],
//?	double	*b,
//?	double	*g
//?)
//?{	double dd=DET(d);
//?	INV(dd,d,b,g);
//?//		dd=d[0][0]*d[1][1]*d[2][2]+d[0][1]*d[1][2]*d[2][0]+d[1][0]*d[2][1]*d[0][2]
//?//	  -d[2][0]*d[1][1]*d[0][2]-d[0][0]*d[2][1]*d[1][2]-d[1][0]*d[0][1]*d[2][2];
//?//		g[0]=(b[0]*d[1][1]*d[2][2]+d[0][1]*d[1][2]*b[2]+b[1]*d[2][1]*d[0][2]
//?//		     -b[2]*d[1][1]*d[0][2]-b[0]*d[2][1]*d[1][2]-b[1]*d[0][1]*d[2][2])/dd;
//?//		g[1]=(d[0][0]*b[1]*d[2][2]+b[0]*d[1][2]*d[2][0]+d[1][0]*b[2]*d[0][2]
//?//		     -d[2][0]*b[1]*d[0][2]-d[0][0]*b[2]*d[1][2]-d[1][0]*b[0]*d[2][2])/dd;
//?//		g[2]=(d[0][0]*d[1][1]*b[2]+d[0][1]*b[1]*d[2][0]+d[1][0]*d[2][1]*b[0]
//?//		     -d[2][0]*d[1][1]*b[0]-d[0][0]*d[2][1]*b[1]-d[1][0]*d[0][1]*b[2])/dd;
//?}
//?void	gradNodeScl2CellVec
//?//Computes gradient of a nodal scalar at the cell centers
//?(
//?	int	nc,//number of cells
//?	struct Cell	*cell,
//?	struct Face	*face,
//?	double	*X,//nodal coordinates
//?	double	*P,//nodal scalar
//?	//RETURNS:
//?	double	*G //cell-center gradient vector
//?)
//?{		//Node-scalar gradient terms
//?		for (int ic=0; ic<nc; ic++)
//?		{	struct Cell	*c=cell+ic;
//?			int
//?				*vert=c->vert,
//?				i0=vert[0];
//?			double
//?				p0=P[i0],//nodal scalar value at the node 0
//?				*x0=X+DIM*i0,//coordinates of node 0
//?				*g=G+DIM*ic,//gradient at cell-center
//?			// Scalar gradient at point y[0:DIM-1]:
//?			/*  Scalar varlue inside a tetrahedron
//?			 *  is represented by a linear function:
//?			 *  P(x0,x1,x2)=p0+g0*x0+g1*x1+g2*x2=p0+Sum[g[j]*x[j],j]
//?			 *  To find p0,g[i], use Nv=4 relations for values of
//?			 *  P at the vertexes of the tetrahedron:
//?			 *  P(x[i][j=0:DIM-1])=p0+Sum[g[j]*x[i][j],j=0:DIM-1]    (1)
//?			 *  i=0:Nv-1
//?			 *  and solve them for p0,g[j] as unknowns.
//?			 *  In fact, since we need only g[j] we can
//?			 *  reduce the system by one equation by 
//?			 *  subtracting the first equation from the rest:
//?			 *  b[i]=Sum[g[j]*d[i][j],j=0:DIM-1]                    (2)
//?			 *  where now i=0:DIM-1, and
//?			 *  b[i]=P(x[i=1:Nv-1][*])-P(x[0][*])
//?			 *  d[i][j]=x[i+1][j]-x[0][j]; i=0:DIM-1, j=0:DIM-1
//?			 *  Solving (2) for g[j] we get:
//?
//?			      / d00  d01  d02 \  g0   b0
//?			     |  d10  d11  d12  | g1 = b1
//?			      \ d20  d21  d22 /  g2   b2
//?
//?			dd = d00*d11*d22+d01*d12*d20+d10*d21*d02
//?			    -d20*d11*d02-d00*d21*d12-d10*d01*d22
//?			g0 = (b0*d11*d22+d01*d12*b2+b1*d21*d02
//?			     -b2*d11*d02-b0*d21*d12-b1*d01*d22)/dd
//?			g1 = (d00*b1*d22+b0*d12*d20+d10*b2*d02
//?			     -d20*b1*d02-d00*b2*d12-d10*b0*d22)/dd
//?			g2 = (d00*d11*b2+d01*b1*d20+d10*d21*b0
//?			     -d20*d11*b0-d00*d21*b1-d10*d01*b2)/dd
//?			 */
//?			dd,b[DIM],d[DIM][DIM];
//?			for (int i=0; i<DIM; i++)
//?			{	int	iv=vert[i+1];
//?				double	*x=X+DIM*iv;
//?				b[i]=P[iv]-p0;
//?				for (int k=0; k<DIM; k++)
//?					d[i][k]=x[k]-x0[k];
//?			}
//?			//grad(d,b,g);
//?			dd=DET(d);
//?			INV(dd,d,b,g);
//?			//Add nodal scalar gradinet to velocity
//?			//for (int j=0; j<DIM; j++)
//?			//	g[j]=a[j];
//?		}
//?}
//?void	normgradCellVec2FaceScl
//?(	//Computes:
//?	//At internal faces: face-normal derivative of cell-center vector 
//?	//At boundary faces: face-normal projection of cell-center vector
//?	// (\nabla\,,vector)
//?	int	nf,
//?	int	nbf,
//?	struct Cell	*cell,
//?	struct Face	*face,
//?	double	*Y,//face-center coordinates
//?	double	*Z,//cell-center coordinates
//?	double	*V,//cell-center vector
//?	//double	*U,//boundar-face vector for V
//?	//RETURNS:
//?	double	*G //face-normal derivative of V (scalar: rank=0)
//?)
//?{
//?
//?#ifdef XXX
//?
//?	for (int iface=0; iface<nbf; iface++)
//?	{	struct Face	*f=face+iface;
//?		int
//?			mface=DIM*iface,
//?			icell=f->cell[0],
//?			mcell=DIM*icell;
//?		struct Cell	*c=cell+icell;
//?		double	d1[DIM],d2,d,//distance and its square 
//?			*y=Y+mface,
//?			*z=Z+mcell,
//?			*v=V+mcell;
//?			//*u;//vector at the boundary face
//?			//dv[DIM];//difference 
//?		//if (U==NULL)u=nullvector;
//?		//else u=U+mface;
//?		for (int i=0; i<DIM; i++) d1[i]=y[i]-z[i];
//?		d2=SCLP(d1,d1);
//?#ifdef DEBUG
//?		if(d2<SMALL)BUG("Distance too small in normgradCellVec2FaceScl");
//?#endif
//?		G[iface]=-(SCLP(v,d1))/d2;
//?	}
//?
//?#else //XXX
//?
//?	//DEBUG: Face-Normal value 
//?	for (int iface=0; iface<nbf; iface++)
//?	{	struct Face	*f=face+iface;
//?		int
//?			icell=f->cell[0],
//?			mcell=DIM*icell;
//?		double
//?			*norm=f->norm,
//?			*v=V+mcell;
//?		G[iface]=SCLP(v,norm);//DEBUG
//?	}
//?#endif //XXX
//?
//?	//Internal faces
//?	for (int iface=nbf; iface<nf; iface++)
//?	{	struct Face	*f=face+iface;
//?		int
//?			mface=DIM*iface,
//?			icell=f->cell[0],
//?			jcell=f->cell[1],
//?			mcell=DIM*icell,
//?			ncell=DIM*jcell;
//?		//struct Cell	*c=cell+icell;
//?		double	d2,//d,c,//distance and its inverse 
//?			r[DIM],//separation vector
//?			*z=Z+mcell,//this cell coordinates
//?			*y=Z+ncell,//neighbor cell coordinates
//?			*v=V+mcell,//this cell vector
//?			*u=V+ncell,//vector at the neighbor cell
//?			dv[DIM];//difference 
//?		d2=0.0;
//?		for (int i=0; i<DIM; i++)
//?		{	double	a=y[i]-z[i];
//?		 	r[i]=a;
//?			d2+=a*a;
//?			dv[i]=u[i]-v[i];
//?		}
//?		//d=sqrt(d2);
//?#ifdef DEBUG
//?		if(d2<SMALL)BUG("Distance too small in normgradCellVec2FaceScl");
//?#endif
//?		//c=1./d;
//?		G[iface]=(SCLP(dv,r))/d2;
//?	}
//?}
//?void	stepPoissonFaceGS
//?//One iteration of cell-centered poisson solver
//?(
//?	int	nf,
//?	int	nbf,
//?	struct Cell	*cell,
//?	struct Face	*face,
//?	double	*Volc,//cell-center volumes
//?	double	*Volf,//face-center volumes
//?	double	*X, //node-coordinates
//?	double	*Y, //face-coordinates
//?	double	*Z, //cell-coordinates
//?	double	*S, //source terms (face scalar)
//?	double	*P  //variable to solve (face scalar)
//?)
//?{	//Poisson solver for face scalar
//?	//Algorithm as in hedra.tex:ccdisc.tex
//?	//Loop over the boundary faces
//?
//?#ifdef XXX
//?		//Normal gradient condition: comes out smooth
//?	for (int iface=0; iface<nbf; iface++)
//?	{	struct Face	*f=face+iface;
//?		int
//?			mface=DIM*iface,
//?			*fc=f->cell,//face-cell index array
//?			*fv=f->vert,//face-vertex index array
//?			boundary_condition=fc[1];
//?			//ncf[2];//neighbor cell-face: local index of this face
//?						 //at neighbor cell
//?		double
//?			volf=Volf[iface],//volume of a face-based control volume
//?			src=S[iface],//source term
//?			*y=Y+mface,//face coordinates
//?			A,B;//\ref{poisitercc}
//?		A=B=0.0;
//?		if 
//?		(//  boundary_condition!=inlet &&
//?		  boundary_condition!=presinlet &&
//?		  boundary_condition!=presoutlet
//?		)
//?		{	int	ic=0,//internal neighbor cell
//?				*cf,//face index array 
//?				icell=fc[ic],
//?				mcell=DIM*icell,
//?				icf=f->cf;//ncf[ic]; //local index of this face at cell icell
//?			double
//?				*arint=f->a+ic*Nfv*DIM,//areas of int.faces on side ic
//?				*z=Z+mcell;//cell-center coordinates
//?			struct Cell	*c=cell+icell;
//?			cf=c->face;
//?			for (int jf=0; jf<Nf1; jf++)
//?			{	int
//?					//direction=1-(((c->mirror.inside>>jf)%2)<<1),//direction of area[*] (see below)
//?					//ncv=0,//number of common vertexes
//?					//nf=(icf+jf+1)%Nf,//local index of neighbor face
//?					jface=cf[(icf+jf+1)%Nf];//global index of neighbbor face
//?				double	a,
//?					*area=arint+DIM*jf,//area vector of the int.face
//?					*yy,//neighbor face coordinates
//?					r[DIM],//vector connecting this face and neighbor
//?					d;  //face centers, the distance between them and its squre
//?				//if(jface<0){jface=-jface-1;direction=-1;}
//?				//else	direction=1;
//?				yy=Y+DIM*jface; //neighbor face coordinates
//?				for (int i=0; i<DIM; i++)
//?					r[i]=yy[i]-y[i];
//?				d=LENGTH(r);
//?#ifdef DEBUG
//?				if(d<SMALL)BUG("DIVISION BY ZERO IN stepPoissonFaceGS.1");
//?#endif
//?				//SUM[(P[jface]-P[iface])/d*area,jface]/vol=src
//?				//P[iface]=(SUM[P[jface]*area/d,jface]-vol*src)/SUM[area/d,jface]
//?				//direction=SCLP(r,area)>=0?1:-1;
//?				//area=direction*LENGTH(area)/d;
//?				a=LENGTH(area)/d;
//?				A+=a*P[jface];
//?				B+=a;
//?			}
//?		}
//?//CH:FLAT_OUTLET		P[iface]=(A-volf*src)/B;
//?		P[iface]=(-A+volf*src)/B;
//?	}
//?
//?#else//XXX
//?
//?	//Face-normal derivative condition
//?	//at the boundary
//?	for (int iface=0; iface<nbf; iface++)
//?	{	struct Face	*f=face+iface;
//?		int
//?			mface=DIM*iface,
//?			*fc=f->cell,//face-cell index array
//?			boundary_condition=fc[1];
//?
//?			//ncf[2];//neighbor cell-face: local index of this face
//?						 //at neighbor cell
//?		if 
//?		(//  boundary_condition!=inlet &&
//?		  boundary_condition!=presinlet &&
//?		  boundary_condition!=presoutlet
//?		)
//?		{	int	ic=0,//internal neighbor cell
//?				*cf,//face index array 
//?				icell=fc[ic],
//?				mcell=DIM*icell,
//?				icf=f->cf;//ncf[ic]; //local index of this face at cell icell
//?			double 
//?				//volc=Volc[icell],
//?				volf=Volf[iface],//volume of a face-based control volume
//?*y=Y+mface,//DEBUG
//?*z=Z+mcell,//DEBUG cell-center coordinates
//?				src=S[iface],//source term
//?				area=f->area,
//?				*norm=f->norm,
//?				*arint=f->a+ic*Nfv*DIM,//areas of int.faces on side ic
//?				A,B;//\ref{poisitercc}
//?			struct Cell	*c=cell+icell;
//?			cf=c->face;
//?			A=B=0.0;
//?			for (int jf=0; jf<Nf1; jf++)
//?			{//Looping through other faces
//?				int
//?					kf=(icf+jf+1)%Nf,//local neighbor face index
//?					//dirj=1-(((c->mirror.inside>>kf)%2)<<1),
//?					jface=cf[kf];//global neighbbor face index
//?				struct Face	*fj=face+jface;
//?				double	a,
//?					*aint=arint+DIM*jf;//area vector of the int.face
//?					//areaj=fj->area,//area vector of the int.face
//?					//*normj=fj->norm;
//?				a=SCLP(norm,aint);
//?
//?				A+=a*P[jface];
//?				B+=a;
//?			}
//?
//?			P[iface]=(volf*src-0.5*A)/(area+0.5*B);
//?		}
//?	}
//?
//?#endif//XXX
//?
//?	//Loop over the internal faces
//?	for (int iface=nbf; iface<nf; iface++)
//?	{	struct Face	*f=face+iface;
//?		int
//?			mface=DIM*iface,
//?			*fc=f->cell,//face-cell index array
//?			ncf[2];//neighbor cell-face: local index of this face
//?			       //at neighbor cell
//?		double
//?			volf=Volf[iface],//volume of a face-based control volume
//?			*y=Y+mface,//face coordinates
//?			src=S[iface],//source term
//?			A,B;//\ref{poisitercc}
//?		ncf[0]=f->cf%Nf;//local index of this face at neighbor cell[0]
//?		ncf[1]=f->cf/Nf;// the same for neighb. cell[1]
//?		A=B=0.0;
//?		for (int ic=0; ic<2; ic++)
//?		{	int
//?				*cf,//face index array 
//?				icell=fc[ic],
//?				mcell=DIM*icell,
//?				icf=ncf[ic]; //local index of this face at cell icell
//?			double
//?				*arint=f->a+ic*Nfv*DIM,//areas of int.faces on side ic
//?				*z=Z+mcell;//cell-center coordinates
//?			struct Cell	*c=cell+icell;
//?			cf=c->face;
//?			for (int jf=0; jf<Nf1; jf++)
//?			{	int
//?					kf=(icf+jf+1)%Nf,//local index of neighbor face
//?					//direction=1-(((c->mirror.inside>>kf)%2)<<1),//direction of area[*] (see below)
//?					//ncv=0,//number of common vertexes
//?					jface=cf[kf];//global index of neighbbor face
//?				double	a,
//?					*area=arint+DIM*jf,//area vector of the int.face
//?					pp,//neighbor face scalar variable values
//?					*yy,//neighbor face coordinates
//?					r[DIM],//vector connecting this face and neighbor
//?					d;  //face centers, the distance between them and its squre
//?				//		*ffv,//neighbor face vertex array
//?				//		a[DIM],//face area vector of the inter-face between 
//?				//		area,  //this face and the neighbor face and its scalar value
//?				//		e[2][DIM],//edges of a face separating this and neighbor face
//?				//		*x[2];//Two common vertexes of faces iface and jface 
//?				//struct Face	*ff;
//?				//if(jface<0){jface=-jface-1;direction=-1;}
//?				//else	direction=1;
//?				//ff=face+jface;
//?				pp=P[jface];//neighbor face scalar variable values
//?				yy=Y+DIM*jface; //neighbor face coordinates
//?				for (int i=0; i<DIM; i++)
//?					r[i]=yy[i]-y[i];
//?				d=LENGTH(r);
//?#ifdef DEBUG
//?				if(d<SMALL)BUG("DIVISION BY ZERO IN stepPoissonFaceGS.1");
//?#endif
//?				a=LENGTH(area)/d;
//?				A+=a*pp;
//?				B+=a;
//?			}
//?		}
//?		P[iface]=(A+volf*src)/B;//BAD-CH: GOOD-POISSON
//?	}
//?}
//?void	stepPoissonCellGS
//?(
//?	int	nc,//number of cells
//?	struct Cell	*cell,
//?	struct Face	*face,
//?	double	*Vol,//cell-volume
//?	double	*Z,//cell-center coordinates
//?	double	*Src, //cell-center source-term
//?	double	*Val //variable to solve 
//?)
//?{
//?static double	onethird=1./3.;
//?	//Poisson solver for central scalar
//?	//Algorithm as in hedra.tex:ccdisc.tex
//?	for (int icell=0; icell<nc; icell++)
//?	{	struct Cell	*c=cell+icell;
//?		int	
//?			*face_index=c->face,
//?			*neighb_cell_index=c->cell;
//?		double
//?			*z=Z+DIM*icell,//cell-center coordinates
//?			*val=Val+icell,//local variable value
//?			src=Src[icell],//local source term
//?		//	voli,//inverse cell volume
//?			vol=Vol[icell],
//?			A,B;//ccdisc.tex:\ref{poisitercc}
//?#ifdef DEBUG
//?				if(vol<SMALL)BUG("DIFFUSION TERM: GRID CELL VOLUME TOO SMALL");
//?#endif
//?	//	voli=1./Vol[icell];
//?		A=B=0.0;
//?		for (int iface=0; iface<Nf; iface++)
//?		{	int
//?				jface=face_index[iface],//face index in face-array 
//?				jcell=neighb_cell_index[iface];//cell index in cell-array
//?			struct Face	*f=face+jface;
//?			double
//?				val1,//local value at the neighbor cell
//?				//*norm,
//?				a,d,d2;
//?			if (jcell<0)
//?			{	
//?				continue;
//?			}
//?			val1=Val[jcell];
//?			//norm=f->norm;
//?			d2=0.0;//squared distance between the neighbor cells
//?			for (int i=0; i<DIM; i++)
//?			{	double	di=Z[DIM*jcell+i]-z[i];
//?				d2+=di*di;
//?			}
//?			d=sqrt(d2);
//?#ifdef DEBUG
//?			if(d<SMALL)BUG("FLOATING OVERFLOW in stepPoissonCellGS");
//?#endif
//?			a=f->area/d;//*onethird;//DEBUG
//?			A+=a*val1;
//?			B+=a;
//?		}
//?		*val=(A-vol*src)/B;
//?	}
//?}
void	subvec
(	//A=(A-B)*d
	int	n,
	double	*A,
	double	*B,
	double	d
)
{
	for (int i=0; i<n; i++)
	{	int	m=DIM*i;
		double
			*a=A+m,
			*b=B+m;
		for (int j=0; j<DIM; j++)
			a[j]-=b[j]*d;
	}
}
void	addvec
(	//A=(A+B)*d
	int	n,
	double	*A,
	double	*B,
	double	d
)
{
	for (int i=0; i<n; i++)
	{	int	m=DIM*i;
		double
			*a=A+m,
			*b=B+m;
		for (int j=0; j<DIM; j++)
			a[j]+=b[j]*d;
	}
}
void	addvec
(	//A=(A+B)*d
	int	n,
	double	*A,
	double	*B
)
{
	for (int i=0; i<n; i++)
	{	int	m=DIM*i;
		double
			*a=A+m,
			*b=B+m;
		for (int j=0; j<DIM; j++)
			a[j]+=b[j];
	}
}
void	subvec
(	//A=(A+B)*d
	int	n,
	double	*A,
	double	*B
)
{
	for (int i=0; i<n; i++)
	{	int	m=DIM*i;
		double
			*a=A+m,
			*b=B+m;
		for (int j=0; j<DIM; j++)
			a[j]-=b[j];
	}
}
void	mul
(	//A=(A+B)*d
	int	n,
	double	*A,
	double	c
)
{
	for (int i=0; i<n; i++)
	{	int	m=DIM*i;
		double
			*a=A+m;
		for (int j=0; j<DIM; j++)
			a[j]+=c;
	}
}
void	div(int n, double *A, double *B)
{
#ifdef DEBUG
	for (int i=0; i<n; i++)
		if (fabs(B[i])<SMALL)BUG("DIVISION BY ZERO in div");
#endif
	for (int i=0; i<n; i++)
		A[i]/=B[i];
}
void	zerovec(int n, double *A)
{
	for (int i=0; i<DIM*n; i++) A[i]=0.0;
}
void	zero(int n, double *A)
{
	for (int i=0; i<n; i++) A[i]=0.0;
}
//	void	getvec///DDD: average
//	(
//		int	vecloc,//offset of vector variable to be interpolated
//		DCell	*cell,//cell where the variable is located
//		double	*X,//point coordinates: X[0:DIM-1]
//		//RETURNS:
//		double	*V //vector at the point: V[0:DIM-1]
//	)
//	{ DNode	**verts;
//		if (cell==NULL) return;
//		verts=cell->vert;
//		for (int i=0; i<DIM; i++)
//		{	int	loc=vecloc+i;
//			double	v=0.0;
//			for (int j=0; j<Nv; j++)
//				v+=verts[j]->var[loc];;
//			V[i]=v/(double)Nv;
//		}
//	}
void	getscl
(
 	int	loc,//offset of vector variable to be interpolated
	DCell	*cell,//cell where the variable is located
	double	*X,//point coordinates: X[0:DIM-1]
	//RETURNS:
	double	*scl //sclar at the point X to be returned
)
{ DNode	**verts;
	double	dd,d[Nv1][DIM],g[DIM],dx[DIM],
		*x0;
	if (cell==NULL) return;
	verts=cell->vert;
	x0=verts[0]->x;
	for (int i=0;i<DIM;i++)
		dx[i]=X[i]-x0[i];
	for (int j=0; j<Nv1; j++)
	{	double	*x=verts[j+1]->x;
		for (int i=0; i<DIM; i++)
			d[j][i]=x[i]-x0[i];
	}
	dd=DET(d);
	{	double	b[DIM],
			s0=verts[0]->var[loc];
		for (int j=0; j<Nv1; j++)
		{	double	s=verts[j+1]->var[loc];
			b[j]=s-s0;
		}
		INV(dd,d,b,g);
		*scl=s0+SCLP(g,dx);
	}
}
void	getvec
(
 	int	vecloc,//offset of vector variable to be interpolated
	DCell	*cell,//cell where the variable is located
	double	*X,//point coordinates: X[0:DIM-1]
	//RETURNS:
	double	*V //vector at the point X: V[0:DIM-1]
)
{ DNode	**verts;
	double	dd,d[Nv1][DIM],g[DIM],dx[DIM],
		*x0;
	if (cell==NULL) return;
	verts=cell->vert;
	x0=verts[0]->x;
	for (int i=0;i<DIM;i++)
		dx[i]=X[i]-x0[i];
	for (int j=0; j<Nv1; j++)
	{	double	*x=verts[j+1]->x;
		for (int i=0; i<DIM; i++)
			d[j][i]=x[i]-x0[i];
	}
	dd=DET(d);
	for (int i=0; i<DIM; i++)
	{	int	loc=vecloc+i;
		double	b[DIM],
			v0=verts[0]->var[loc];
		for (int j=0; j<Nv1; j++)
		{	double	v=verts[j+1]->var[loc];
			b[j]=v-v0;
		}
		INV(dd,d,b,g);
		V[i]=v0+SCLP(g,dx);
	}
}
void	getvarbuf
(
 	int	nbuf,//offset of vector variable to be interpolated
	DCell	*cell,//cell where the variable is located
	double	*X,//point coordinates: X[0:DIM-1]
	//RETURNS:
	double	*V //vector at the point X: V[0:DIM-1]
)
{ DNode	**verts;
	double	dd,d[Nv1][DIM],g[DIM],dx[DIM],
		*x0;
	if (cell==NULL) return;
	verts=cell->vert;
	x0=verts[0]->x;
	for (int i=0;i<DIM;i++)
		dx[i]=X[i]-x0[i];
	for (int j=0; j<Nv1; j++)
	{	double	*x=verts[j+1]->x;
		for (int i=0; i<DIM; i++)
			d[j][i]=x[i]-x0[i];
	}
	dd=DET(d);
	for (int i=0; i<nbuf; i++)
	{	int	loc=i;
		double	b[DIM],
			v0=verts[0]->var[loc];
		for (int j=0; j<Nv1; j++)
		{	double	v=verts[j+1]->var[loc];
			b[j]=v-v0;
		}
		INV(dd,d,b,g);
		V[i]=v0+SCLP(g,dx);
	}
}
//?void	getvec
//?(
//?	Cell	*cell,//cell where the variable is located
//?	double	*xp,//point coordinates vector: xp[0:DIM-1]
//?	double	*X,//node coordinates array: X[0:DIM*nn-1]
//?	double	*V,//node vectors array: V[0:DIM*nn-1]
//?	//RETURNS:
//?	double	*v //vector at the point X: V[0:DIM-1]
//?)
//?{ int	*verts;
//?	double	dd,d[Nv1][DIM],g[DIM],dx[DIM],
//?		*x0;//zero-th vertex
//?	if (cell==NULL) return;
//?	verts=cell->vert;
//?	x0=X+DIM*verts[0];
//?	for (int i=0;i<DIM;i++)
//?		dx[i]=xp[i]-x0[i];
//?	for (int j=0; j<Nv1; j++)
//?	{	double	*x=X+DIM*verts[j+1];
//?		for (int i=0; i<DIM; i++)
//?			d[j][i]=x[i]-x0[i];
//?	}
//?	dd=DET(d);
//?	for (int i=0; i<DIM; i++)
//?	{	double	b[DIM],
//?			v0=V[DIM*verts[0]+i];
//?		for (int j=0; j<Nv1; j++)
//?		{	double	v1=V[DIM*verts[j+1]+i];
//?			b[j]=v1-v0;
//?		}
//?		INV(dd,d,b,g);
//?		v[i]=v0+SCLP(g,dx);
//?	}
//?}
