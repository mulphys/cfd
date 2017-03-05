/* 
 * Differential Algebra
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

///#define PRESSURE_VELOCITY_COUPLING
#define HYBRID
//#define CENTRAL_DIFFERENCING

void	Domain::VijUj
(	//Differentiating Vij with respect to Xj: Nabla_j.V_ij
	int	iVij,    // memory offset of a node-stored tensor Vij
	int	iUj,
	int	iVijUj // memery offset of a node-stored vector dVijdXj
)
{	DNode	*node;
	if(dnode_root==NULL)return;
	node=dnode_root;
	do
	{	double	*var=node->var,
			*vij=var+iVij,
			*uj=var+iUj,
			*vijuj=var+iVijUj;
		for(int i=0;i<DIM;i++)
		{	double 
				*vj=vij+i*DIM;
///				s=0.0;
///			for(int j=0;j<DIM;j++)
///				s+=vj[j]*uj[j];
			vijuj[i]=SCLP(vj,uj);
		}
		node=node->next;
	}	while(node!=dnode_root);
}
void	Domain::dVdXi
(	//Differentiating V with respect to Xi
	int	ivoln, // memory offset of a nodal volume variable
	int	iV,    // memory offset of a node-stored scalar V
	int idVdXi // memory offset of a node-stored vector dVdXi
)
{	// Gradient of a scalar:
	// Input: scalar; Output: vector
	DNode	*node;
	DCell	*cell;
	// Set dVdXi to zero:
	if(dnode_root==NULL)return;
	node=dnode_root;
	do
	{	double	*var=node->var+idVdXi;
		for(int i=0;i<DIM;i++)var[i]=0.0;
		node=node->next;
	}	while(node!=dnode_root);
	//Loop through the cells and 
	// assemble derivatives
	if(dcell_root==NULL)return;
	cell=dcell_root;
	do
	{	DNode	**verts=cell->vert;
		DCell	**neibs=(DCell**)cell->neib;
		ElementStatus	*facetypes=cell->facetype;
		double	///volc=cell->var[ivolc],volci=0.5/volc,
			e[Nv][DIM];//tetrahedral edge-vectors
		EDGES(e,verts);//computing edge-vectors
		for(int iv=0;iv<Nv;iv++)
		{	double	a[DIM],Vf,
				*var=verts[iv]->var,
				*dVdX=var+idVdXi;
			DCell	*neib=neibs[iv];
			ElementStatus	facetype=facetypes[iv];
			AREA2(iv,e,a);//computing face area vector (actually 2*area)
			//Compute V at the center of the face
			Vf=0.0;
			for(int jv=1;jv<Nv;jv++)
				Vf+=verts[(iv+jv)%Nv]->var[iV];	
			Vf/=Nv1;
#ifdef PRESSURE_VELOCITY_COUPLING
			if(facetype==internal)
///			if(neib!=NULL)
			{	int	jv=(int)(mask[iv]&cell->index)>>((uint)iv<<1);
					//jv is the index of the opposite node to the face iv
					// as addressed in the neighbor cell
				double	*var_neib=neibs[iv]->vert[jv]->var;

///DCell	**neibneibs=(DCell**)neibs[iv]->neib;
///for(jv=0;jv<Nv;jv++)if(neibneibs[jv]==cell)break;
///if(jv>=Nv){fprintf(stderr,"0:jv=%d\n",jv);exit(1);}
///var_neib=neibs[iv]->vert[jv]->var;

				Vf=0.5*(Vf+0.5*(var[iV]+var_neib[iV]));
			}
#else
#ifdef CENTRAL_DIFFERENCING
			if(facetype==internal)
///			if(neib!=NULL)
			{	int	jv=(int)(mask[iv]&cell->index)>>((uint)iv<<1);
					//jv is the index of the opposite node to the face iv
					// as addressed in the neighbor cell
				double	*var_neib=neibs[iv]->vert[jv]->var;

///DCell	**neibneibs=(DCell**)neibs[iv]->neib;
///for(jv=0;jv<Nv;jv++)if(neibneibs[jv]==cell)break;
///if(jv>=Nv){fprintf(stderr,"0:jv=%d\n",jv);exit(1);}
///var_neib=neibs[iv]->vert[jv]->var;

				Vf=0.5*(var[iV]+var_neib[iV]);
			}
			else
				Vf=var[iV];
#endif
#endif
///				if(facetype==outlet)///DDD
///					Vf=var[iV];
			//Contribution to dVdX:
			for(int i=0;i<DIM;i++)
				dVdX[i]+=Vf*a[i];
			if(facetype!=internal)
///			if(neib==NULL)
			{	//Add contribution to the boundary nodes
				for(int jv=1;jv<Nv;jv++)
				{	double	*dVdXv=verts[(iv+jv)%Nv]->var+idVdXi;
					for(int i=0;i<DIM;i++)
						dVdXv[i]+=Vf*a[i];
				}
			}
		}
		cell=cell->next;
	}	while(cell!=dcell_root);
	node=dnode_root;
	do//Normalize by volume to get the derivatives right
	{	double
			*var=node->var,
			vol=var[ivoln],voli=0.5/vol,//coeff 0.5 accounts for 2*area
			*dVdX=var+idVdXi;
			for(int i=0;i<DIM;i++)
				dVdX[i]*=voli;
		node=node->next;
	}	while(node!=dnode_root);
}
void	Domain::dVidXi
(	//Differentiating Vi with respect to Xi: Nabla_i.V_i
	int	ivoln,  // memory offset of a nodal volume variable
	int	iVi,    // memory offset of a node-stored vector Vi
	int idVidXi // memory offset of a node-stored scalar dVdXi
)
{	// Delatation of a vector:
	// Input: vector; Output: scalar
	DNode	*node;
	DCell	*cell;
	// Set dVdXi to zero:
	if(dnode_root==NULL)return;
	node=dnode_root;
	do
	{	node->var[idVidXi]=0.0;
		node=node->next;
	}	while(node!=dnode_root);
	//Loop through the cells and 
	// assemble derivatives
	if(dcell_root==NULL)return;
	cell=dcell_root;
	do
	{	DNode	**verts=cell->vert;
		DCell	**neibs=(DCell**)cell->neib;
		ElementStatus	*facetypes=cell->facetype;
		double	///volc=cell->var[ivolc],volci=0.5/volc,
			e[Nv][DIM];//tetrahedral edge-vectors
		EDGES(e,verts);//computing edge-vectors
		for(int iv=0;iv<Nv;iv++)
		{	double	a[DIM],Vf[DIM],d,
				*var=verts[iv]->var,
				*dVdX=var+idVidXi;
			DCell	*neib=neibs[iv];
			ElementStatus	facetype=facetypes[iv];
			AREA2(iv,e,a);//computing face area vector (actually 2*area)
			//Compute V at the center of the face
			for(int i=0;i<DIM;i++)
			{	int	j=iVi+i;
				double	v=0.0;
				for(int jv=1;jv<Nv;jv++)
					v+=verts[(iv+jv)%Nv]->var[j];
				Vf[i]=v/Nfv;
			}
#ifdef CENTRAL_DIFFERENCING
			if(facetype==internal)
///			if(neib!=NULL)
			{	int	jv=(int)(mask[iv]&cell->index)>>((uint)iv<<1);
					//jv is the index of the opposite node to the face iv
					// as addressed in the neighbor cell
				double	*var_neib=neibs[iv]->vert[jv]->var;

///DCell	**neibneibs=(DCell**)neibs[iv]->neib;
///for(jv=0;jv<Nv;jv++)if(neibneibs[jv]==cell)break;
///if(jv>=Nv){fprintf(stderr,"2:jv=%d\n",jv);exit(1);}
///var_neib=neibs[iv]->vert[jv]->var;

				for(int i=0;i<DIM;i++)
					Vf[i]=0.5*(var[iVi+i]+var_neib[iVi+i]);
			}
			else///DDD_0
				for(int i=0;i<DIM;i++)
					Vf[i]=var[iVi+i];
#endif
///				if(facetype==outlet)///DDD
///					for(int i=0;i<DIM;i++)
///						Vf[i]=var[iVi+i];
			//Contribution to dVdX:
			d=SCLP(Vf,a);
			*dVdX+=d;
///			if(neib==NULL)
			if(facetype!=internal)
			{	//Contribution to the boundary nodes
				for(int jv=1;jv<Nv;jv++)
					verts[(iv+jv)%Nv]->var[idVidXi]+=d;
			}
		}
		cell=cell->next;
	}	while(cell!=dcell_root);
	node=dnode_root;
	do//Normalize by volume to get the derivatives right
	{	double
			*var=node->var,
			vol=var[ivoln],voli=0.5/vol;//coeff 0.5 accounts for 2*area
			var[idVidXi]*=voli;
		node=node->next;
	}	while(node!=dnode_root);
}
void	Domain::dVidXj
(	//Differentiating Vi with respect to Xj
	int	ivoln, // memory offset of a nodal volume variable
	int	iVi,    // memory offset of a node-stored vector Vi
	int idVidXj // momory offset of a node-stored tensor dVidXj
)
{	// Differentiation of a vector
	// Input: vector; Ouput: tensor of rank 2
	const int	dim2=DIM*DIM;
	DNode	*node;
	DCell	*cell;
	// Set dVdXi to zero:
	if(dnode_root==NULL)return;
	node=dnode_root;
	do
	{	double	*var=node->var+idVidXj;
		for(int i=0;i<dim2;i++)
				var[i]=0.0;
		node=node->next;
	}	while(node!=dnode_root);
	//Loop through the cells and 
	// assemble derivatives
	if(dcell_root==NULL)return;
	cell=dcell_root;
	do
	{	DNode	**verts=cell->vert;
		DCell	**neibs=(DCell**)cell->neib;
		ElementStatus	*facetypes=cell->facetype;
		double	///volc=cell->var[ivolc],volci=0.5/volc,
				e[Nv][DIM];//tetrahedral edge-vectors
		EDGES(e,verts);//computing edge-vectors
		for(int iv=0;iv<Nv;iv++)
		{	double	a[DIM],Vf[DIM],
				*var=verts[iv]->var,
				*dVdX=var+idVidXj;
			DCell	*neib=neibs[iv];
			ElementStatus	facetype=facetypes[iv];
			AREA2(iv,e,a);//computing face area vector (actually 2*area)
			//Compute V at the center of the face
///				if(facetypes[iv]==outlet)
///				{	double	*V=var+iVi;
///					if(SCLP(V,a)>0.0)
///						for(int i=0;i<DIM;i++) Vf[i]=V[i];
///					else
///						for(int i=0;i<DIM;i++) Vf[i]=0.0;
///				}
///				else
///				{				
				for(int i=0;i<DIM;i++)
				{	int	j=iVi+i;
					double	v=0.0;
					for(int jv=1;jv<Nv;jv++)
						v+=verts[(iv+jv)%Nv]->var[j];
					Vf[i]=v/Nv1;
				}
///				}
#ifdef CENTRAL_DIFFERENCING_DIFFUSION
	//BAD FOR DIFFUSION
			if(facetype==internal)
///			if(neib!=NULL)
			{	int	jv=(int)(mask[iv]&cell->index)>>((uint)iv<<1);
					//jv is the index of the opposite node to the face iv
					// as addressed in the neighbor cell
				double	*var_neib=neibs[iv]->vert[jv]->var;

///DCell	**neibneibs=(DCell**)neibs[iv]->neib;
///for(jv=0;jv<Nv;jv++)if(neibneibs[jv]==cell)break;
///if(jv>=Nv){fprintf(stderr,"3:jv=%d\n",jv);exit(1);}
///var_neib=neibs[iv]->vert[jv]->var;


				for(int i=0;i<DIM;i++)
					Vf[i]=0.5*(var[iVi+i]+var_neib[iVi+i]);
			}
#endif
			//Contribution to dVdX:
			for(int i=0;i<DIM;i++)
			{	int	k=i*DIM;
				for(int j=0;j<DIM;j++)
					dVdX[k+j]+=Vf[i]*a[j];
			}
			if(facetype!=internal)
///			if(neib==NULL)
			{	//Contribution to the boundary nodes
				for(int jv=1;jv<Nv;jv++)
				{	double	*dVdXv=verts[(iv+jv)%Nv]->var+idVidXj;
					for(int i=0;i<DIM;i++)
					{	int	k=i*DIM;
						for(int j=0;j<DIM;j++)
							dVdXv[k+j]+=Vf[i]*a[j];
					}
				}
			}
		}
		cell=cell->next;
	}	while(cell!=dcell_root);
	node=dnode_root;
	do//Normalize by volume to get the derivatives right
	{	double
			*var=node->var,
			vol=var[ivoln],voli=0.5/vol,//coeff 0.5 accounts for 2*area
			*dVdX=var+idVidXj;
		for(int i=0;i<dim2;i++)
			dVdX[i]*=voli;
		node=node->next;
	}	while(node!=dnode_root);
}
void	Domain::dVijdXj
(	//Differentiating Vij with respect to Xj: Nabla_j.V_ij
	int	ivoln, // memory offset of a nodal volume variable
	int	iVij,    // memory offset of a node-stored tensor Vij
	int	idVijdXj // memery offset of a node-stored vector dVijdXj
)
{	// Delatation of a second rank tensor:
	// Input: tensor of rank 2; Output: vector
	const	int	dim2=DIM*DIM;
	DNode	*node;
	DCell	*cell;
	// Set dVdXi to zero:
	if(dnode_root==NULL)return;
	node=dnode_root;
	do
	{	double	*var=node->var+idVijdXj;
		for(int i=0;i<DIM;i++)
			var[i]=0.0;
		node=node->next;
	}	while(node!=dnode_root);
	//Loop through the cells and 
	// assemble derivatives
	if(dcell_root==NULL)return;
	cell=dcell_root;
	do
	{	DNode	**verts=cell->vert;
		DCell	**neibs=(DCell**)cell->neib;
		ElementStatus	*facetypes=cell->facetype;
		double	e[Nv][DIM];//tetrahedral edge-vectors
		EDGES(e,verts);//computing edge-vectors
		for(int iv=0;iv<Nv;iv++)
		{	double	a[DIM],Vf[dim2],d[DIM],
				*var=verts[iv]->var,
				*dVdX=var+idVijdXj;
			DCell	*neib=neibs[iv];
			ElementStatus	facetype=facetypes[iv];
			AREA2(iv,e,a);//computing face area vector (actually 2*area)
			//Compute V at the center of the face
			for(int i=0;i<DIM;i++)
			{	int	k=DIM*i,ii=iVij+k;
				for(int j=0;j<DIM;j++)
				{	double	v=0.0;
					for(int jv=1;jv<Nv;jv++)
						v+=verts[(iv+jv)%Nv]->var[ii+j];
					Vf[k+j]=v/Nfv;
				}
			}
#ifdef CENTRAL_DIFFERENCING_DIFFUSION
			//BAD FOR DIFFUSION
			if(facetype==internal)
///			if(neib!=NULL)
			{	int	jv=(int)(mask[iv]&cell->index)>>((uint)iv<<1);
					//jv is the index of the opposite node to the face iv
					// as addressed in the neighbor cell
				double	*var_neib=neibs[iv]->vert[jv]->var;

///DCell	**neibneibs=(DCell**)neibs[iv]->neib;
///for(jv=0;jv<Nv;jv++)if(neibneibs[jv]==cell)break;
///if(jv>=Nv){fprintf(stderr,"4:jv=%d\n",jv);exit(1);}
///var_neib=neibs[iv]->vert[jv]->var;


				for(int i=0;i<DIM;i++)
				{	int	k=DIM*i,ii=iVij+k;
					for(int j=0;j<DIM;j++)
					{	double	v=0.0;
						Vf[k+j]=0.5*(var[ii+j]+var_neib[ii+j]);
					}
				}
			}
#endif
			//Contribution to dVdX:
			for(int i=0;i<DIM;i++)
			{	double	
///					*v=Vf+DIM*i,
///					b=SCLP(v,a);
					b=0.0;
				for(int j=0;j<DIM;j++)
					b+=a[j]*(Vf[DIM*i+j]+Vf[DIM*j+i]);
				b*=0.5;
				d[i]=b;
				dVdX[i]+=b;
			}
			if(facetype!=internal)
///			if(neib==NULL)
			{	//Contribution to the boundary nodes
				for(int jv=1;jv<Nv;jv++)
				{	double	*dVdXv=verts[(iv+jv)%Nv]->var+idVijdXj;
					for(int i=0;i<DIM;i++)
						dVdXv[i]+=d[i];
				}
			}
		}
		cell=cell->next;
	}	while(cell!=dcell_root);
	node=dnode_root;
	do//Normalize by volume to get the derivatives right
	{	double
			*var=node->var,
			vol=var[ivoln],voli=0.5/vol,//coeff 0.5 accounts for 2*area
			*dVdX=var+idVijdXj;
		for(int i=0;i<DIM;i++)
			dVdX[i]*=voli;
		node=node->next;
	}	while(node!=dnode_root);
}
void	Domain::dVidXjXj
(	//Laplace operator: diffusion of a vector
	int	ivoln,
	int	iVi,   // index to the variable to be differentiated
	int	idVidXjXj
)
{
	DNode	*node;
	DCell	*cell;
	// Set dVidXjXj to zero:
	if(dnode_root==NULL)return;
	node=dnode_root;
	do
	{	double	*var=node->var+idVidXjXj;
		for(int i=0;i<DIM;i++)
			var[i]=0.0;
		node=node->next;
	}	while(node!=dnode_root);
	cell=dcell_root;
	do//Assemble the derivatives
	{	DNode	**verts=cell->vert;
		DCell	**neibs=(DCell**)cell->neib;
		ElementStatus	*facetypes=cell->facetype;
		double	e[Nv][DIM];//tetrahedral edge-vectors
		EDGES(e,verts);//computing edge-vectors
		for(int iv=0;iv<Nv;iv++)
		{	DNode	*vert=verts[iv];
			DCell	*neib=neibs[iv];
			double
				*xv=vert->x,//coordinates of the vertex iv
				a[DIM],dd,b,///DDD d[DIM],
				*var=verts[iv]->var,
				*vv=var+iVi,//value of Vi at the vertex
				*dV=var+idVidXjXj,//value of Vi_jj
				dv[DIM];//First derivative of Vi at the face
			ElementStatus	facetype=facetypes[iv];
			AREA2(iv,e,a);//computing face area vector (actually 2*area)
			if(facetype!=internal)
///			if(neib==NULL)
			{	//Use backward derivative approximation
				//Compute face center coordinates and variable 
				double	xf[DIM],vf[DIM],
					bnd=facetype==inlet&&vert->state.boundary==0?1.0:0.0;//DDD
				for (int i=0; i<DIM; i++)
				{	int	ii=iVi+i;
					double	xcenter=0.0,vcenter=0.0;
					for (int jv=1; jv<Nv; jv++)
					{	DNode	*vertj=verts[(iv+jv)%Nv];
						xcenter+=vertj->x[i];
						vcenter+=vertj->var[ii];
					}
					xf[i]=xcenter/(double)Nfv;
					vf[i]=vcenter/(double)Nfv;
					dv[i]=bnd*(vf[i]-vv[i]);
				}
			}
			else//neib!=NULL
			{	//Use central differencing for the derivative
				int	jv=(int)(mask[iv]&cell->index)>>((uint)iv<<1);
					//jv is the index of the opposite node to the face iv
					// as addressed in the neighbor cell
				DNode	*neibvert=neibs[iv]->vert[jv];
				double	*xn=neibvert->x,
					*var_neib=neibvert->var,
					*vn=var_neib+iVi;
				int	nodetype=vert->type;
				dd=0.0;
				for (int i=0; i<DIM; i++)
				{	double	r=xn[i]-xv[i];
					dd+=r*r;
					dv[i]=vn[i]-vv[i];
				}
			}
b=LENGTH(a)/sqrt(dd);///DDD
			for(int i=0;i<DIM;i++)
			{///DDD	d[i]=b*dv[i];
				dV[i]+=b*dv[i];///DDD d[i];
			}
///DDD				if(neib==NULL)
///DDD				{	//Contribution to the boundary nodes
///DDD					for(int jv=1;jv<Nv;jv++)
///DDD					{	double	*dVj=verts[(iv+jv)%Nv]->var+idVidXjXj;
///DDD						for(int i=0;i<DIM;i++)
///DDD							dVj[i]+=d[i];
///DDD					}
///DDD				}
		}
		cell=cell->next;
	}	while(cell!=dcell_root);
	node=dnode_root;
	do//Normalize by volume to get the derivatives right
	{	double
			*var=node->var,
			vol=var[ivoln],voli=0.5/vol,//coeff 0.5 accounts for 2*area
			*dVdX=var+idVidXjXj;
///if(node->type<0)///DDD_2
///for(int i=0;i<DIM;i++)
///dVdX[i]=0.0;
///else
		for(int i=0;i<DIM;i++)
			dVdX[i]*=voli;
		node=node->next;
	}	while(node!=dnode_root);
}
void	Domain::dViUjdXj
(	// Differentiation of ViUj over Xj: Nabla_j.(ViUj)
	// Vi is approximated using the upwind scheme
	int	ivoln, // memory offset of a nodal volume variable
	int	iVi,   // Convected vector
	int iUi,   // Velocity
	int	idViUjdXj
)
{	// Convection of a vector:
	// Input: two vectors; Ouput: vector
const double onethird=1./3.;///DDD
double	vol=0.0,aout[3]={0.0,0.0,0.0};///DDD
	DNode	*node;
	DCell	*cell;
	// Set dVdXi to zero:
	if(dnode_root==NULL)return;
	node=dnode_root;
	do
	{	double	*var=node->var+idViUjdXj;
		for(int i=0;i<DIM;i++)
			var[i]=0.0;
		node=node->next;
	}	while(node!=dnode_root);
	//Loop through the cells and 
	// assemble derivatives
	if(dcell_root==NULL)return;
	cell=dcell_root;
	do
	{	DNode	**verts=cell->vert;
		DCell	**neibs=(DCell**)cell->neib;
		ElementStatus	*facetypes=cell->facetype;
		double	e[Nv][DIM];//tetrahedral edge-vectors
		EDGES(e,verts);//computing edge-vectors
		for(int iv=0;iv<Nv;iv++)
		{	double	a[DIM],Vf[DIM],Uf[DIM],b,d[DIM],
				*var=verts[iv]->var,
				*U=var+iUi,
				*V=var+iVi,
				*dVdX=var+idViUjdXj;
			ElementStatus	facetype=facetypes[iv];
			DNode	*vert=verts[iv];
			DCell	*neib=neibs[iv];
			int	nodetype=vert->type;
			AREA2(iv,e,a);//computing face area vector (actually 2*area)
			if(facetype==internal)
///			if(neib!=NULL)
			{	int	jv=(int)(mask[iv]&cell->index)>>((uint)iv<<1);
				//jv is the index of the opposite node to the face iv
				// as addressed in the neighbor cell
				DNode	*neibvert=neib->vert[jv];


///DCell	**neibneibs=(DCell**)neibs[iv]->neib;
///for(jv=0;jv<Nv;jv++)if(neibneibs[jv]==cell)break;
///if(jv>=Nv){fprintf(stderr,"6:jv=%d\n",jv);exit(1);}
///neibvert=neibs[iv]->vert[jv];


///					if
///					(	  (!vert->state.boundary||nodetype==(int)outlet)
///						&&(!neibvert->state.boundary||neibvert->type==(int)inlet||neibvert->type==(int)outlet)
///	///					&&(!(vert->state.boundary==1&&neibvert->state.boundary==1))
///					)
//				if((vert->state.boundary==0)&&(neibvert->state.boundary==0||neibvert->type==inlet||neibvert->type==outlet))
				{
#ifdef HYBRID
					const	double	hybrid_threshold=0.5;
					double	exit_flux,//=SCLP(Uf,a);
						*var=neibvert->var,
						*Uneib=var+iUi,
						*Vneib=var+iVi,
						Uc[DIM],
						exit_vel,threshold_vel;
					for(int i=0;i<DIM;i++)
						Uc[i]=U[i]+Uneib[i];
					exit_flux=SCLP(a,Uc);
///					exit_flux=SCLP(a,Uf);
					exit_vel=exit_flux/LENGTH(a);
					threshold_vel=hybrid_threshold*LENGTH(Uc);
///					threshold_vel=hybrid_threshold*LENGTH(Uf);
					if(exit_vel>threshold_vel)
					{	for(int i=0;i<DIM;i++)
						{	Uf[i]=U[i];
							Vf[i]=V[i];
						}
					}
					else
					if(exit_vel<-threshold_vel)
					{
						for(int i=0;i<DIM;i++)
						{	Uf[i]=Uneib[i];
							Vf[i]=Vneib[i];
						}
					}
#ifdef CENTRAL_DIFFERENCING
					else
					{	for(int i=0;i<DIM;i++)
						{	Uf[i]=0.5*Uc[i];
							Vf[i]=0.5*(V[i]+Vneib[i]);
						}
					}
#endif//CENTRAL_DIFFERENCING
#else //UPWIND
#ifdef UPWIND
					if(exit_flux>0.0)
						for(int i=0;i<DIM;i++)
						{	Uf[i]=U[i];
							Vf[i]=V[i];
						}
					else
					{	double	*var=neibvert->var,
							*Uneib=var+iUi,
							*Vneib=var+iVi;
						for(int i=0;i<DIM;i++)
						{	Uf[i]=Uneib[i];
							Vf[i]=Vneib[i];
						}
					}
#endif///UPWIND
#endif///UPWIND/HYBRID
				}
			}
			else// neib==NULL
			{
				if
				(	(vert->state.boundary&&nodetype!=outlet)||
					(facetype!=(int)outlet&&facetype!=(int)inlet)
				)
				{	for(int i=0;i<DIM;i++)
						Uf[i]=0.0;
				}
				else
				{	//Compute U and V at the center of the face as face average
					for(int i=0;i<DIM;i++)
					{	double	u=0.0,v=0.0;
						for(int jv=1;jv<Nv;jv++)
						{	double	*var=verts[(iv+jv)%Nv]->var;
							u+=var[iUi+i];
							v+=var[iVi+i];
						}
						Uf[i]=u/Nfv;
						Vf[i]=v/Nfv;
					}
				}
			}
			b=SCLP(Uf,a);
			//Contribution to dVdX:
			for(int i=0;i<DIM;i++)
			{	d[i]=Vf[i]*b;
				dVdX[i]+=d[i];
			}	
///			if(neib==NULL)
			if(facetype!=internal)
			{	//Contribution to the boundary nodes
				for(int jv=1;jv<Nv;jv++)
				{	DNode	*nbvert=verts[(iv+jv)%Nv];
					double	*dVdXv=nbvert->var+idViUjdXj;
					if(nbvert->type==(int)outlet)
					for(int i=0;i<DIM;i++)
						dVdXv[i]+=d[i];
				}
			}
		}
		cell=cell->next;
	}	while(cell!=dcell_root);
	node=dnode_root;
	do//Normalize by volume to get the derivatives right
	{	double
			*var=node->var,
			vol=var[ivoln],voli=0.5/vol,//coeff 0.5 accounts for 2*area
			*dVdX=var+idViUjdXj;
		for(int i=0;i<DIM;i++)
			dVdX[i]*=voli;
		node=node->next;
	}	while(node!=dnode_root);
}
