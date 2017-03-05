#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "main.h"
#include "io.h"
#include "geom.h"
#include "func.h"
#include "var.h"
#include "vecalg.h"
//#include "particles.h"
#include "tool.h"
#include "domain.h"

int	ConnectivitySize[DIM+1][DIM+1]=
{
	//node  edge  face  cell
	{ -1  , -1  , -1  , -1  }, //node
	{  2  , -1  , -1  , -1  }, //edge
	{  5  ,  3  ,  6  ,  2  }, //face
	{  4  ,  4  ,  4  ,  4  }  //cell
};
char
	*elementype[]=
	{
		(char*)"points",
		(char*)"nodes",
		(char*)"edges",
		(char*)"faces",
		(char*)"cells",
		(char*)"boundary nodes",
		(char*)"boundary faces",
		(char*)"boundary cells"
	},
	*spantype[]=  {(char*)"element",(char*)"neighbors",(char*)"grid"};

//?double	nullvector[DIM]={0.0,0.0,0.0};
//?void	Domain::sort()
//?{	int
//?		nf, //total number of faces
//?		nif,//number of internal (fluid) faces	
//?		nc=Ne[cells],
//?		i,*I;
//?//printf("DEBUG:sort skipped\n");return;//DEBUG
//?	if(option.debug)
//?	{	printf("Sorting all cells: from boundary inward\n");FLUSH;
//?		printf
//?		(	"Allocating temporal data storage of %d*%d=%d bytes ...",
//?			nc,sizeof(int),nc*sizeof(int)
//?		);FLUSH;
//?	}
//?	ALLOC(I,int,nc);//cell-flags: 0-not scanned, 1-scanned
//?	if(option.debug){printf(" done\nInitializing ...");FLUSH;}
//?	for (i=0; i<nc; i++)I[i]=0;
//?	//Count number of boundary faces and
//?	//Collect all boundary cells at the beginning of cell-array
//?	if(option.debug){printf(" done\nSorting boundary cells\n");FLUSH;}
//?	nbf=nbc=0;
//?	for (i=0; i<nc; i++)
//?	{	//Skip first boundary cells if any
//?		int	boundarycell=0;
//?		for (int j=0; j<Nf; j++)
//?			if (cell[i].cell[j]<0)
//?			{	nbf++;boundarycell=1;
//?			}
//?		if (!boundarycell) break;
//?		nbc++;
//?		I[i]=1;
//?	}
//?	if (i<nc-1)
//?	for (int j=i+1;j<nc;j++)
//?	{	int	boundarycell;
//?		boundarycell=0;
//?		for (int k=0; k<Nf; k++)
//?			if (cell[j].cell[k]<0)
//?			{	nbf++;boundarycell=1;
//?			}
//?		if (boundarycell)
//?		{	nbc++;
//?			swapCells(i,j);
//?			I[i++]=1;
//?		}
//?	}
//?	if(option.debug)
//?	{	printf("Number of boundary cells=%d\n",nbc);
//?		printf("Number of faces: boundary=%d,",nbf);
//?	}
//?	if((nif=4*nc-nbf)%2!=0)
//?	{	fprintf
//?		(	stderr,
//?			"ERROR: 4 * no_of_cells (%d) - no. boundary faces (%d) = %d is not a multiple of 2\n",
//?			nc,nbf,nif
//?		);
//?		exit(1);
//?	}
//?	nif/=2;//number of internal faces
//?	nf=nif+nbf;
//?	if (Ne[faces]<=0)
//?		Ne[faces]=nf;//Total number of faces
//?	else
//?		if (Ne[faces]!=nf)
//?		{	fprintf
//?			(	stderr,
//?				"Sort: Computed number of faces (%d) does not agree with the earlier value (%d)",
//?				nf,Ne[faces]
//?			);
//?		}
//?	if(option.debug)
//?	{	printf("internal=%d, all=%d\n",nif,Ne[faces]);
//?		printf("Sorting all cells ...",nif,Ne[faces]);FLUSH;
//?	}
//?	//Sort all cells from the boundary inward
//?	// suits well for an elliptic solver
//?	// for a parabolic solver only inlet boundary - invard sort
//?	// should be considered
//?	i=0;
//?	for (int j=nbf; i<j&&j<nc; i++)
//?	{	//move the neighbors of cell i towards the beginning
//?		// of cell-array 
//?		for (int k=0; k<Nf; k++)
//?		{	int	n=cell[i].cell[k];
//?			if (n>=0&&I[n]==0)
//?			{	swapCells(j,n);
//?				I[j]=1;
//?				if (++j==nc)break;
//?			}
//?		}
//?	}
//?	if(option.debug)
//?	{	printf(" done.\n",nif,Ne[faces]);
//?		printf
//?		(	"Releasing temporal data storage of %d*%d=%d bytes ...",
//?			nc,sizeof(int),nc*sizeof(int)
//?		);FLUSH;
//?	}
//?	free(I);
//?}
//?void	Domain::setFaces()
//?{	/* Face-cell connectivity */
//?	int	iface,
//?		nf=Ne[faces],
//?		nc=Ne[cells];
//?	//Count number of boundary faces and
//?	//Collect all boundary cells at the beginning of cell-array
//?	if (type==dynamic)
//?	{	printf("WARNING: setFaces not implemented for domain.type=dynamic\n");
//?		fflush(stdout);
//?		return;
//?	}
//?	if (option.debug){printf("setFaces\n");FLUSH;}	
//?	sort();
//?	if(option.debug)
//?	{	printf
//?		(	"Allocating %d*%d=%d bytes for faces\n",
//?			nf,sizeof(struct Face),
//?			nf*sizeof(struct Face)
//?		);FLUSH;
//?	}
//?	ALLOC(face,struct Face,Ne[faces]);
//?	if(option.debug){printf("Constructing %d boundary faces ...",nbf);FLUSH;}
//?///{//count boundary cells
//?///	 int nbc=0;
//?///printf("Checking %d cells:\n",nc);
//?///for (int ic=0; ic<nc; ic++)
//?///for (int iv=0; iv<Nv; iv++)
//?///if (cell[ic].cell[iv]<0)nbc++;
//?///printf("Number of boundary cells: nbc=%d\n",nbc);
//?///}
//?	for (int ic=0; ic<nc; ic++)
//?		for (int i=0; i<Nf; i++) cell[ic].face[i]=-1;
//?		//all cell[].face[] should be >=0 in the end
//?	//Setup the boundary faces
//?	iface=0;
//?	for (int ic=0; ic<nbc; ic++)
//?	{//Boundary cells first,
//?	 // so that boundary faces are all
//?	 // aligned at the beginning of face-array
//?		struct Cell *c=cell+ic;
//?		int	*v=c->vert,*cc=c->cell;
//?		double
//?			e[Nv][DIM],//edges of a tetrahedron
//?			*X=coordinates[nodes].val,
//?			*Y=coordinates[faces].val,
//?			*x=X+DIM*v[0];
//?		for (int i=0; i<Nf; i++)
//?		{	double
//?				*y=X+DIM*v[(i+1)%Nv];
//?			for (int k=0; k<DIM; k++)
//?				e[i][k]=y[k]-x[k];
//?			x=y;
//?		}
//?		c->mirror.inside=0;
//?		for (int i=0; i<Nf; i++)
//?		{	int
//?				jc=cc[i], //neighbor cell index
//?				i1=(i+1)%Nv,
//?				i2=(i+2)%Nv;
//?			double a[DIM],area,areai,//face area vector
//?				b,*y,*z;
//?			struct Face	*f;
//?			if (jc>=0) continue;//Internal faces are setup later
//?			if (iface>=nbf) BUG("Grid inconsistency in the number of faces");
//?			f=face+iface;
//?			y=Y+DIM*iface;
//?			f->cell[0]=ic; //0-cell points to ic
//?			f->cf=i;//local index of this face at cell 0
//?			f->cell[1]=jc; //1-cell points to jc => boundary conditions
//?			//Face-vertex connectivity:
//?			for (int j=0; j<Nfv; j++)
//?				f->vert[j]=v[(i+j+1)%Nv];
//?			//Face-center coordinates:
//?			for (int k=0; k<DIM; k++)
//?			{	b=0.0;
//?				for (int j=0; j<Nfv; j++)
//?					b+=X[DIM*f->vert[j]+k];
//?				y[k]=b/(double)Nfv;
//?			}
//?			//Compute face area vector
//?			VECP(a,e[i1],e[i2]);
//?			b = SCLP(a,e[i])<0 ? -0.5 : 0.5;
//?			for (int j=0; j<DIM; j++) a[j]*=b;
//?			area=LENGTH(a);
//?			areai=1./area;
//?			f->area=area;
//?			for (int j=0; j<DIM; j++) f->norm[j]=a[j]*areai;
//?			/* now f->norm is the face area normal vector 
//?			 * directed away from the 0-cell
//?			 */
//?			if (iface>=nbf) BUG("Face numbers don't match");
//?			c->face[i]=iface++;//cell-face connectivity
//?			//c->mirror.inside+=0<<i;
//?		}
//?	}
//?	coordinates[boundary_faces].val=coordinates[faces].val;
//?	if(option.debug){printf(" done\n",iface);FLUSH;}
//?	//Setup the internal faces
//?	//ALLOC(I,struct CellFace,nc);
//?	if(option.debug){printf("Computing internal faces ...");FLUSH;}
//?	for (int ic=0; ic<nc; ic++)
//?	{	struct Cell *c=cell+ic;
//?		int	*v=c->vert,*cc=c->cell;
//?		double
//?			d[Nv][DIM],//edges of a tetrahedron
//?			*X=coordinates[nodes].val,
//?			*Y=coordinates[faces].val,
//?			*x=X+DIM*v[0];
//?		for (int i=0; i<Nf; i++)
//?		{	double
//?				*y=X+DIM*v[(i+1)%Nv];
//?			for (int k=0; k<DIM; k++)
//?				d[i][k]=y[k]-x[k];
//?			x=y;
//?		}
//?		c->mirror.vv=0;
//?		c->mirror.inside=0;
//?		for (int i=0; i<Nf; i++)
//?		{	int	jf,//the local neighbor-face index: jf=0:Nf
//?				      // the index to the current face i as reffered
//?				      // from the neighbor cell jc
//?				jface,//the global face index of the neighbor cell: jface=0:Ne[faces]
//?				jc=cc[i], //neighbor cell index
//?				i1,i2;
//?			double a[DIM],area,areai,g,*y;//face area vector
//?			struct Cell	*b;//neighbor cell
//?			struct Face	*f;
//?			if (jc<0) continue;//boundary faces were processed earlier.
//?			i1=(i+1)%Nv;
//?			i2=(i+2)%Nv;
//?			b=cell+jc;
//?			//Retrieving the local face index (j=0...Nf-1)
//?			//of face i as addressed on the neighbor cell
//?			/* Set neighbor face indexes: c->f */
//?			for (jf=0; jf<Nf&&b->cell[jf]!=ic; jf++);
//?//-					for (j=0; j<Nv; j++) 
//?//-					{//algorithm valid only for a tetrahedron
//?//-						int	k;
//?//-						for (k=0; k<Nv&&b->vert[j]!=vert[k]; k++);
//?//-						if (k==Nv) break;
//?//-					}
//?			//Check if the face is already registered
//?			if (jf<Nv)
//?			{	//set the lower to bits to vertex number
//?				double
//?					a1,a2,e[DIM],*x;
//?//					c->mirror.vv=mask[i]&j<<(i<<1);
//?				c->mirror.vv+=jf<<(i<<1);
//?				// The neigbor face index is now stored in vv
//?				// and can be retrieved for every face i as
//?				//j=(int)(mask[i]&c->vv)>>((uint)i<<1);
//?				if ((jface=b->face[jf])>=0) 
//?				{	//int	jface=cell[jc].face[j];
//?					//Face was already visited from cell jc
//?					//refer to it with a negative index:
//?//-					c->face[i]=-jface-1;//which designate that 
//?					c->face[i]=jface;//which designate that 
//?												// area vectors point inside the cell
//?					c->mirror.inside+=1<<i;
//?					continue;
//?				}
//?				if (iface>=nf) BUG("Grid inconsistency in the number of faces");
//?				f=face+iface;
//?				y=Y+DIM*iface;
//?				//Fece-node connectivity
//?				for (int k=0; k<Nfv; k++)
//?					f->vert[k]=v[(i+k+1)%Nv];
//?				//Face-center coordinates
//?				for (int k=0; k<DIM; k++)
//?				{	g=0.0;
//?					for (int j=0; j<Nfv; j++)
//?					g+=X[DIM*f->vert[j]+k];
//?//					g+=X[DIM*v[(i+j+1)%Nv]+k];
//?					y[k]=g/(double)Nfv;
//?				}
//?				//face area vector
//?				VECP(a,d[i1],d[i2]);
//?				g = SCLP(a,d[i])<0 ? -.5 : .5;
//?//			if (SCLP(a,e[i])<0)
//?				for (int k=0; k<DIM; k++) a[k]*=g;
//?				area=LENGTH(a);
//?				areai=1./area;
//?				f->area=area;
//?				for (int k=0; k<DIM; k++)
//?					f->norm[k]=a[k]*areai;
//?				//	for (double *p=a; p-a<DIM; p++) *p*= b;
//?				/* now f->norm is the face area unit normal vector 
//?				 * and a is the face area normal vector
//?				 * directed away from the 0-cell
//?				 */
//?				//Face-cell connectivity:
//?				f->cell[0]=ic; //0-cell points to ic
//?				f->cell[1]=jc; //1-cell points to jc
//?				f->cf=jf*Nf+i; //local indexes of this face at neighbor cells
//?				c->face[i]=iface++;//cell-face connectivity
//?			  //Check for overlapping cells
//?//-				z=node[b->vert[j]].x;//neighbor's vertex opposite to the face
//?				x=X+DIM*b->vert[jf];//neighbor's vertex opposite to the face
//?					//Find normal vector to the face i
//?//-				x=f->x;//face-center vector
//?				for (int  k=0; k<DIM; k++)
//?				{  //coordinates of the vertex facing face i
//?					//g[k]=y[k]-x[k];=d[2]
//?				  //coordinates of the opposite vertex to face i
//?					/* x[0:2]= coordinates of any vertex on face i 
//?						(set in the beginning of the i-loop) */
//?//-					e[k]=z[k]-x[k];
//?					e[k]=x[k]-y[k];
//?				}
//?				a1=(SCLP(d[i],a));
//?				a2=(SCLP(e,a));
//?				if (a1*a2<0.0)
//?				{	fprintf
//?					(	stderr,
//?						"\nError: The grid has overlapping cells at: %g, %g, %g\n",
//?						x[0],x[1],x[2]
//?					);
//?					exit(1);
//?				}
//?			}
//?			else
//?			{	int	*cv=c->vert, *bv=b->vert;
//?				fprintf
//?				(	stderr,
//?					"GRID ERROR: Neighbor cells %d,%d don't point to each other:\n"
//?					"\t%d,%d,%d,%d : %d,%d,%d,%d:\n",ic,jc,
//?					cv[0],cv[1],cv[2],cv[3],
//?					bv[0],bv[1],bv[2],bv[3]
//?				);
//?				exit(1);
//?			}
//?		}
//?	}
//?	if (option.debug)printf(" %d faces done\n",iface);FLUSH;
//?#ifdef DEBUG
//?	for (int iface=0; iface<nbf; iface++)
//?	if (face[iface].cell[1]>=0) 
//?		BUG("setFaces: non-boundary face in the wrong part of the array");
//?#endif
//?}
void	Domain::setBoundaryCellFlags()
{
	DCell	*cell=dcell_root;
	if (dcell_root==NULL)return;
	do
	{	int	i;
		DCell	**neibs=(DCell**)cell->neib;
		ElementStatus	*facetype=cell->facetype;
		for (i=0; i<Nv; i++)
		{	if(neibs[i]==NULL||facetype[i]!=internal)
			{	cell->state.boundary=1;
				break;
			}
		}
		if (i==Nv) cell->state.boundary=0;
		cell=cell->next;
	}	while (cell!=dcell_root);
}
void	Domain::setBoundaryVertexes()
//For static domains should be called after sort (nbc must be known)
{	int	
		nn=Ne[nodes],
		ibv,*I;
//?	if (type==dynamic)
//?	{	
		DNode	*node_root=dnode_root,*node=node_root;
		DCell	*cell_root=dcell_root,*cell=cell_root;
		if (node_root==NULL)return;
		do //SET THE BOUNDARY FLAGS TO ZERO
		{	node->state.boundary=0;
			node=node->next;
		}	while(node!=node_root);
		if (cell_root==NULL)return;
		cell=cell_root;
		do //SET THE BOUNDARY FLAGS
		{	DNode	**verts=cell->vert;
			void	**neibs=cell->neib;
			int 	*facetypes=(int*)cell->facetype;
			for (int iv=0; iv<Nv; iv++)
			{	if(neibs[iv]==NULL)
				{	int	facetype=(int)facetypes[iv];
					if(facetype<0)//this is a real boundary 
					//	- not the outlet into a parallel domain
					for(int jv=0;jv<Nv1;jv++)
					{	int	jvert=(iv+jv+1)%Nv;
						DNode	*node=verts[jvert];
						node->state.boundary=1;
						if(facetype<boundary)//in/out-let
							node->state.fixed=1;
						if(node->type==0||node->type<facetype)
							node->type=facetype;
					}
				}
			}
			cell=cell->next;
		}	while(cell!=cell_root);
		return;
//?	}
//?	//Static domain
//?	if(option.debug) 
//?		printf("Allocating %d*%d bytes of temporal storage\n",nn,sizeof(int));
//?	ALLOC0(I,int,nn);
//?	//Count boundary vertexes
//?	ibv=0;
//?	for (int ic=0; ic<nbc; ic++)
//?	{	struct Cell	*c=cell+ic;
//?		int
//?			*cv=c->vert,
//?			*cc=c->cell;
//?		for (int jf=0; jf<Nf; jf++)//scan the faces
//?		if (cc[jf]<0)
//?		{	//Boundary face
//?			for (int iv=0; iv<Nfv; iv++)
//?			{	int	ivert=cv[(iv+jf+1)%Nv];
//?				if (I[ivert]==1) continue;
//?				I[ivert]=1;
//?				ibv++;
//?			}
//?		}
//?	}
//?	nbv=ibv;
//?	Ne[boundary_nodes]=nbv;
//?	if (option.verbose|option.debug)
//?	{	printf("Number of boundary vertexes: nbv=%d\n",nbv);FLUSH;
//?	}
//?	if (option.debug)
//?	{	printf
//?		(	"Initializing %d*%d=%d bytes for boundary vertex data\n",
//?			nbv,sizeof(struct BoundaryVertex),
//?			nbv*sizeof(struct BoundaryVertex)
//?		);FLUSH;
//?	}
//?	ALLOC(Bv,struct BoundaryVertex,nbv);
//?	//Assign BoundaryVertex struct
//?	if(option.debug) 
//?		printf
//?		(	"Assigning conditions on %d boundary vertexes from %d boundary cells\n",
//?			nbv,nbc
//?		);	
//?	for (int i=0; i<nbv; i++) Bv[i].b=0;
//?	for (int i=0; i<nn; i++) I[i]=0;
//?	ibv=0;
//?	for (int ic=0; ic<nbc; ic++)
//?	{	struct Cell	*c=cell+ic;
//?		int
//?			*cv=c->vert,
//?			*cf=c->face,	
//?			*cc=c->cell;
//?		for (int jf=0; jf<Nf; jf++)//scan the faces
//?		if (cc[jf]<0)
//?		{	//Boundary face
//?			int	face_boundary_condition=cc[jf];
//?			for (int iv=0; iv<Nfv&&ibv<nbv; iv++)
//?			{	int	ivert=cv[(iv+jf+1)%Nv],
//?					nodetype=node[ivert].type;
//?				if (I[ivert]==1) continue;
//?				if (nodetype!=internal)
//?				{	if(face_boundary_condition<nodetype)
//?					{
//?						Bv[ibv].b=face_boundary_condition;
//?						node[ivert].type=face_boundary_condition;
//?					}
//?					else
//?						Bv[ibv].b=nodetype;
//?				}
//?				I[ivert]=1;
//?				Bv[ibv].i=ivert;
//?				node[ivert].ind=ibv;//connect node to the boundary 
//?				for (int k=0; k<DIM; k++) Bv[ibv].norm[k]=0.0;
//?				ibv++;
//?				if(ibv>nbv)BUG("Domain::setBoundaryVertexes: mismatch in the number of boundary vertexes");
//?			}
//?		}
//?	}
//?#ifdef DEBUG
//?	//Check 1
//?	for (int i=0; i<nbv; i++)
//?	if (node[Bv[i].i].ind!=i) printf("ERROR: vertex: %d, node=%d\n",i, node[Bv[i].i].ind);
//?	//Check 2
//?	for (int i=0; i<nbv; i++) I[i]=0;
//?	for (int iface=0; iface<nbf; iface++)
//?	{//compute boundary normal vectors at the boundary nodes
//?		struct Face	*f=face+iface;
//?		int	*fv=f->vert;
//?		for (int iv=0; iv<Nfv; iv++)
//?		{	int	inode=node[fv[iv]].ind;
//?			I[inode]=1;
//?			if (inode>=nbv)
//?			{	printf("inode=%d >= nbv=%d; iface=%d, fv[%d]=%d\n",inode,nbv,iface,iv,fv[iv]);FLUSH;//DEBUG
//?			}
//?		}
//?	}
//?	for (int i=0; i<nbv; i++)
//?	if (I[i]==0) 
//?	{	int	ibv=Bv[i].i;
//?		double	*x=coordinates[nodes].val+DIM*ibv;
//?		printf("Missed point: %d: %g,%g,%g\n",i,x[0],x[1],x[2]);FLUSH;
//?	}
//?	//Check 3
//?
//?#endif
//?	if(option.debug||option.verbose)
//?	{	printf
//?		(	"done\nComputing boundary normal vectors by scanning %d face elements ... ", 
//?			 nbf
//?		);FLUSH;
//?	}
//?	for (int iface=0; iface<nbf; iface++)
//?	{//compute boundary normal vectors at the boundary nodes
//?		struct Face	*f=face+iface;
//?		int
//?			*fv=f->vert;
//?		double
//?			*norm=f->norm;
//?		for (int iv=0; iv<Nfv; iv++)
//?		{	int	inode=node[fv[iv]].ind;
//?			for (int k=0; k<DIM; k++)
//?				Bv[inode].norm[k]+=norm[k];
//?			I[inode]++;
//?		}
//?	}
//?	for (int ibv=0; ibv<nbv; ibv++)
//?	{	//Nomalize boundary norm vectors
//?		if (I[ibv]<=0)
//?		{	fprintf(stderr,"I[%d]=%d\n",ibv,I[ibv]);
//?			BUG("DIVISION BY ZERO IN Domain::setBoundaryVertexes");
//?		}
//?		for (int i=0; i<DIM; i++) Bv[ibv].norm[i]/=(double)I[ibv];
//?	}
//?	if (option.debug)
//?	{	//Check unity of normals
//?		for (int ibv=0; ibv<nbv; ibv++)
//?		{	double	*norm=Bv[ibv].norm;
//?			if (sqrt(norm[0]*norm[0]+norm[1]*norm[1]+norm[2]*norm[2])<SMALL)
//?				BUG("UNITY OF NORMALS VIOLATED IN Domain::setBoundaryVertexes");
//?		}
//?		printf(" done\n");FLUSH;
//?	}
//?	if (option.debug) printf("Deallocating temp storage\n");
//?	free (I);
//?}
//?void	Domain::setVolumes()
//?{
//?	static const double	onethird=1./3.;
//?	int
//?		idomain=this-domain_root,
//?		nn=Ne[nodes],
//?		nf=Ne[faces],
//?		nc=Ne[cells];
//?	double
//?		*X=coordinates[nodes].val,
//?		*Y=coordinates[faces].val,
//?		*Z=coordinates[cells].val;
//?if(type==dynamic)return;//DEBUG
//?	//Compute cell-volumes (as in CV cell-centered scheme)
//?	fluxFaceVec2CellScl
//?	(	Ne[cells],//number of cells
//?		cell,face,
//?		coordinates[faces].val,//face-vector variable which gradient is to be computed
//?		volumes[cells].val //scalar node-variable representing the delatation of V
//?	);
//?	if(option.debug)//Test the total volume:
//?	{	int	nc=Ne[cells];
//?		double
//?			*volc=volumes[cells].val,
//?			totvol=0.0;
//?		for (double *v=volc; v-volc<nc; v++)totvol+=*v;
//?		printf
//?		(	"Total volume of domain %d computed by cell scheme = %g\n",
//?			idomain,totvol
//?		);
//?	}
//?	//Compute nodal volumes (as in FEM overlapping CV schame)
//?	fluxFaceVec2NodeScl
//?	//Computes delatation of a face-vector at the nodes
//?	// delatation=(\nabla\,,variable)
//?	(	Ne[cells],//number of cells
//?		Ne[nodes],//number of nodes
//?		cell,face,
//?		coordinates[faces].val,//face-vector variable which gradient is to be computed
//?		volumes[nodes].val //scalar node-variable representing the delatation of V
//?	);
//?	if(option.debug)//Test the total volume:
//?	{	int	nn=Ne[nodes];
//?		double
//?			*voln=volumes[nodes].val,
//?			totvol=0.0;
//?		for (double *v=voln; v-voln<nn; v++)totvol+=*v;
//?		printf
//?		(	"Total volume of domain %d computed by nodal scheme = %g\n",
//?			idomain,0.25*totvol
//?		);
//?	}
//?	//Compute cell-center coordinates
//?	for (int ic=0; ic<nc; ic++)
//?	{	int	*cv=cell[ic].vert;
//?		double	*z=Z+DIM*ic;//Cell-center coordinates
//?		for (int i=0; i<DIM; i++)
//?		{	double a=0.0;
//?			for (int j=0; j<Nv; j++)
//?				a+=X[DIM*cv[j]+i];
//?			z[i]=a/(double)Nv;
//?		}
//?	}	
//?	//Compute face-based volumes and areas
//?	{	double
//?			*vol=volumes[faces].val,
//?			totvol=0.0;//DEBUG
//?		for (int iface=0; iface<nf; iface++)
//?			vol[iface]=0.0;
//?		for (int icell=0; icell<nc; icell++)
//?		{//Looping over all cells
//?			struct Cell	*c=cell+icell;
//?			int
//?				*cv=c->vert,//cell-vertex connectivity
//?				*cf=c->face;//cell-face connectivity
//?			double
//?				*z=Z+DIM*icell;//cell-center coordinates
//?			for (int jf=0; jf<Nf; jf++)
//?			{//Looping over cell-faces
//?				int	jface=cf[jf],//global face index
//?					*fv,//face-vertex connectivity array
//?					jdirection=1-(((c->mirror.inside>>jf)%2)<<1);
//?					//1: norm pointed outside; -1: pointed inside
//?				double	dvol,
//?					e0[DIM],//one of the two edges defining the int.face
//?					*norm,//face normal vector
//?					area,//face area
//?					*y,//face-center vector
//?					*x0;//first vertex of the face
//?				struct Face	*f;
//?				//if(jface<0){jface=-jface-1;jdirection=-1;}
//?				//else	jdirection=1;
//?				f=face+jface;
//?				y=Y+DIM*jface;
//?				norm=f->norm;//face normal vector
//?				area=f->area;//face area
//?				fv=f->vert;//face-vertexes
//?				x0=X+DIM*fv[1];//cv[(jf+2)%Nf];//coords of the first vertex of the int.face
//?				for (int i=0; i<DIM; i++)
//?					e0[i]=x0[i]-z[i];//first edge of int.face
//?				dvol=onethird*(double)jdirection*SCLP(norm,y)*area;
//?				for (int iv=0; iv<Nfv; iv++)
//?				{//Looping over face vertexes
//?				 //to find intermediate faces
//?					int	idirection;
//?					double
//?						a[DIM],//int.face area vector
//?						d[DIM],//distance between the int.face center and face jface
//?						yy[DIM],//int.face-center vector
//?						e[DIM],//second edge of the face
//?						Rfv=1./(double)Nfv,
//?						*x=X+DIM*fv[(iv+2)%Nfv];//cv[(jf+iv+3)%Nf];//coords of the next vertex of the face
//?					for (int i=0; i<DIM; i++)
//?					{	e[i]=x[i]-z[i];//second edge of the int.face
//?						yy[i]=Rfv*(x[i]+x0[i]+z[i]);//middle point coords of the int.face
//?						d[i]=yy[i]-y[i];
//?					}
//?					VECP(a,e0,e);//area vector of the int.face
//?					idirection=SCLP(d,a)>=0.0?1:-1;
//?					for (int i=0; i<DIM; i++) a[i]*=0.5*(double)idirection;
//?					//idirection==-1: a points inside; 1: outside
//?					dvol+=onethird*SCLP(a,yy);//DEBUG
//?					for (int i=0; i<DIM; i++)
//?						f->a[(((1-jdirection)>>1)*Nfv+iv)*DIM+i]=a[i];
//?					x0=x;
//?					for (int i=0; i<DIM; i++) e0[i]=e[i];
//?
//?//{//DEBUG
//?//double	r[DIM],b,
//?//*area=f->a+(((1-jdirection)>>1)*Nfv+iv)*DIM;
//?//for(int i=0; i<DIM; i++)
//?//r[i]=y[i]-z[i];
//?////if((b=SCLP(area,r))>0.0)
//?//printf("+icell=%d,jface=%3d(%1d),a[%d]=%g,%g,%g,y=%6.2f%6.2f%6.2f,z=%6.2f%6.2f%6.2f\n",icell,jface,iv,(((1-jdirection)>>1)*Nfv+iv)*DIM,area[0],area[1],area[2],y[0],y[1],y[2],z[0],z[1],z[2]);
//?//}
//?
//?				}
//?				totvol+=dvol;//DEBUG
//?				{	double
//?						d[DIM];//distance from the face-center to the cell-center
//?					for (int i=0; i<DIM; i++)
//?						d[i]=y[i]-z[i];
//?					vol[jface]+=onethird*(double)jdirection*SCLP(d,norm)*area;
//?				}
//?			}
//?		}
//?		if(option.debug)
//?		{	double	v=0.0;
//?			printf("Total face-based volume = %lg\n",totvol);
//?			for (int i=0; i<nf; i++)
//?				v+=vol[i];
//?			printf("Checked face volume = %lg\n",v);
//?			printf("volume1-volume2=%lg\n",v-totvol);
//?			if(fabs(v-totvol)>0.01*totvol/(double)nc)BUG("Volumes don't agree in setVolumes\n");
//?		}
//?	}
//?	//Compute c1ii coefficients: hedra.tex:fem.tex
//?	for (int i=0; i<nn; i++)
//?	{//	node[i].vol=0.0;
//?		node[i].c1ii=0.0;
//?	}
//?	for (int ic=0; ic<nc; ic++)
//?	{	int	*cv,*cf,*cc;
//?		double
//?			*z=Z+DIM*ic,
//?			area[Nf],
//?			vol9i,vol=0.0;//cell-volume,vol9, 
//?		struct Cell *c=cell+ic;
//?		cv=c->vert;
//?		cf=c->face;
//?		cc=c->cell;
//?		for (int i=0; i<Nf; i++)
//?		{
//?			int
//?				jv=cv[i],
//?				jf=cf[i],
//?				jc=cc[i],
//?				direction=1-(((c->mirror.inside>>i)%2)<<1);
//?			double
//?				dVol,
//?				//area,//face areas
//?				*y,
//?				*norm;//face normal vector (direction determined by direction)
//?			struct Face	*f;
//?			/* compute the face area vector a[0..DIM] */
//?			f=face+jf;
//?			y=Y+DIM*jf;//face-center coordinates
//?			norm=f->norm;
//?			area[i]=f->area;
//?			dVol=(double)direction*onethird*area[i]*(SCLP(y,norm));//volume flux
//?			vol+=dVol;
//?		}
//?#ifdef DEBUG
//?		if (vol<SMALL)BUG("Volume too small in setVolumes");
//?#endif
//?		vol9i=1./(9*vol);
//?		for (int i=0; i<Nf; i++)
//?		{
//?			int
//?				jv=cv[i],
//?				jf=cf[i];
//?				//jc=cc[i];
//?			double	dc=area[i]*area[i]*vol9i;
//?			node[cv[i]].c1ii+=dc;
//?			//	if (jc<0)
//?			//	{	//boundary cell
//?			//		for (int j=0; j<Nfv; j++)
//?			//			node[cv[(i+j+1)%Nv]].c1ii+=dc;
//?			//	}
//?		}
//?	}
//?//		if (option.debug)
//?//		{//Compute the total volume
//?//			double	vol=0.0;
//?//			for (int ic=0; ic<nc; ic++) vol+=cell[ic].vol;
//?//			printf("Total volume based on cell-elements: vol=%g\n",vol);
//?//			vol=0.0;
//?//			for (int in=0; in<nn; in++) vol+=node[in].vol;
//?//			printf("Total volume based on node-elements: vol=%g\n",vol/4.);
//?//		}
}
DCell::DCell(void)
{
	for (int i=0; i<Nv; i++)
		vert[i]=NULL;
	for (int i=0; i<Nf; i++)
	{	neib[i]=NULL;
		facetype[i]=boundary;
	}
	next=prev=NULL;
	state.relax=0;
	var=NULL;
}
DCell::DCell(int varbufsize)
{
	for (int i=0; i<Nv; i++)
		vert[i]=NULL;
	for (int i=0; i<Nf; i++)
	{	neib[i]=NULL;
		facetype[i]=boundary;
	}
	next=prev=NULL;
	state.relax=0;
	var=new double[varbufsize];
}
DCell::DCell(DCell *cell)
{
	for (int i=0; i<Nv; i++)
		vert[i]=cell->vert[i];
	for (int i=0; i<Nf; i++)
	{	neib[i]=cell->neib[i];
		facetype[i]=cell->facetype[i];
	}
	next=cell->next;
	prev=cell->prev;
	state.relax=0;
	var=NULL;
}
DCell::DCell(DCell *cell, int varbufsize)
{
	for (int i=0; i<Nv; i++)
		vert[i]=cell->vert[i];
	for (int i=0; i<Nf; i++)
	{	neib[i]=cell->neib[i];
		facetype[i]=cell->facetype[i];
	}
	next=cell->next;
	prev=cell->prev;
	state.relax=0;
	var=new double[varbufsize];
	for (int i=0; i<varbufsize; i++)
		var[i]=cell->var[i];
}
DNode::DNode()
{
	var=NULL;
	bond=NULL;
}
DNode::DNode(int mvar)
{	if(mvar>0)
		var=new double[mvar];
	else
		var=NULL;
	bond=NULL;
}
DNode::DNode(int mvar, double y[])
{
	var=new double[mvar];
	state.relax=0;
	state.insidetool=0;
	state.toolsurface=0;
	state.fixed=0;
	bond=NULL;
	for (int i=0; i<DIM; i++)
		x[i]=y[i];
}
DNode::~DNode()
{
	delete var;
}
DCell::~DCell()
{
	delete var;
}
BFaceList::BFaceList()
{
	var=NULL;
}
BFaceList::BFaceList(int mvar)
{
	var=new double[mvar];
}
BFaceList::~BFaceList()
{
	delete var;
}
void	Domain::checkGrid
(	DNode	*node_root,
	DCell *cell_root
)
{//Grid consistency check:
	int	nc=0,nbf=0;
	printf("Checking grid consistency ...");
	DNode	*node=node_root;
	DCell	*cell=cell_root;
	do
	{	DNode	**verts=cell->vert;
		DCell	**neibs=(DCell**)cell->neib;
		///printf("nc=%d\n",nc);
		for(int iv=0;iv<Nv;iv++)
		{	int	jv;
			DCell	*neib=neibs[iv],
				**neibneibs;
			if(neib!=NULL)
			{	neibneibs=(DCell**)neib->neib;
				///printf("iv=%d\n",iv);FLUSH;
				for(jv=0;jv<Nv;jv++)
				{	///printf("\tj=%d\n",jv);FLUSH;
					if(neibneibs[jv]==cell)break;
				}
				///printf("\tjv=%d\n",jv);FLUSH;
				if(jv==Nv)ERROR("Cell-cell connectivity error");
				{//Cell-proximity check
					//Check cell-node connectivity
					DNode	**neibnodes=neib->vert;
					for(int ivv=1;ivv<Nv;ivv++)
					{	int	jvv;
						DNode	*vert=verts[(iv+ivv)%Nv];
						for(jvv=1;jvv<Nv;jvv++)
						if(neibnodes[(jv+jvv)%Nv]==vert)break;
						if(jvv==Nv)ERROR("Cell-node connectivity error");
					}
				}
			}
			else nbf++;
		}
		cell=cell->next;
		nc++;
	}	while(cell!=cell_root);
	printf(" passed: nc=%d, nbf=%d\n",nc,nbf);
///		node=node_root;
///		do
///		{	double	*x=node->x;
///			if(node->state.boundary)
///			if(fabs(x[0])<2.9&&fabs(x[1])<2.9&&fabs(x[2])<2.9)
///			{
///			printf("neibvert->x=%g,%g,%g\n",x[0],x[1],x[2]);
///			ERROR("BOUNDARY FLAGS");
///			}	
///			node=node->next;
///		}	while(node!=node_root);
}
