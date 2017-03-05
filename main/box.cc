#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "main.h"
#include "io.h"
#include "geom.h"
#include "vecalg.h"
#include "var.h"
//#include "particles.h"
#include "tool.h"
#include "domain.h"
#include "box.h"

#define	J(i,j,k)	(i)+(j)*2+(k)*4

#define	CASE	4
#define	Mv	5

int	Iv[Mv][4]=
//			{1,2,4,7},
//			{0,1,2,4},
//			{1,3,2,7},
//			{1,4,5,7},
//			{4,6,7,2}
	{
		{1,2,4,7},//0
		{0,1,2,4},//1
		{1,2,3,7},//2
		{1,4,5,7},//3
		{2,4,6,7} //4
	};
/*	0		{*,1,2,*,4,*,*,7},something to think about
	1		{0,1,2,*,4,*,*,*},
	2		{*,1,2,3,*,*,*,7},
	3		{*,1,*,*,4,5,*,7},
	4		{*,*,2,*,4,*,6,7}
*/
/*****************************************************************
                           0, 0,+1
      iid=+1:                ^        0,+1,0
                             |        /
                             |       /
                   6_______________________7
                   /                      /|
                  / |                    / |
                 /______________________/  |
                4|  |                   |5 |
 -1, 0, 0 <---   |                      |  |    ---> +1, 0, 0
                 |  |2 _ _ _ _ _ _ _ _ _| _|3
                 |  /                   |  /
                 | /                    | /
                 |______________________|/
                0                        1
              T0:T          |          T0:
                            |
                            V
                          0,0,-1

       iid=-1:     7_______________________6
                   /                      /|
                  / |                    / |
                 /______________________/  |
                5|  |                   |4 |
                 |                      |  | 
                 |  |3 _ _ _ _ _ _ _ _ _| _|2
                 |  /                   |  /
                 | /                    | /
                 |______________________|/
                1                        0



              X(3)
              ^
              |    X(2)
              |   /
              |  /
              | /
              |/
              o---------> X(1)


 Connectivity:

Local hexagon has coordinates: 0,0,0
Neighbors: i,j,k; where {i,j,k}=-1,+1

Iv[0][*] -> Zero tetrahedron is in the center of the local hexagon
     => it only has neighbors in the same hexagon

**************************************************/

#define Jc(i0,i1,i2,m)	Mv*(i0+n0*(i1+n1*(i2)))+m
#define N(i,j,k)	(i)+(j)*ni+(k)*ni*nj

//?void	Box::connect(struct Cell *C)
//?{
//?	int
//?		n,
//?		i0,i1,i2,
//?		n0=nx[0],
//?		n1=nx[1],
//?		n2=nx[2],
//?		ni=n0+1,
//?		nj=n1+1,
//?		nk=n2+1;
//?#ifdef PARALLEL
//?	uint	*b;
//?#endif
//?	if(option.debug){printf("Connecting box domain ...");FLUSH;}
//?	/* Cell-node connectivity */
//?	for (i0=0; i0<n0; i0++)
//?	for (i1=0; i1<n1; i1++)
//?	for (i2=0; i2<n2; i2++)
//?	{
//?		int
//?			j,j0,j1,j2,jd,ja,jb,
//?//	   	in=N(i0,i1,i2), //primary node index
//?//			iv=in*DIM,      //vector index
//?			ic0=Jc(i0,i1,i2,0),//hexahedral cell index
//?			ic1=ic0+1,//hexahedral cell index
//?			ic2=ic0+2,//hexahedral cell index
//?			ic3=ic0+3,//hexahedral cell index
//?			ic4=ic0+4,//hexahedral cell index
//?			m[8];
//?		/* Checkerboard fashion to get connectivity right */
//?		j=(i0+i1+i2)%2;
//?		j0=j;j2=2-3*j;jd=1-2*j;
//?		j*=n0+1;
//?		ja=j -1;
//?		jb=n0-j;
//?		j1=-j0+1;
//?		n=0;
//?		for (int k=0; k<2; k++)
//?		for (int j=0; j<2; j++)
//?		for (int i=j0; i!=j2; i+=jd,n++)
//?			m[n]=N(i0+i,i1+j,i2+k);
//?		//CELL-VERTEX CONNECTIVITY
//?		//	Pick up Mv tetrahedral subelements 
//?		for (int iv=0; iv<Mv; iv++)
//?			for (int j=0; j<Nv; j++)
//?				C[ic0+iv].vert[j]=m[Iv[iv][j]];
//?  		// 0-tetra (central: {1,2,4,7}) has only neighbors within the same hexa
//?		//CELL-CELL CONNECTIVITY
//?		C[ic0].cell[3]=ic0+1; //	 3->124->1
//?		C[ic0].cell[0]=ic0+4; //	 0->247->4
//?		C[ic0].cell[1]=ic0+3; //	 1->471->3
//?		C[ic0].cell[2]=ic0+2; //	 2->712->2
//?		// 1-tetra {0,1,2,4}:                    this.face->neighbor.face
//?		C[ic1].cell[3]=(j=i2- 1)!=-1?Jc(i0,i1, j,3):bc[2*2+0];//012->457
//?		C[ic1].cell[0]=ic0;                                   //124->124
//?		C[ic1].cell[1]=(j=i0-jd)!=ja?Jc( j,i1,i2,1):bc[2*0+j0];//j0,240->240
//?		C[ic1].cell[2]=(j=i1- 1)!=-1?Jc(i0, j,i2,2):bc[2*1+0];//401->237
//?		// 2-tetra {1,2,3,7} 
//?		C[ic2].cell[3]=(j=i2- 1)!=-1?Jc(i0,i1, j,4):bc[2*2+0];//123->467	
//?		C[ic2].cell[0]=(j=i1+ 1)!=n1?Jc(i0, j,i2,1):bc[2*1+1];//237->104
//?		C[ic2].cell[1]=(j=i0+jd)!=jb?Jc( j,i1,i2,2):bc[2*0+j1];//j0,371->371
//?		C[ic2].cell[2]=ic0;                                    //712->712
//?		// 3-tetra {1,4,5,7}
//?		C[ic3].cell[3]=(j=i1- 1)!=-1?Jc(i0, j,i2,4):bc[2*1+0];//145->276
//?		C[ic3].cell[0]=(j=i2+ 1)!=n2?Jc(i0,i1, j,1):bc[2*2+1];//457->102
//?		C[ic3].cell[1]=(j=i0+jd)!=jb?Jc( j,i1,i2,3):bc[2*0+j1];//j0,571->571
//?		C[ic3].cell[2]=ic0;                                    //714->714
//?		// 4-tetra {2,4,6,7} -> 4230
//?		C[ic4].cell[3]=(j=i0-jd)!=ja?Jc( j,i1,i2,4):bc[2*0+j0];//j0,246->246
//?		C[ic4].cell[0]=(j=i2+ 1)!=n2?Jc(i0,i1, j,2):bc[2*2+1];//467->132
//?		C[ic4].cell[1]=(j=i1+ 1)!=n1?Jc(i0, j,i2,3):bc[2*1+1];//672->541
//?		C[ic4].cell[2]=ic0;                                   //724->724
//?	}
//?	if(option.debug){printf(" done\n");FLUSH;}
//?}
//?void	Box::create
//?(
//?	double xmin[], 
//?	double xmax[], 
//?	double *X,//coordinates array
//?	struct Node *node
//?)
//?{	/*	Create a 1x1x1 cube */
//?	int	n0=nx[0]+1,n1=nx[1]+1,n2=nx[2]+1,
//?		ni=n0,nj=n1,nk=n2;
//?	double
//?		x0=xmin[0],
//?		y0=xmin[1],
//?		z0=xmin[2],
//?		x1=xmax[0],
//?		y1=xmax[1],
//?		z1=xmax[2];
//?	if (option.debug)
//?	{	printf
//?		("Creating a %d x %d x %d box at (%g,%g) (%g,%g) (%g,%g)\n",
//?			nx[0],nx[1],nx[2],
//?			x0,x1,y0,y1,z0,z1
//?		);FLUSH;
//?	}
//?	for (int i2 = 0; i2 < n2; i2++)
//?	for (int i1 = 0; i1 < n1; i1++)
//?	for (int i0 = 0; i0 < n0; i0++)
//?	{	int	i=N(i0,i1,i2);
//?		double	*x=X+i*DIM;
//?		x[0]=x0+i0/(double)(n0-1)*(x1-x0);
//?		x[1]=y0+i1/(double)(n1-1)*(y1-y0);
//?		x[2]=z0+i2/(double)(n2-1)*(z1-z0);
//?		node[i].type=0;
//?	}	
//?	if (option.debug)
//?	{	//DEBUG OUTPUT
//?		double	totvol=0;
//?		printf("Total volume = ");FLUSH;
//?		for (int k=0; k<nk-1; k++)
//?		for (int j=0; j<nj-1; j++)
//?		for (int i=0; i<ni-1; i++)
//?		{	int	n[4];
//?			double	vol=1.0;
//?			n[0]=N(i,j,k);
//?			n[1]=N(i+1,j  ,k  );
//?			n[2]=N(i  ,j+1,k  );
//?			n[3]=N(i  ,j  ,k+1);
//?			for (int l=0; l<3; l++)
//?			{	double a=0.;
//?				for (int m=0; m<3; m++)
//?				{
//?				//	double c=node[n[l+1]].x[m]-node[n[0]].x[m];
//?					double c=X[DIM*n[l+1]+m]-X[DIM*n[0]+m];
//?					a+=c*c;
//?				}
//?				vol*=sqrt(a);
//?			}
//?			if (vol<SMALL)
//?			{
//?				printf("%d,%d,%d: vol=%g\n",i,j,k,vol);
//?				for (int l=0; l<4; l++)
//?				{	for (int m=0; m<3; m++)
//?						printf("%g\t",X[DIM*n[l]+m]);
//?					printf("\n");
//?				}
//?			}
//?			totvol+=vol;
//?		}
//?		printf("%g\n",totvol);FLUSH;
//?	}
//?}
void	Box::connect
(
	DNode	*node_root,
	DCell *cell_root
)
{	int
		n,
		i0,i1,i2,
		n0=nx[0],
		n1=nx[1],
		n2=nx[2],
		ni=n0+1,
		nj=n1+1,
		nk=n2+1,
		nc=5*n0*n1*n2,
		nn=ni*nj*nk;
	DNode	*node,**I;//Node index
	DCell	*cell,**C;//Cell index
#ifdef PARALLEL
	uint	*b;
#endif
	if(option.debug){printf("Connecting box domain ...");FLUSH;}
	if(Nv!=4)ERROR("Number of vertexes should be equal to 4\n");
	if(Nf!=4)ERROR("Number of faces should be equal to 4\n");
	/* Index all cells and nodes */
	I=new DNode*[nn];
	node=node_root;
	for(int i=0;i<nn;i++)
	{	I[i]=node;
		node=node->next;
	}
	C=new DCell*[nc];
	cell=cell_root;
	for(int i=0;i<nc;i++)
	{	C[i]=cell;
		cell->state.boundary=0;
		for(int j=0;j<Nf;j++)
			cell->facetype[j]=boundary;
		cell=cell->next;
	}
	/* Cell-node connectivity */
	for (i0=0; i0<n0; i0++)
	for (i1=0; i1<n1; i1++)
	for (i2=0; i2<n2; i2++)
	{
		int	
			j,j0,j1,j2,jd,ja,jb,
			ic0=Jc(i0,i1,i2,0),//hexahedral cell index
			ic1=ic0+1,//hexahedral cell index
			ic2=ic0+2,//hexahedral cell index
			ic3=ic0+3,//hexahedral cell index
			ic4=ic0+4,//hexahedral cell index
			m[8],ic[Nf];
		/* Checkerboard fashion to get connectivity right */
		j=(i0+i1+i2)%2;
		j0=j;j2=2-3*j;jd=1-2*j;
		j*=n0+1;
		ja=j -1;
		jb=n0-j;
		j1=-j0+1;
		n=0;
		for (int k=0; k<2; k++)
		for (int j=0; j<2; j++)
		for (int i=j0; i!=j2; i+=jd,n++)
			m[n]=N(i0+i,i1+j,i2+k);
		//CELL-VERTEX CONNECTIVITY
		//	Pick up Mv tetrahedral subelements 
		for (int iv=0; iv<Mv; iv++)
			for (int j=0; j<Nv; j++)
				C[ic0+iv]->vert[j]=I[m[Iv[iv][j]]];
			// 0-tetra (central: {1,2,4,7}) has only neighbors within the same hexa
		//CELL-CELL CONNECTIVITY
		C[ic0]->neib[3]=C[ic0+1]; //	 3->124->1
		C[ic0]->neib[0]=C[ic0+4]; //	 0->247->4
		C[ic0]->neib[1]=C[ic0+3]; //	 1->471->3
		C[ic0]->neib[2]=C[ic0+2]; //	 2->712->2
		for(int i=0;i<Nf;i++)
		{
			C[ic0]->facetype[i]=internal;
		}
		// 1-tetra {0,1,2,4}:                    this.face->neighbor.face
		ic[3]=(j=i2- 1)!=-1?Jc(i0,i1, j,3):bc[2*2+0];//012->457
		ic[0]=ic0;                                   //124->124
		ic[1]=(j=i0-jd)!=ja?Jc( j,i1,i2,1):bc[2*0+j0];//j0,240->240
		ic[2]=(j=i1- 1)!=-1?Jc(i0, j,i2,2):bc[2*1+0];//401->237
		for(int i=0;i<Nf;i++)
		if(ic[i]>=0)
		{	C[ic1]->neib[i]=C[ic[i]];
			C[ic1]->facetype[i]=internal;
		}
		else
		{	C[ic1]->neib[i]=NULL;
			C[ic1]->state.boundary=1;
			C[ic1]->facetype[i]=(ElementStatus)ic[i];//012->457
		}
		// 2-tetra {1,2,3,7} 
		ic[3]=(j=i2- 1)!=-1?Jc(i0,i1, j,4):bc[2*2+0];//123->467	
		ic[0]=(j=i1+ 1)!=n1?Jc(i0, j,i2,1):bc[2*1+1];//237->104
		ic[1]=(j=i0+jd)!=jb?Jc( j,i1,i2,2):bc[2*0+j1];//j0,371->371
		ic[2]=ic0;                                    //712->712
		for(int i=0;i<Nf;i++)
		if(ic[i]>=0)
		{	C[ic2]->neib[i]=C[ic[i]];
			C[ic2]->facetype[i]=internal;
		}
		else
		{	C[ic2]->neib[i]=NULL;
			C[ic2]->state.boundary=1;
			C[ic2]->facetype[i]=(ElementStatus)ic[i];//012->457
		}
		// 3-tetra {1,4,5,7}
		ic[3]=(j=i1- 1)!=-1?Jc(i0, j,i2,4):bc[2*1+0];//145->276
		ic[0]=(j=i2+ 1)!=n2?Jc(i0,i1, j,1):bc[2*2+1];//457->102
		ic[1]=(j=i0+jd)!=jb?Jc( j,i1,i2,3):bc[2*0+j1];//j0,571->571
		ic[2]=ic0;                                    //714->714
		for(int i=0;i<Nf;i++)
		if(ic[i]>=0)
		{	C[ic3]->neib[i]=C[ic[i]];
			C[ic3]->facetype[i]=internal;
		}
		else
		{	C[ic3]->neib[i]=NULL;
			C[ic3]->state.boundary=1;
			C[ic3]->facetype[i]=(ElementStatus)ic[i];//012->457
		}
		// 4-tetra {2,4,6,7} -> 4230
		ic[3]=(j=i0-jd)!=ja?Jc( j,i1,i2,4):bc[2*0+j0];//j0,246->246
		ic[0]=(j=i2+ 1)!=n2?Jc(i0,i1, j,2):bc[2*2+1];//467->132
		ic[1]=(j=i1+ 1)!=n1?Jc(i0, j,i2,3):bc[2*1+1];//672->541
		ic[2]=ic0;                                   //724->724
		for(int i=0;i<Nf;i++)
		if(ic[i]>=0)
		{	C[ic4]->neib[i]=C[ic[i]];
			C[ic4]->facetype[i]=internal;
		}
		else
		{	C[ic4]->neib[i]=NULL;
			C[ic4]->state.boundary=1;
			C[ic4]->facetype[i]=(ElementStatus)ic[i];//012->457
		}
	}
	delete I;
	delete C;
	if(option.debug){printf(" done\n");FLUSH;}
}
void	Box::create
(
	double xmin[], 
	double xmax[], 
	struct DNode *root
)
{	/*	Create a 1x1x1 cube */
	int	n0=nx[0]+1,n1=nx[1]+1,n2=nx[2]+1,
		ni=n0,nj=n1,nk=n2,nn=ni*nj*nk;
	double
		x0=xmin[0],
		y0=xmin[1],
		z0=xmin[2],
		x1=xmax[0],
		y1=xmax[1],
		z1=xmax[2];
	DNode	*node=root,**I;//Node index
	if (option.debug)
	{	printf
		("Creating a %d x %d x %d box at (%g,%g) (%g,%g) (%g,%g)\n",
			nx[0],nx[1],nx[2],
			x0,x1,y0,y1,z0,z1
		);FLUSH;
	}
	/* Index all nodes */
	I=new DNode*[nn];
	node=root;
	for(int i=0;i<nn;i++)
	{	I[i]=node;
		node=node->next;
	}
	for (int i2 = 0; i2 < n2; i2++)
	for (int i1 = 0; i1 < n1; i1++)
	for (int i0 = 0; i0 < n0; i0++)
	{	int	i=N(i0,i1,i2);
		double	*x=I[i]->x;
		x[0]=x0+i0/(double)(n0-1)*(x1-x0);
		x[1]=y0+i1/(double)(n1-1)*(y1-y0);
		x[2]=z0+i2/(double)(n2-1)*(z1-z0);
	}
	if (option.debug)
	{	//DEBUG OUTPUT
		double	totvol=0;
		printf("Total volume = ");FLUSH;
		for (int k=0; k<nk-1; k++)
		for (int j=0; j<nj-1; j++)
		for (int i=0; i<ni-1; i++)
		{	int	n[4];
			double	vol=1.0;
			n[0]=N(i,j,k);
			n[1]=N(i+1,j  ,k  );
			n[2]=N(i  ,j+1,k  );
			n[3]=N(i  ,j  ,k+1);
			for (int l=0; l<3; l++)
			{	double a=0.;
				for (int m=0; m<3; m++)
				{
				//	double c=node[n[l+1]].x[m]-node[n[0]].x[m];
					double c=I[n[l+1]]->x[m]-I[n[0]]->x[m];
					a+=c*c;
				}
				vol*=sqrt(a);
			}
			if (vol<SMALL)
			{	printf("%d,%d,%d: vol=%g\n",i,j,k,vol);
				for (int l=0; l<4; l++)
				{	for (int m=0; m<3; m++)
						printf("%g\t",I[n[l]]->x[m]);
					printf("\n");
				}
			}
			totvol+=vol;
		}
		delete I;
		printf("%g\n",totvol);FLUSH;
	}
}	

