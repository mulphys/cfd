/* 
 * Generates cylindrical grid
 */
#include	<stdio.h>
#include	<math.h>
#include	"main.h"
#include	"io.h"
#include	"geom.h"
#include	"vecalg.h"
#include	"var.h"
#include	"tool.h"
#include	"domain.h"
#include	"cylinder.h"

//?void	Cylinder::ring
//?(
//?	int	i,//index of the ring
//?	int	ip,//index of the cross-sectional plane
//?	double	*GX,//global coordinates
//?	Cell	*C //connectivity
//?)
//?{	int	nnd,inpr,inpr0,inpr1,nnr,
//?	inp=nnp*ip;//start index of the nodes on the plane
//?	double	pi=4.*atan(1.),
//?			pio3=pi/3.,//sector
//?			dr,r,dz,z,
//?			*X;
//?	dr=R/(double)nr;
//?	r=i*dr;
//?	dz=L/(double)ns;
//?	z=ip*dz-L/2.;
//?	//Set coordinates
//?	nnd=NNR(i-1);//nodes on the disk of radius i-1
//?	inpr=inp+nnd;//start node-index of the ring
//?	nnr=i>0?6*i:1;//number of nodes in the ring
//?	inpr1=inpr+nnr;//end node-index of the ring
//?	X=GX+inpr*DIM;
//?	if(i==0)
//?	{
//?		X[0]=X[1]=0.0;X[2]=z;
//?	}
//?	else
//?	{	int	j,isec;
//?		double	*x,dt,t;
//?		t=0.0;
//?		dt=pio3/(double)i;
//?		for (isec=0; isec<6; isec++)
//?		for (j=0; j<i; j++,t+=dt)
//?		{
//?			x=X+DIM*(isec*i+j);
//?			x[0]=r*sin(t);x[1]=r*cos(t);x[2]=z;
//?		}
//?		if (ip>0)
//?		{//Create cells
//?			int
//?				isec,it,j,
//?				ittop,itbot,itout,itclo,itcnt,
//?				in0,in1,in2,in3,in4,in5,
//?				ic,ictop,icbot,icsid,itint,
//?				itthis,
//?				its,itsr,itsr0,itsr1,
//?				ics,icsr,
//?				itop=0,ibot=1,isid=2,
//?				nnr0,
//?				nnd0=NNR(i-2),//nodes on the disk of radius i-2
//?				inpr0=inp+nnd0,//start node-index of the ring i-1
//?				is=ip-1;//index to the current disk-slice
//?			nnr0=i>1?6*(i-1):1;//number of nodes in the inner-ring
//?			ntr=(2*i-1)*6; //number of trapezoids in the ring
//?			its=nts*is;//start trapeziod index of the slice
//?			itsr=its+NTD(i-1);//start-trapeziod index of this ring
//?			itsr0=i>1?its+NTD(i-2):0;//start-trapeziod index of the inner ring
//?			itsr1=its+NTD(i);//start-trapeziod index of the outter ring
//?			ics=ncs*is;//start-cell index of the slice
//?			icsr=ics+NCR(i-1);//start-cell index of the ring
//?			t=0.0;
//?			for (isec=0; isec<6; isec++)
//?			{	for (j=0; j<i; j++,t+=dt)
//?				{	//outer subset
//?					//tetra-types: top,bot,side
//?					//Preffixes:
//?					//it=index to the trapezoid
//?					//ic=index to the cell
//?					//in=index to the node
//?					//Suffixes:
//?					//clo=clockwise
//?					//cnt=counterclocwise
//?					it=isec*(2*i-1)+2*j;//trapez index inside the ring
//?					itthis=itsr+it;//inside the cylinder
//?					//Neighbor trapeziods:
//?					ittop=itthis+nts;
//?					itbot=itthis-nts;
//?					itout=i<nr?itsr1+isec*(2*i+1)+2*j+1:-1;
//?					itclo=itsr+(it+ntr-1)%ntr;
//?					itcnt=itsr+(it+ntr+1)%ntr;
//?
//?					in0=inpr0+(isec*(i-1)+j)%nnr0;
//?					in1=inpr+isec*i+j;
//?					in2=inpr+(isec*i+j+1)%nnr;
//?					in3=in0-nnp;
//?					in4=in1-nnp;
//?					in5=in2-nnp;
//?
//?					ic=icsr+3*it;
//?					ictop=ic+itop;
//?					icbot=ic+ibot;
//?					icsid=ic+isid;
//?
//?					C[ictop].vert[0]=in0;
//?					C[ictop].vert[1]=in1;
//?					C[ictop].vert[2]=in2;
//?					C[ictop].vert[3]=in3;
//?
//?					C[ictop].cell[0]=icsid;
//?					C[ictop].cell[1]=NB(itcnt,itop);
//?					C[ictop].cell[2]=NB(itclo,j>0?isid:itop);
//?					C[ictop].cell[3]=NB(ittop,ibot);
//?
//?					C[icbot].vert[0]=in2;
//?					C[icbot].vert[1]=in3;
//?					C[icbot].vert[2]=in4;
//?					C[icbot].vert[3]=in5;
//?	
//?					C[icbot].cell[0]=NB(itbot,itop);
//?					C[icbot].cell[1]=NB(itout,isid);
//?					C[icbot].cell[2]=NB(itcnt,j<i-1?ibot:isid);
//?					C[icbot].cell[3]=icsid;
//?	
//?					C[icsid].vert[0]=in2;
//?					C[icsid].vert[1]=in3;
//?					C[icsid].vert[2]=in4;
//?					C[icsid].vert[3]=in1;
//?	
//?					C[icsid].cell[0]=NB(itclo,ibot);
//?					C[icsid].cell[1]=NB(itout,itop);
//?					C[icsid].cell[2]=ictop;
//?					C[icsid].cell[3]=icbot;
//?				}
//?				for (j=0; j<i-1; j++,t+=dt)
//?				{//inner subset
//?					it=isec*(2*i-1)+2*j+1;//trapez index inside the ring
//?					itthis=itsr+it;//inside the cylinder
//?					//Neighbor trapeziods:
//?					ittop=itthis+nts;
//?					itbot=itthis-nts;
//?					itint=itsr0+isec*(2*(i-1)-1)+2*j;
//?					itclo=itsr+(it+ntr-1)%ntr;
//?					itcnt=itsr+(it+ntr+1)%ntr;
//?	
//?					in0=inpr+(isec*i+j+1)%nnr;
//?					in1=inpr0+(isec*(i-1)+j);
//?					in2=inpr0+(isec*(i-1)+j+1)%nnr0;
//?					in3=in0-nnp;
//?					in4=in1-nnp;
//?					in5=in2-nnp;
//?	
//?					ic=icsr+3*it;
//?					ictop=ic+itop;
//?					icbot=ic+ibot;
//?					icsid=ic+isid;
//?	
//?					C[ictop].vert[0]=in0;
//?					C[ictop].vert[1]=in1;
//?					C[ictop].vert[2]=in2;
//?					C[ictop].vert[3]=in4;
//?	
//?					C[ictop].cell[0]=NB(itint,isid);
//?					C[ictop].cell[1]=icsid;
//?					C[ictop].cell[2]=NB(itclo,itop);
//?					C[ictop].cell[3]=NB(ittop,ibot);
//?	
//?					C[icbot].vert[0]=in0;
//?					C[icbot].vert[1]=in3;
//?					C[icbot].vert[2]=in4;
//?					C[icbot].vert[3]=in5;
//?	
//?					C[icbot].cell[0]=NB(itbot,itop);
//?					C[icbot].cell[1]=icsid;
//?					C[icbot].cell[2]=NB(itcnt,isid);
//?					C[icbot].cell[3]=NB(itclo,ibot);
//?	
//?					C[icsid].vert[0]=in0;
//?					C[icsid].vert[1]=in2;
//?					C[icsid].vert[2]=in4;
//?					C[icsid].vert[3]=in5;
//?	
//?					C[icsid].cell[0]=NB(itint,ibot);
//?					C[icsid].cell[1]=icbot;
//?					C[icsid].cell[2]=NB(itcnt,itop);
//?					C[icsid].cell[3]=ictop;
//?				}
//?			}
//?		}
//?	}
//?}
void	Cylinder::ring
(	//CREATES CYLINDRICAL MESH 
	int	ir,//index of the ring
	int	ip,//index of the cross-sectional plane
	DNode	**I,
	DCell	**C
)
{	int	nnd,inpr,inpr0,inpr1,nnr,
	inp=nnp*ip;//start index of the nodes on the plane
	double	pi=4.*atan(1.),
			pio3=pi/3.,//sector
			dr,r,dz,z;
	DNode	**X;
	dr=R/(double)nr;
	r=ir*dr;
	dz=L/(double)ns;
	z=ip*dz-L/2.;
	//Set coordinates
	nnd=NNR(ir-1);//nodes on the disk of radius i-1
	inpr=inp+nnd;//start node-index of the ring
	nnr=ir>0?6*ir:1;//number of nodes in the ring
	inpr1=inpr+nnr;//end node-index of the ring
	X=I+inpr;
	if(ir==0)
	{
		X[0]->x[0]=X[0]->x[1]=0.0;X[0]->x[2]=z;
	}
	else
	{	int	j,isec;
		double	*x,dt,t;
		t=0.0;
		dt=pio3/(double)ir;
		for (isec=0; isec<6; isec++)
		for (j=0; j<ir; j++,t+=dt)
		{
			x=X[isec*ir+j]->x;
			x[0]=r*sin(t);x[1]=r*cos(t);x[2]=z;
		}
		if (ip>0)
		{//Create cells
			int
				isec,it,j,
				ittop,itbot,itout,itclo,itcnt,
				in0,in1,in2,in3,in4,in5,
				ic,ictop,icbot,icsid,itint,
				itthis,
				its,itsr,itsr0,itsr1,
				ics,icsr,
				itop=0,ibot=1,isid=2,
				nnr0,
				nnd0=NNR(ir-2),//nodes on the disk of radius i-2
				inpr0=inp+nnd0,//start node-index of the ring i-1
				is=ip-1;//index to the current disk-slice
			nnr0=ir>1?6*(ir-1):1;//number of nodes in the inner-ring
			ntr=(2*ir-1)*6; //number of trapezoids in the ring
			its=nts*is;//start trapeziod index of the slice
			itsr=its+NTD(ir-1);//start-trapeziod index of this ring
			itsr0=ir>1?its+NTD(ir-2):0;//start-trapeziod index of the inner ring
			itsr1=its+NTD(ir);//start-trapeziod index of the outter ring
			ics=ncs*is;//start-cell index of the slice
			icsr=ics+NCR(ir-1);//start-cell index of the ring
			t=0.0;
			for (isec=0; isec<6; isec++)
			{	for (j=0; j<ir; j++,t+=dt)
				{	//outer subset
					//tetra-types: top,bot,side
					//Preffixes:
					//it=index to the trapezoid
					//ic=index to the cell
					//in=index to the node
					//Suffixes:
					//clo=clockwise
					//cnt=counterclocwise
					int	vert[Nv],neib[Nf];
					it=isec*(2*ir-1)+2*j;//trapez index inside the ring
					itthis=itsr+it;//inside the cylinder
					//Neighbor trapeziods:
					ittop=itthis+nts;
					itbot=itthis-nts;
					itout=ir<nr?itsr1+isec*(2*ir+1)+2*j+1:-1;
					itclo=itsr+(it+ntr-1)%ntr;
					itcnt=itsr+(it+ntr+1)%ntr;

					in0=inpr0+(isec*(ir-1)+j)%nnr0;
					in1=inpr+isec*ir+j;
					in2=inpr+(isec*ir+j+1)%nnr;
					in3=in0-nnp;
					in4=in1-nnp;
					in5=in2-nnp;

					ic=icsr+3*it;
					ictop=ic+itop;
					icbot=ic+ibot;
					icsid=ic+isid;

					vert[0]=in0;
					vert[1]=in1;
					vert[2]=in2;
					vert[3]=in3;
					for(int jv=0;jv<Nv;jv++)
						C[ictop]->vert[jv]=I[vert[jv]];

					neib[0]=icsid;
					neib[1]=NB(itcnt,itop);
					neib[2]=NB(itclo,j>0?isid:itop);
					neib[3]=NB(ittop,ibot);
					for(int jf=0;jf<Nf;jf++)
					{	DCell	*c=C[ictop];
						if(neib[jf]>=0&&neib[jf]<nc)
						{	c->neib[jf]=C[neib[jf]];
							c->facetype[jf]=internal;
						}
						else
						{
							c->neib[jf]=NULL;
							c->facetype[jf]=(ElementStatus)neib[jf];
							c->state.boundary=1;
						}
					}
					vert[0]=in2;
					vert[1]=in3;
					vert[2]=in4;
					vert[3]=in5;
					for(int jv=0;jv<Nv;jv++)
						C[icbot]->vert[jv]=I[vert[jv]];

					neib[0]=NB(itbot,itop);
					neib[1]=NB(itout,isid);
					neib[2]=NB(itcnt,j<ir-1?ibot:isid);
					neib[3]=icsid;
					for(int jf=0;jf<Nf;jf++)
					{	DCell	*c=C[icbot];
						if(neib[jf]>=0&&neib[jf]<nc)
						{	c->neib[jf]=C[neib[jf]];
							c->facetype[jf]=internal;
						}
						else
						{
							c->neib[jf]=NULL;
							c->facetype[jf]=(ElementStatus)neib[jf];
							c->state.boundary=1;
						}
					}
					vert[0]=in2;
					vert[1]=in3;
					vert[2]=in4;
					vert[3]=in1;
					for(int jv=0;jv<Nv;jv++)
						C[icsid]->vert[jv]=I[vert[jv]];

					neib[0]=NB(itclo,ibot);
					neib[1]=NB(itout,itop);
					neib[2]=ictop;
					neib[3]=icbot;
					for(int jf=0;jf<Nf;jf++)
					{	DCell	*c=C[icsid];
						if(neib[jf]>=0&&neib[jf]<nc)
						{	c->neib[jf]=C[neib[jf]];
							c->facetype[jf]=internal;
						}
						else
						{
							c->neib[jf]=NULL;
							c->facetype[jf]=(ElementStatus)neib[jf];
							c->state.boundary=1;
						}
					}
				}
				for (j=0; j<ir-1; j++,t+=dt)
				{//inner subset
					int	vert[Nv],neib[Nf];
					it=isec*(2*ir-1)+2*j+1;//trapez index inside the ring
					itthis=itsr+it;//inside the cylinder
					//Neighbor trapeziods:
					ittop=itthis+nts;
					itbot=itthis-nts;
					itint=itsr0+isec*(2*(ir-1)-1)+2*j;
					itclo=itsr+(it+ntr-1)%ntr;
					itcnt=itsr+(it+ntr+1)%ntr;
	
					in0=inpr+(isec*ir+j+1)%nnr;
					in1=inpr0+(isec*(ir-1)+j);
					in2=inpr0+(isec*(ir-1)+j+1)%nnr0;
					in3=in0-nnp;
					in4=in1-nnp;
					in5=in2-nnp;
	
					ic=icsr+3*it;
					ictop=ic+itop;
					icbot=ic+ibot;
					icsid=ic+isid;
	
					vert[0]=in0;
					vert[1]=in1;
					vert[2]=in2;
					vert[3]=in4;
					for(int jv=0;jv<Nv;jv++)
						C[ictop]->vert[jv]=I[vert[jv]];
	
					neib[0]=NB(itint,isid);
					neib[1]=icsid;
					neib[2]=NB(itclo,itop);
					neib[3]=NB(ittop,ibot);
					for(int jf=0;jf<Nf;jf++)
					{	DCell	*c=C[ictop];
						if(neib[jf]>=0&&neib[jf]<nc)
						{	c->neib[jf]=C[neib[jf]];
							c->facetype[jf]=internal;
						}
						else
						{
							c->neib[jf]=NULL;
							c->facetype[jf]=(ElementStatus)neib[jf];
							c->state.boundary=1;
						}
					}
					vert[0]=in0;
					vert[1]=in3;
					vert[2]=in4;
					vert[3]=in5;
					for(int jv=0;jv<Nv;jv++)
						C[icbot]->vert[jv]=I[vert[jv]];
					neib[0]=NB(itbot,itop);
					neib[1]=icsid;
					neib[2]=NB(itcnt,isid);
					neib[3]=NB(itclo,ibot);
					for(int jf=0;jf<Nf;jf++)
					{	DCell	*c=C[icbot];
						if(neib[jf]>=0&&neib[jf]<nc)
						{	c->neib[jf]=C[neib[jf]];
							c->facetype[jf]=internal;
						}
						else
						{
							c->neib[jf]=NULL;
							c->facetype[jf]=(ElementStatus)neib[jf];
							c->state.boundary=1;
						}
					}
					vert[0]=in0;
					vert[1]=in2;
					vert[2]=in4;
					vert[3]=in5;
					for(int jv=0;jv<Nv;jv++)
						C[icsid]->vert[jv]=I[vert[jv]];

					neib[0]=NB(itint,ibot);
					neib[1]=icbot;
					neib[2]=NB(itcnt,itop);
					neib[3]=ictop;
					for(int jf=0;jf<Nf;jf++)
					{	DCell	*c=C[icsid];
						if(neib[jf]>=0&&neib[jf]<nc)
						{	c->neib[jf]=C[neib[jf]];
							c->facetype[jf]=internal;
						}
						else
						{
							c->neib[jf]=NULL;
							c->facetype[jf]=(ElementStatus)neib[jf];
							c->state.boundary=1;
						}
					}
				}
			}
		}
	}
}
void	Cylinder::init
(	int	*nx,  //nx[0]: number of rings; nx[1]: number of slices
	double	*dd //dd[0]: cylinder radius; dd[1]: cylinder length
)
{
	R=dd[0];
	L=dd[1];
	nr=nx[0];//number of rings in the cross-sectional plane
	ni=nr+1;
	ns=nx[1];//number of disks in the length of the cylinder
	np=ns+1;//number of cross-sectional planes
	nts=NTD(nr);//trapezoids in a cross-section slice
	ncs=NCR(nr);//cells in a cross-section slice
	nnp=NNR(nr);//nodes on the cross-sectional plane
	nt=nts*ns;//trapezoids in the cylinder
	nn=np*nnp;//nodes in the cylinder
	nc=ncs*ns;//cells in the cylinder
}
//?void	Cylinder::create
//?(
//?	double	*X,
//?	struct Node	*node_root,
//?	struct Cell	*cell_root
//?)
//?{
//?	if(option.verbose)
//?	{	printf("Creating cylinder block ... ");FLUSH;
//?	}
//?	for (int ip=0; ip<np; ip++)
//?	for (int i=0; i<ni; i++)
//?		ring(i,ip,X,cell_root); 
//?	if(option.verbose)
//?	{	printf("done\n");FLUSH;
//?	}
//?	//Set boundary conditions
//?	if(node_root!=NULL)
//?	{	double
//?			dr=R/(double)nr,
//?			dz=L/(double)ns,
//?			zmin=0.5*(-L+dz),
//?			zmax=0.5*( L-dz),
//?			rmax=R-0.5*dr;
//?		for(int inode=0;inode<nn;inode++)
//?		{	int	jnode=DIM*inode;
//?			double
//?				x=X[jnode+0],
//?				y=X[jnode+1],
//?				r=sqrt(x*x+y*y),
//?				z=X[jnode+2];
//?			Node	*node=node_root+inode;
//?			if(z<zmin)node->type=bc[2];
//?			else
//?			if(z>zmax)node->type=bc[0];
//?			else
//?			if(r>rmax)node->type=bc[1];
//?			else	node->type=internal;
//?		}
//?	}
//?	if(cell_root!=NULL)
//?	{
//?		for(int icell=0;icell<nc;icell++)
//?		{	Cell	*cell=cell_root+icell;
//?			int
//?				*verts=cell->vert,
//?				*neibs=cell->cell;
//?			for(int iface=0;iface<Nf;iface++)
//?			{	if(neibs[iface]<0)
//?				{	int	facetype=neibs[iface];
//?					for(int iv=0;iv<Nfv;iv++)
//?					{	int
//?							ivert=(iface+iv+1)%Nv,//valid only if Nf==Nv
//?							inode=verts[ivert],jnode=DIM*inode;
//?						Node	*node=node_root+inode;
//?						if(node->type!=internal)neibs[iface]=node->type;
//?					}
//?				}
//?			}
//?		}
//?	}
//?}
void	Cylinder::create
(
	DNode	*node_root,
	DCell	*cell_root
)
{
	DNode	**I,*node;
	DCell	**C,*cell;
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
		{	//cell->facetype[j]=internal;
			cell->facetype[j]=boundary;
			cell->neib[j]=NULL;
		}
		cell=cell->next;
	}
	if(option.verbose)
	{	printf("Creating cylinder block ... ");FLUSH;
	}
	for (int ip=0; ip<np; ip++)
	for (int i=0; i<ni; i++)
		ring(i,ip,I,C); 
	if(option.verbose)
	{	printf("done\n");FLUSH;
	}
	delete I;
	delete C;
	//Set boundary conditions
	if(node_root!=NULL)
	{
		double
			dz=L/(double)ns,
			dr=R/(double)nr,
			rmax=R-0.5*dr,
			zmin=0.5*(-L+dz),
			zmax=0.5*( L-dz);
		DNode	*node=node_root;
		do
		{
			double
				x=node->x[0],
				y=node->x[1],
				z=node->x[2],
				r=sqrt(x*x+y*y);
			node->state.boundary=1;
			if(z<zmin)node->type=bc[2];
			else
			if(z>zmax)node->type=bc[0];
			else
			if(r>rmax)node->type=bc[1];
			else
			{	node->type=internal;
				node->state.boundary=0;
			}
			node=node->next;
		}	while(node!=node_root);
	}
}
Cylinder::~Cylinder()
{
}
int	Cylinder::getnodes()
{
	return nn;
}
int	Cylinder::getcells()
{
	return nc;
}
void	Cylinder::setbc(int b[])
{
	for (int i=0; i<3; i++)
		bc[i]=b[i];
	//TODO:
	//Construct the boundary face list
	//For each boundary face do:
	//IF the boundary face normal vector is up -> set bc[0]
	//IF ...                              down -> set bc[2]
	//ELSE                                     -> set bc[1]
	//...
}
void	Cylinder::setgrid(int n[])
{
	for (int i=0; i<2; i++)
		nx[i]=n[i];
}

