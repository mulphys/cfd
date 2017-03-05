#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
///#include <bzlib.h>
#include "vecalg.h"
#include "io.h"
#include "main.h"
#include "geom.h"
#include "var.h"
#include "tool.h"
#include "domain.h"
#include "output.h"
///#include "box.h"
///#include "cylinder.h"
#include "templates.cc"

#ifdef WITH_MPI
#include "mpi.h"
using namespace MPI;
#endif

#define MAXGRIDCELLS	1000000000

Domain::Domain()
{
	model=-1;
	fixedlimits=0;
	nvar=0;
	np=mp=0;
	nbv=0;//number of boundary vertexes
	nbf=0;//number of boundary faces
	nbc=0;//number of boundary cells
	dnode_root=NULL;
	bnode_root=NULL;
	dcell_root=NULL;
	bcell_root=NULL;
	bface_root=NULL;
#ifdef WITH_DOVE
//	dove_root=NULL;
	dove_ring = new Ring<Dove*>();
#endif
	tool=NULL;
///	outype=default_output;
	for (int i=0;i<maxelements;i++)
		Ne[i]=0;
}
Domain::~Domain()
{
	if (type==dynamic)
	{	//Deleting dynamic grid
		deleteGrid(dcell_root);
		deleteGrid(dnode_root);
	}
	delete variable;
///	delete segm;
}
int	Domain::getnvar()
{
	return nvar;
}
void	Domain::allocvar(int type, Variable *var)
{
	int
		dim=var[type].dimension,
		size=Ne[type];
	ALLOC(var[type].val,double,size*dim);
	var[type].size=size;
printf("Domain::allocvar:size=%d\n",size);fflush(stdout);//-
}
void	Domain::allocgeovar(Variable	*var)
{
	for (int type=0; type<maxelements; type++)
	{	int
			dim=var[type].dimension,
			size=Ne[type];
		ALLOC(var[type].val,double,size*dim);
		var[type].size=size;
printf("Domain::allocgeovar:size=%d\n",size);fflush(stdout);//-
	}
}
void	Domain::setDisplay()
{
	for (int i=0; i<maxelements; i++)
	{	disp[i].lighting=1;
	}
	dispvar=-1;
	dispvarcomp=0;
}
void	Domain::setDomain()
{	int
		varbufsize[maxelements],
		idomain=this-domain_root;
	if(option.debug){printf("Setting Domain %d:",idomain+1);FLUSH;}
	setGeom(configfile);
printf("Domain::setDomain:0:variable=%d\n",variable);//-
	setModel();
printf("Domain::setDomain:1:variable=%d\n",variable);//-
//?	if (type!=dynamic)
//?	{
//?		setFaces();
//?		setBoundaryVertexes();
//?		setVolumes();
//?	}
	setDisplay();
printf("Domain::setDomain:2:variable=%d\n",variable);//-

	//Set variables
	///nvar=maxvar[idomain];
	if(option.debug|option.verbose)
	{	printf("Domain %s: no. of variables: %d\n",name,nvar);
		printf
		(	"Allocating %d bytes for %d variables headers\n",
			sizeof(class Variable)*nvar,nvar
		);FLUSH;
	}
	variable=new Variable[nvar];
printf("Domain::setDomain:3:nvar=%d\n",nvar);fflush(stdout);//-
	//	domain->nproc=countChildren(i,DOMAIN_KEYWORD,PROCEDURE_KEYWORD,configfile);
	//	domain->procedure=new Procedure[domain->nproc];
	//Setup variables
	for (int i=0; i<maxelements; i++)varbufsize[i]=0;
	for (int ivar=0; ivar<nvar; ivar++)
	{	Variable	*var=variable+ivar;
		defvar(ivar,var->name,var->type,var->dimension,this);
		var->setinit(var->dimension);
		if (type==dynamic&&var->type!=points) 
		{	int	dim=var->dimension;
			var->size=0;
			var->loc=varbufsize[var->type];
			varbufsize[var->type]+=dim;
		}
		else
		{	var->size=Ne[var->type];
			if (option.verbose|option.debug)
			{	printf
				(	"Allocating %d*%d*%d=%d bytes for variable '%s', type=%s",
					var->size,var->dimension,sizeof(double),
					var->size*var->dimension*sizeof(double),
					var->name,elementype[var->type]
				);
				if (type!=points) 
					printf(", offset: %d",var->loc);
				printf("\n");
			}
			ALLOC(var->val,double,var->size*var->dimension);
printf("Domain::setDomain:3:ivar=%d, variable=%ld, size=%d * dim=%d = %d\n",ivar,variable,var->size, var->dimension, var->size*var->dimension);//-
		}
	}
//?	if (type==dynamic)
//?	{	//Allocating dynamic variables
		deleteRing(bface_root);
		createBoundaryFaceList(bface_root);
		for (int vartype=0; vartype<maxelements; vartype++)
		{
			if(option.verbose&&varbufsize[vartype]>0)
				printf
				(	"Allocating %d bytes for variables at each element of type '%s'\n",
					varbufsize[vartype],elementype[vartype]
				);
			switch(vartype)
			{
				case nodes:
					allocElementVars(varbufsize[vartype],dnode_root);
					break;
				case cells:
					allocElementVars(varbufsize[vartype],dcell_root);
	///				allocCellVars(varbufsize[vartype],varbufsize[boundary_faces],dcell_root);
					break;
	//			case boundary_nodes:
	//				break;
				case boundary_faces:
					allocElementVars(varbufsize[vartype],bface_root);
					break;
				default:
					if (varbufsize[vartype]>0)
					{	fprintf(stderr,"Can't have dynamic variable of type %d: %s\n",vartype,elementype[vartype]);
						exit(1);
					}
			}
		}
//?	}
printf("Domain::setDomain:4:variable=%ld\n",variable);//-
	//Initialize all variables (user-routines)
	for (int ivar=0; ivar<nvar; ivar++)
	{	Variable	*var=variable+ivar;
		int	size=type==dynamic&&var->type!=points?var->loc:var->size;
		initvar(ivar,var->name,var->type,size,var->dimension,var->val,this);
		if(option.verbose|option.debug)
		printf
		(	"Variable %d initialized: %s, type=%d, dim=%d, size=%d\n",
			ivar+1,var->name,var->type,var->dimension,var->size
		);fflush(stdout);
	}
printf("Domain::setDomain:END:variable=%d\n",variable);//-
}
int	Domain::Steps //Advance all variables by niter iterations 
(
	int niter
)
{
	for (int i=0; i<niter; i++)
	{
 		(*step)(runtime.step,this);
		/******************/

#ifdef WITH_MPI
		//double x = fabs(runtime.current-(int)runtime.current);
		double x = fmod(runtime.current, 1.0);
		/*too much communication overhead here: 
		if each step the variables are sent to proc 0 
		if use the following if statement, then overhead is less,
		but the animation may not so smooth: not real time */
		if(x<runtime.step)
		{
			int ivar=getDisplayVariable();
			if(iproc==0)
			{	//recv variables from other processes
				for(int j=0; j<ndomains; j++)
				{
					Domain* domj= domain_root+j;
					if(domj->iproc!=0)
					{	
						int ivar=domj->getDisplayVariable();
						int tag= j;
						DNode	*node_root=domj->dnode_root,
								*node=node_root;
	
						double *var=node->var;
						int	varloc=domj->variable[ivar].loc,
							lenvar=domj->variable[ivar].dimension;
						int n;
						COMM_WORLD.Recv(&n, 1, MPI_INT, domj->iproc, tag); 
						double* recvVarBuffer=new double[n*lenvar];
						COMM_WORLD.Recv(recvVarBuffer, n*lenvar, MPI_DOUBLE, domj->iproc, tag); 
						int k=0;
						do
						{	double
								*v=node->var + varloc;
							for(int j=0; j<lenvar; j++)
							{
								v[j]=recvVarBuffer[k*lenvar+j];
							}
							node=node->next;
							k++;
						}	while(node!=node_root);
						delete recvVarBuffer;
					}
				}
			}
			else //send variables to proc0
				sendDispVar2Proc0(ivar);
		}
#endif

	/***************************/

	}
	return	niter;
}


/************************ added by zhang *************/

#ifdef WITH_MPI
void	Domain::sendDispVar2Proc0(int ivar)
{
	DNode	*node_root=dnode_root,
			*node=node_root;

	int	n=n_nodes;
	if(ivar<0) return;
	if(node_root==NULL) return;
	int	varloc=variable[ivar].loc,
		lenvar=variable[ivar].dimension;
	double* sendVarBuffer=new double[n*lenvar];
	int i=0;
	do
	{	double
			*v=node->var + varloc;
		for(int j=0; j<lenvar; j++)
		{
			sendVarBuffer[i*lenvar+j]=v[j];
		}
		node=node->next;
		i++;
	}	while(node!=node_root);

	int tag = this-domain_root;
	COMM_WORLD.Send(&n, 1, MPI_INT, 0, tag); 
	COMM_WORLD.Send(sendVarBuffer, n*lenvar, MPI_DOUBLE, 0, tag); 
//printf(" sendVarBuffer[0]=%g sendVarBuffer[n*lenvar-1]=%g \n", sendVarBuffer[0],sendVarBuffer[n*lenvar-1]);
	delete sendVarBuffer;
}
#endif

/************************************************/

//?void	Domain::run()
//?{
//?	//for (int i=0; i<ndomains; i++)
//?	// Process the domains
//?	{	//Pressure solver
//?		int	type=cells;
//?		for (int ie=0; ie<Ne[type]; ie++)
//?		{	
//?
//?			
//?		}
//?	}		
//?}
void	Domain::getGridLimits
(
	double	ymin[],
	double	ymax[]
)
{	int	ne=Ne[nodes];
	double	*X=coordinates[nodes].val;
	if (type==dynamic) 
	{	if (!fixedlimits)
		{	for (int i=0; i<DIM; i++)
			{	xmin[i]=+LARGE;
				xmax[i]=-LARGE;
			}
			DNode	*node=dnode_root;
			do
			{	double	*x=node->x;
				for (int i=0; i<DIM; i++)
					{	if (xmin[i]>x[i]) xmin[i]=x[i];
					if (xmax[i]<x[i]) xmax[i]=x[i];
				}
				node=node->next;
			}	while(node!=dnode_root);
		}
		for (int i=0; i<DIM; i++)
		{	ymin[i]=xmin[i];
			ymax[i]=xmax[i];
		}
		return;
	}
	for (int i=0; i<DIM; i++)
	{	ymin[i]=LARGE; ymax[i]=-LARGE;
	}
	for (int i=0; i<ne; i++)
	{	double	*x=X+DIM*i;
		for (int j=0; j<DIM; j++)
		{	if (ymin[j]>x[j])ymin[j]=x[j];
			if (ymax[j]<x[j])ymax[j]=x[j];
		}
	}
	for (int i=0; i<DIM; i++)
	{	xmin[i]=ymin[i];
		xmax[i]=ymax[i];
	}
}
template <class Anything>
void	Domain::swap(Anything *&a, Anything *&b)
{
	int	size=sizeof(Anything);
	Anything c;
	memcpy(&c,a,size);
	memcpy(a,b,size);
	memcpy(b,&c,size);
//	c=a;a=b;b=c;
}
//?void	Domain::swapCells
//?(	//Swap two cells in the grid
//?	int	ia, 
//?	int	ib
//?)
//?{	//Face-cell and face-face connectivities
//?	// will be corrupt after this procedure
//?	// and will have to be regenerated by
//?	// setFaces();
//?	static struct Cell	tmp;
//?	int	cellsize=sizeof(struct Cell),ic[2];
//?	struct Neighbor
//?	{	int
//?			ind,
//?			face,
//?			cell;
//?	}	neighbor[2][Nf],*nb;
//?	struct Cell	*a[2];
//?	ic[0]=ia;ic[1]=ib;
//?	for (int i=0;i<2;i++)a[i]=cell+ic[i];
//?	//Reconnect the neighbors
//?	//NOTE: neighbor face index, vv is not updated
//?	for (int i=0; i<Nf; i++)
//?	{	for (int k=0;k<2;k++)
//?		{	int	inb;
//?			nb=neighbor[k]+i;
//?			if ((inb=a[k]->cell[i])>=0)
//?			{	for (int j=0; j<Nf; j++)
//?					if (cell[inb].cell[j]==ic[k]) 
//?					{	nb->ind=inb;
//?						nb->face=j;
//?						nb->cell=ic[(k+1)%2];
//?						break;
//?					}
//?			}
//?			else	nb->ind=-1; 
//?		}
//?	}
//?	for (int k=0; k<2; k++)
//?	for (int i=0; i<Nf; i++)
//?	{	nb=neighbor[k]+i;
//?		if (nb->ind<0) continue;
//?		cell[nb->ind].cell[nb->face]=nb->cell;
//?	}
//?	swap(a[0],a[1]);
//?}
void	Domain::setElementColor(int element_type, float rgbcolor[])
{
	if (element_type<maxelements)
	for (int irgb=0; irgb<3; irgb++)
		disp[element_type].rgbcolor[irgb]=rgbcolor[irgb];
}
void	Domain::getElementColor(int element_type, float rgbcolor[])
{
	if (element_type<maxelements)
		for (int irgb=0; irgb<3; irgb++)
			rgbcolor[irgb]=disp[element_type].rgbcolor[irgb];
}
void	Domain::setDisplayVariable(int ivar)
{
	dispvar=ivar;
}
void	Domain::setDisplayVariableComp(int icomp)
{
	dispvarcomp=icomp;
}
int	Domain::getLighting          (int type){return disp[type].lighting;}
int	Domain::getDisplayVariable    (){return dispvar        ;}
int	Domain::getDisplayVariableComp(){return dispvarcomp    ;}
int	Domain::getNoBoundaryVertexes (){return nbv            ;}
int	Domain::getNoBoundaryFaces    (){return nbf            ;}
int	Domain::getNoBoundaryCells    (){return nbc            ;}
int	Domain::getNoVariables        (){return nvar           ;}
int	Domain::getNoVariables(Elements vartype)
{
	int	n=0;
	for (int i=0; i<nvar; i++)
		if (variable[i].type==vartype) n++;
	return n;
}
int	Domain::getNoVariables(Elements vartype, int dim)
{
	int	n=0;
	for (int i=0; i<n; i++)
		if (variable[i].type==vartype&&variable[i].dimension==dim) n++;
	return n;
}
int	Domain::getVarBufSize(Elements vartype)
{	int n=0;
	if(variable!=NULL)
	{	
		for (int i=0; i<nvar; i++)
		if (variable[i].type==vartype)
			n+=variable[i].dimension;
	}
	return n;
}
///	void	Domain::loadGrid(char *filename)
///	{	int	nn,nc,//number of nodes and cells
///			gotnodes=0,
///			gotcells=0,
///			gotnncon=0,
///			gotcncon=0,
///			gotcccon=0;
///		char	word[MAXLINLEN];
///		FILE	*inp;
///		deleteRing(dnode_root);
///		deleteRing(dcell_root);
///		nn=nc=0;
///		OPENREAD(filename,inp);
///		GETWORD(word,inp);
///		while (!feof(inp))
///		{	if(!gotnodes&&strcmp(word,"NODES")==0)
///			{	DNode	*current;
///				if (option.verbose|option.debug)
///				{	printf("Reading nodes ...");FLUSH;}
///				while (!feof(inp))
///				{	int i,type;
///					double	x[DIM];
///					for (i=0; i<DIM; i++)
///					{	GETWORD(word,inp);
///						if (isalpha(*word))break;
///						sscanf(word,"%lg",x+i);
///					}
///					if (i<DIM)
///					{	if (i==0) break;
///						fprintf
///						(	stderr,
///							"ERROR IN %s: SAW A WORD '%s' WHEN A NUMBER WAS EXPECTED\n",
///							filename, word
///						);	exit(1);
///					}
///					if (dnode_root==NULL)
///					{	dnode_root=new DNode;
///						current=dnode_root;
///						current->next=current->prev=current;
///					}
///					else
///					{	DNode	*old=current;
///						current->next=new DNode;
///						current=current->next;
///						current->prev=old;
///						current->next=dnode_root;
///						dnode_root->prev=current;
///					}
///					///current->nneib=0;
///					current->neib=NULL;
///					for (i=0; i<DIM; i++)
///						current->x[i]=x[i];
///					GETWORD(word,inp);
///					sscanf(word,"%d",&type);
///					if (type==1)current->type.boundary=1;
///					else	current->type.boundary=0;
///					current->type.insidetool=0;
///					nn++;
///				}
///				if (option.verbose|option.debug)
///				{	printf(" %d node coordinates read\n",nn);FLUSH;}
///				gotnodes=1;
///			}
///			else
///			if(!gotnncon&&strcmp(word,"NODE-NODE-CONNECTIVITY")==0)
///			{	DNode
///					*node,
///					**NP;//node-pointer array
///				if (nn==0)
///				{	fprintf
///					(	stderr,
///						"ERROR IN %s: CAN'T READ NODE CONNECTIVITY BEFORE NODES ARE DEFINED\n",
///						filename
///					);	exit(1);
///				}
///				//Create a node pointer array
///				if (option.debug|option.verbose)
///				{	printf("Creating node-node connectivity for %d nodes ... ",nn);FLUSH;}
///				ALLOC(NP,DNode*,nn);
///				node=dnode_root;
///				for (int i=0; i<nn; i++)
///				{
///					NP[i]=node;
///					node=node->next;
///				}
///				node=dnode_root;
///				for (int i=0; i<nn&&!feof(inp); i++)
///				{	int	nneib;
///					DNodeList	*neib;
///					GETWORD(word,inp);
///					if (isalpha(*word))
///					{	fprintf
///						(	stderr,
///							"ERROR IN %s: SAW A WORD '%s' WHEN A NUMBER WAS EXPECTED\n",
///							filename,word
///						);	exit(1);
///					}
///					sscanf(word,"%d",&nneib);
///	#ifdef WTC
///					//-ALLOC(node->neib,DNode*,nneib);
///					//-ALLOC(node->dist,double,nneib);
///					//-ALLOC(node->e,double,nneib*DIM);
///					//-node->nneib=nneib;
///	#endif
///					for (int j=0; j<nneib; j++)
///					{	int	ineib;
///						GETWORD(word,inp);
///						sscanf(word,"%d",&ineib);
///						node->neib[j]=NP[ineib];
///					}
///					node=node->next;
///				}
///				delete(NP);
///				gotnncon=1;
///				if (option.verbose|option.debug)
///				{	printf("done\n");FLUSH;}
///			}
///			else
///			if(!gotcells&&strcmp(word,"CELLS")==0)
///			{	while (!feof(inp))
///				{	GETWORD(word,inp);
///					if (isalpha(*word))break;
///	
///					nc++;
///				}
///				gotcells=1;
///			}
///			else
///			if(!gotcncon&&strcmp(word,"CELL-NODE-CONNECTIVITY")==0)
///			{	if (nn=0||nc==0)
///				{	fprintf
///					(	stderr,
///						"ERROR IN %s: CAN'T READ CELL-NODE CONNECTIVITY BEFORE NODES OR CELLS ARE DEFINED\n",
///						filename
///					);	exit(1);
///				}
///				while (!feof(inp))
///				{	GETWORD(word,inp);
///					if (isalpha(*word))break;
///				}
///				gotcncon=1;
///			}
///			else
///			if(!gotcccon&&strcmp(word,"CELL-CELL-CONNECTIVITY")==0)
///			{	if (nc=0)
///				{	fprintf
///					(	stderr,
///						"ERROR IN %s: CAN'T READ CELL-CELL CONNECTIVITY BEFORE CELLS ARE DEFINED\n",
///						filename
///					);	exit(1);
///				}
///				while (!feof(inp))
///				{	GETWORD(word,inp);
///					if (isalpha(*word))break;
///				}
///				gotcccon=1;
///			}
///			else	break;
///		}
///		fclose(inp);
///	}
void	Domain::printVariables()
{
	for (int i=0; i<nvar; i++)
		printf
		(
			"%d. %s: type=%s, dim=%d, no.elements=%d\n",
			i+1,variable[i].name,
			elementype[variable[i].type],
			variable[i].dimension,variable[i].size
		);
}
int	Domain::ihostcell
(	double *z, //z[DIM] coordinates of a point 
	DCell *cell//initial cell
)
{//If x is inside cell returns -1
 //otherwise returns local index of the neighbor cell in the direction to x
	DNode	**vert;
	int	ineib=-1;
	double	a,b0,b1,*x,*y,
		b[DIM],c[DIM],d[Nv][DIM];
	if (cell==NULL) return -1;
	vert=cell->vert;
	x=vert[0]->x;
	for (int i=0; i<Nv; i++)
	{
	/*	Find edges of the tetrahedral sub-lement:
	 *	d[k]=edge-vectors of a sub-element
	 */
		y=vert[(i+1)%Nv]->x;
		for (int j=0; j<DIM; j++)
			d[i][j]=y[j]-x[j];
		x=y;
	}
	//Loop throught the vertexes
	for(int i=0;i<Nv;i++)
	{	int	i1=(i+1)%Nv,i2=(i+2)%Nv;
		y=vert[i1]->x;
		for(int j=0;j<DIM;j++)c[j]=y[j]-z[j];
		VECP(b,d[i1],d[i2]);
		b0=SCLP(b,c);/* distance of Z from the face formed by d[i1],d[i2] */
		b1=SCLP(b,d[i]);
		if (b0*b1<0.0)
			return i;/* point is outside of the tetrahedron 
		                           on the side of face i */
	}
	return -1; /* point is inside */
}
//?Cell	*Domain::hostcell
//?(	double *z, //z[DIM] coordinates of a point 
//?	Cell *cell//initial cell
//?)
//?{//If x is inside cell returns cell
//? //otherwise returns neighbor cell closer to x
//?	int	*vert,*neib;
//?	double	*X=coordinates[nodes].val;
//?	register	double	a,b0,b1,*x,*y;
//?	double	b[DIM],c[DIM],d[Nv][DIM];
//?	if (cell==NULL) return NULL;
//?	vert=cell->vert;
//?	neib=cell->cell;
//?	x=X+DIM*vert[0];
//?	for (int i=0; i<Nv; i++)
//?	{
//?	/*	Find edges of the tetrahedral sub-lement:
//?	 *	d[k]=edge-vectors of a sub-element
//?	 */
//?		y=X+DIM*vert[(i+1)%Nv];
//?		for (int j=0; j<DIM; j++)
//?			d[i][j]=y[j]-x[j];
//?		x=y;
//?	}
//?	//Loop throught the vertexes
//?	for(int i=0;i<Nv;i++)
//?	{	int	i1=(i+1)%Nv,i2=(i+2)%Nv;
//?		y=X+DIM*vert[i1];
//?		for(int j=0;j<DIM;j++)c[j]=y[j]-z[j];
//?		VECP(b,d[i1],d[i2]);
//?		b0=SCLP(b,c);/* distance of Z from the face formed by d[i1],d[i2] */
//?		b1=SCLP(b,d[i]);
//?		if (b0*b1<0.0)
//?			return cell+neib[i];/* point is outside of the tetrahedron 
//?		                           on the side of face i */
//?	}
//?	return cell; /* point is inside */
//?}
DCell	*Domain::hostcell
(	double *z, //z[DIM] coordinates of a point 
	DCell *cell//initial cell
)
{//If x is inside cell returns cell
 //otherwise returns neighbor cell closer to x
	DNode	**vert;
	DCell	**neib;
	register	double	a,b0,b1,*x,*y;
	double	b[DIM],c[DIM],d[Nv][DIM];
	if (cell==NULL) return NULL;
	vert=cell->vert;
	neib=(DCell**)cell->neib;
	x=vert[0]->x;
	for (int i=0; i<Nv; i++)
	{
	/*	Find edges of the tetrahedral sub-lement:
	 *	d[k]=edge-vectors of a sub-element
	 */
		y=vert[(i+1)%Nv]->x;
		for (int j=0; j<DIM; j++)
			d[i][j]=y[j]-x[j];
		x=y;
	}
	//Loop throught the vertexes
	for(int i=0;i<Nv;i++)
	{	int	i1=(i+1)%Nv,i2=(i+2)%Nv;
		y=vert[i1]->x;
		for(int j=0;j<DIM;j++)c[j]=y[j]-z[j];
		VECP(b,d[i1],d[i2]);
		b0=SCLP(b,c);/* distance of Z from the face formed by d[i1],d[i2] */
		b1=SCLP(b,d[i]);
		if (b0*b1<0.0)
			return neib[i];/* point is outside of the tetrahedron 
		                           on the side of face i */
	}
	return cell; /* point is inside */
}
//?ElementStatus	Domain::hostcell//returns either bc or a face number
//?(
//?	double *z, //z[DIM] coordinates of a point 
//?	Cell *current_cell,//initial cell
//?	Cell	*&newcell,
//?	int	*iface//cell face number
//?)
//?{//If x is inside cell returns cell
//? //otherwise returns neighbor cell closer to x
//?	int	*vert,*facetype=current_cell->cell,*neib;
//?	ElementStatus	bc;
//?	double	a,b0,b1,*x,*y,
//?		b[DIM],c[DIM],d[Nv][DIM],
//?		*X=coordinates[nodes].val;
//?	Cell	*root_cell=cell;
//?	if (root_cell==NULL) {newcell=NULL; *iface=-1; return dead;}
//?	vert=current_cell->vert;
//?	neib=current_cell->cell;
//?	x=X+DIM*vert[0];
//?	for (int i=0; i<Nv; i++)
//?	{
//?	/*	Find edges of the tetrahedral sub-lement:
//?	 *	d[k]=edge-vectors of a sub-element
//?	 */
//?	//	y=X+DIM*v[(i+1)%Nv];
//?		y=X+DIM*vert[(i+1)%Nv];
//?		for (int j=0; j<DIM; j++)
//?			d[i][j]=y[j]-x[j];
//?		x=y;
//?	}
//?	//Loop throught the vertexes
//?	for(int i=0;i<Nv;i++)
//?	{	int	i1=(i+1)%Nv,i2=(i+2)%Nv;
//?		y=X+DIM*vert[i1];
//?		for(int j=0;j<DIM;j++)c[j]=y[j]-z[j];
//?		VECP(b,d[i1],d[i2]);
//?		b0=SCLP(b,c);/* distance of Z from the face formed by d[i1],d[i2] */
//?		b1=SCLP(b,d[i]);
//?		if (b0*b1<0.0)
//?		{	newcell=root_cell+neib[i];
//?			*iface=i;
//?			return facetype[i]<0?(ElementStatus)facetype[i]:internal;
//?		}
//?	}
//?	newcell=current_cell;
//?	*iface=-1;
//?	return internal; /* point is inside */
//?}
///	ElementStatus	Domain::hostcell//returns either bc or a face number
///	(
///		double *z, //z[DIM] coordinates of a point 
///		DCell *cell,//initial cell
///		DCell	*&newcell,
///		int	*iface//cell face number
///	)
///	{//If x is inside cell returns cell
///	 //otherwise returns neighbor cell closer to x
///		DNode	**vert;
///		DCell	**neib;
///		ElementStatus	bc,*facetype=cell->facetype;
///		register	double	a,b0,b1,*x,*y;
///		double	b[DIM],c[DIM],d[Nv][DIM];
///		if (cell==NULL) {newcell=NULL; *iface=-1; return dead;}
///		vert=cell->vert;
///		neib=(DCell**)cell->neib;
///		x=vert[0]->x;
///		for (int i=0; i<Nv; i++)
///		{
///		/*	Find edges of the tetrahedral sub-lement:
///		 *	d[k]=edge-vectors of a sub-element
///		 */
///		//	y=X+DIM*v[(i+1)%Nv];
///			y=vert[(i+1)%Nv]->x;
///			for (int j=0; j<DIM; j++)
///				d[i][j]=y[j]-x[j];
///			x=y;
///		}
///		//Loop throught the vertexes
///		for(int i=0;i<Nv;i++)
///		{	int	i1=(i+1)%Nv,i2=(i+2)%Nv;
///			y=vert[i1]->x;
///			for(int j=0;j<DIM;j++)c[j]=y[j]-z[j];
///			VECP(b,d[i1],d[i2]);
///			b0=SCLP(b,c);/* distance of Z from the face formed by d[i1],d[i2] */
///			b1=SCLP(b,d[i]);
///			if (b0*b1<0.0)
///			{	newcell=neib[i];
///				*iface=i;
///				return facetype[i];
///			}
///		}
///		newcell=cell;
///		*iface=-1;
///		return internal; /* point is inside */
///	}
ElementStatus	Domain::hostcell//returns either bc or a face number
(
	double *z, //z[DIM] coordinates of a point 
	DCell *cell,//initial cell
	DCell	*&newcell,
	double	*norm, //norm[DIM]: face-normal vector
	int	*iface//cell face number
)
{//If x is inside cell returns cell
 //otherwise returns neighbor cell closer to x
	DNode	**vert;
	DCell	**neib;
	ElementStatus	bc,*facetype=cell->facetype;
	register	double	a,b0,b1,*x,*y;
	double	b[DIM],c[DIM],d[Nv][DIM];
	if (cell==NULL) {newcell=NULL; *iface=-1; return dead;}
	vert=cell->vert;
	neib=(DCell**)cell->neib;
	x=vert[0]->x;
	for (int i=0; i<Nv; i++)
	{
	/*	Find edges of the tetrahedral sub-lement:
	 *	d[k]=edge-vectors of a sub-element
	 */
	//	y=X+DIM*v[(i+1)%Nv];
		y=vert[(i+1)%Nv]->x;
		for (int j=0; j<DIM; j++)
			d[i][j]=y[j]-x[j];
		x=y;
	}
	//Loop throught the vertexes
	for(int i=0;i<Nv;i++)
	{	int	i1=(i+1)%Nv,i2=(i+2)%Nv;
		y=vert[i1]->x;
		for(int j=0;j<DIM;j++)c[j]=y[j]-z[j];
		VECP(b,d[i1],d[i2]);
		b0=SCLP(b,c);/* distance of Z from the face formed by d[i1],d[i2] */
		b1=SCLP(b,d[i]);
		if (b0*b1<0.0)
		{	newcell=neib[i];
			*iface=i;
			if(facetype[i]!=internal)
			{//return the boundary-normal vector
				b0=1./LENGTH(b);
				if(b1<0)b0=-b0;
				for(int i=0;i<DIM;i++)norm[i]=b0*b[i];
			}
			return facetype[i];
		}
	}
	newcell=cell;
	*iface=-1;
	return internal; /* point is inside */
}
ElementStatus	Domain::hostcell
(	double *z, //z[DIM] coordinates of a point 
	DCell *cell,//initial cell
	DCell	*&newcell
)
{//If x is inside cell returns cell
 //otherwise returns neighbor cell closer to x
	DNode	**vert;
	DCell	**neib;
	ElementStatus	bc,*facetype=cell->facetype;
	register	double	a,b0,b1,*x,*y;
	double	b[DIM],c[DIM],d[Nv][DIM];
	if (cell==NULL) {newcell=NULL; return dead;}
	vert=cell->vert;
	neib=(DCell**)cell->neib;
	x=vert[0]->x;
	for (int i=0; i<Nv; i++)
	{
	/*	Find edges of the tetrahedral sub-lement:
	 *	d[k]=edge-vectors of a sub-element
	 */
	//	y=X+DIM*v[(i+1)%Nv];
		y=vert[(i+1)%Nv]->x;
		for (int j=0; j<DIM; j++)
			d[i][j]=y[j]-x[j];
		x=y;
	}
	//Loop throught the vertexes
	for(int i=0;i<Nv;i++)
	{	int	i1=(i+1)%Nv,i2=(i+2)%Nv;
		y=vert[i1]->x;
		for(int j=0;j<DIM;j++)c[j]=y[j]-z[j];
		VECP(b,d[i1],d[i2]);
		b0=SCLP(b,c);/* distance of Z from the face formed by d[i1],d[i2] */
		b1=SCLP(b,d[i]);
		if (b0*b1<0.0)
		{	newcell=neib[i];
			return facetype[i];
//				if (newcell!=NULL)
//					return newcell->facetype[i];/* point is outside of the tetrahedron 
//																 on the side of face i */
//				else
//					return dead;
		}
	}
	newcell=cell;
	return internal; /* point is inside */
}
//?Cell	*Domain::findcell
//?(	double *x,//x[DIM] coordinates of a point 
//?	Cell *cell//current cell
//?)
//?{//Locates cell containing point x
//? // works only if 'x' is visible from 'cell'
//?	int	counter=0;
//?	Cell	*next,*prev;
//?	next=prev=cell;
//?	while ((next=hostcell(x,cell))!=NULL)
//?	{	if (next==cell) return cell;
//?		if (next==prev) return NULL;
//?		prev=cell;
//?		cell=next;//To escape from endless looping
//?		if (++counter>MAXGRIDCELLS) 
//?		{	fprintf
//?			(	stderr,
//?				"WARNING: Endless loop in locating the cell\n"
//?			);
//?			return NULL;
//?		}
//?	}
//?	return NULL;
//?}
DCell	*Domain::findcell
(	double *x,//x[DIM] coordinates of a point 
	DCell *cell//current cell
)
{//Locates cell containing point x
 // works only if 'x' is visible from 'cell'
	int	counter=0;
	DCell	*next,*prev;
	next=prev=cell;
	while ((next=hostcell(x,cell))!=NULL)
	{	if (next==cell) return cell;
		if (next==prev) return NULL;
		prev=cell;
		cell=next;//To escape from endless looping
		if (++counter>MAXGRIDCELLS) 
		{	fprintf
			(	stderr,
				"WARNING: Endless loop in locating the cell\n"
			);
			return NULL;
		}
	}
	return NULL;
}
DCell	*Domain::findcell
(	double *x//x[DIM] coordinates of a point 
)
{
	return	findcell(x,dcell_root);
}
ElementStatus	Domain::findcell
(	double *x,//x[DIM] coordinates of a point 
	DCell *cell,//current cell
	DCell	*&newcell //new cell
)
{//Locates cell containing point x
 // works only if 'x' is visible from 'cell'
	int	counter=0;
	ElementStatus	bc;
	DCell	*next,*prev;
	next=prev=cell;
	while ((bc=hostcell(x,cell,next))==internal)
	{	if (next==cell) {newcell=cell;return bc;}
		if (next==prev) 
		{	fprintf
			(	stderr,
				"WARNING: Cell loop on themselves while locating the cell\n"
			);
			newcell=NULL;
			return dead;
		}
		prev=cell;
		cell=next;
		if (++counter>MAXGRIDCELLS) 
		{	fprintf
			(	stderr,
				"WARNING: Loop iterations exceeded %d in locating the cell\n",
				MAXGRIDCELLS
			);
			newcell=NULL;
			return dead;
		}
	}
	newcell=NULL;
	return bc;
}
///	ElementStatus	Domain::findcell//returns both cell and face number
///	(	double *x,//x[DIM] coordinates of a point 
///		DCell *cell,//current cell
///		DCell	*&newcell, //new cell
///		int	*iface //-1 if inside the cell, 0:Nf if outside the cell
///	)
///	{//Locates cell containing point x
///	 // works only if 'x' is visible from 'cell'
///		int	counter=0;
///		ElementStatus	bc;
///		DCell	*next,*prev;
///		next=prev=cell;
///		while ((bc=hostcell(x,cell,next,iface))==internal)
///		{
///			if (next==cell) {newcell=cell;return bc;}
///			if (next==prev) 
///			{	fprintf
///				(	stderr,
///					"WARNING: Cell loop on themselves while locating the cell\n"
///				);
///				newcell=NULL;
///				*iface=-1;
///				return dead;
///			}
///			prev=cell;
///			cell=next;
///			if (++counter>MAXGRIDCELLS) 
///			{	fprintf
///				(	stderr,
///					"WARNING: Loop iterations exceeded %d in locating the cell\n",
///					MAXGRIDCELLS
///				);
///				newcell=NULL;
///				*iface=-1;
///				return dead;
///			}
///		}
///		newcell=NULL;
///		return bc;
///	}
ElementStatus	Domain::findcell
(	double *x,//x[DIM] coordinates of a point 
	DCell *cell,//current cell
	DCell	*&newcell, //new cell
	double	*norm //norm[DIM]: face-normal vector
)
{//Locates cell containing point x
 // works only if 'x' is "visible" from 'cell'
	int	counter=0,iface;
	ElementStatus	bc;
	DCell	*next,*prev;
	next=prev=cell;
	while ((bc=hostcell(x,cell,next,norm,&iface))==internal)
	{	if (next==cell) {newcell=cell;return bc;}
		if (next==prev) 
		{	fprintf
			(	stderr,
				"WARNING: Cell loop on themselves while locating the cell\n"
			);
			newcell=NULL;
			return dead;
		}
		prev=cell;
		cell=next;
		if (++counter>MAXGRIDCELLS) 
		{	fprintf
			(	stderr,
				"WARNING: Loop iterations exceeded %d in locating the cell\n",
				MAXGRIDCELLS
			);
			newcell=NULL;
			return dead;
		}
	}
	newcell=NULL;
	return bc;
}
//?ElementStatus	Domain::findcell//returns both cell and face number
//?(	double *x,//x[DIM] coordinates of a point 
//?	Cell *cell,//current cell
//?	Cell	*&newcell, //new cell
//?	int	*iface //-1 if inside the cell, 0:Nf if outside the cell
//?)
//?{//Locates cell containing point x
//? // works only if 'x' is visible from 'cell'
//?	int	counter=0;
//?	ElementStatus	bc;
//?	Cell	*next,*prev;
//?	next=prev=cell;
//?	while ((bc=hostcell(x,cell,next,iface))==internal)
//?	{
//?		if (next==cell) {newcell=cell;return bc;}
//?		if (next==prev) 
//?		{	fprintf
//?			(	stderr,
//?				"WARNING: Cell loop on themselves while locating the cell\n"
//?			);
//?			newcell=NULL;
//?			*iface=-1;
//?			return dead;
//?		}
//?		prev=cell;
//?		cell=next;
//?		if (++counter>MAXGRIDCELLS) 
//?		{	fprintf
//?			(	stderr,
//?				"WARNING: Loop iterations exceeded %d in locating the cell\n",
//?				MAXGRIDCELLS
//?			);
//?			newcell=NULL;
//?			*iface=-1;
//?			return dead;
//?		}
//?	}
//?	newcell=NULL;
//?	return bc;
//?}
void	Domain::createBoundaryNodeList(BNodeList *&root)
{
	if (dnode_root!=NULL)
	{	DNode	*node=dnode_root;
		BNodeList	*bp;
		if (root!=NULL) deleteRing(root);
		do
		{	if (node->state.boundary)
			{//Create a new cell-face pointer
				if (root==NULL)
				{	root=new BNodeList;
					bp=root;
				}
				else
				{	bp->next=new BNodeList;
					bp=bp->next;
				}
				bp->node=node;
			}
			node=node->next;
		}	while(node!=dnode_root);
		if (root!=NULL)bp->next=root;
	}
}
void	Domain::createBoundaryFaceList(BFaceList *&root)
{	int nbf=0;	
	if (dcell_root!=NULL)
	{	DCell	*cell=dcell_root;
		BFaceList	*bp;
		do
		{	DNode **vert=cell->vert;
			DCell	**neib=(DCell**)cell->neib;
			ElementStatus	*facetype=cell->facetype;
			for (int iv=0; iv<Nv; iv++)
			{	if (neib[iv]==NULL)
				//if (facetype[iv]==boundary||neib[iv]==NULL)
				//if (facetype[iv]!=internal)
				{//Create a new cell-face pointer
					int
						i1=(iv+1)%Nv,
						i2=(iv+2)%Nv,
						i3=(iv+3)%Nv;
					double	a[DIM],area,
						e0[DIM],e1[DIM],e2[DIM];//two edges of the face
					if (root==NULL)
					{	root=new BFaceList;
						bp=root;
					}
					else
					{	bp->next=new BFaceList;
						bp=bp->next;
					}
					bp->cell=cell;
					bp->iface=iv;
					//Face center
					for (int i=0; i<DIM; i++)
					{	double	a=0.0;
						for (int j=0; j<Nfv; j++)
						{	int	k=(iv+j+1)%Nv;
							a+=vert[k]->x[i];
						}
						bp->x[i]=a/(double)Nfv;
						e0[i]=vert[i1]->x[i]-vert[iv]->x[i];
						e1[i]=vert[i2]->x[i]-vert[i1]->x[i];
						e2[i]=vert[i3]->x[i]-vert[i2]->x[i];
					}
					//Face area vector
					VECP(a,e1,e2);
					if (SCLP(e0,a)<0.0)
						for (double *p=a; p-a<DIM; p++) *p*=-.5;
					else
						for (double *p=a; p-a<DIM; p++) *p*= .5;
					area=LENGTH(a);
					for (int i=0; i<DIM; i++) bp->norm[i]=a[i]/area;
					bp->area=area;
					nbf++;
				}
			}
			cell=cell->next;
		}	while(cell!=dcell_root);
		bp->next=root;
	}
}
void	Domain::updateBoundaryFaceList(BFaceList *&root)
{
	if (root!=NULL)
	{	BFaceList	*bface=root;
		do
		{	DCell
				*cell=bface->cell,
				**neib=(DCell**)cell->neib;
			DNode **vert=cell->vert;
			ElementStatus	*facetype=cell->facetype;
			int
				iv=bface->iface,
				i1=(iv+1)%Nv,
				i2=(iv+2)%Nv,
				i3=(iv+3)%Nv;
			double	a[DIM],area,areai,
				e0[DIM],e1[DIM],e2[DIM];//two edges of the face
			//Face center
			for (int i=0; i<DIM; i++)
			{	double	a=0.0;
				for (int j=0; j<Nfv; j++)
				{	int	k=(iv+j+1)%Nv;
					a+=vert[k]->x[i];
				}
				bface->x[i]=a/(double)Nfv;
				e0[i]=vert[i1]->x[i]-vert[iv]->x[i];
				e1[i]=vert[i2]->x[i]-vert[i1]->x[i];
				e2[i]=vert[i3]->x[i]-vert[i2]->x[i];
			}
			//Face area vector
			VECP(a,e1,e2);
			if (SCLP(e0,a)<0.0)
				for (double *p=a; p-a<DIM; p++) *p*=-.5;
			else
				for (double *p=a; p-a<DIM; p++) *p*= .5;
			area=LENGTH(a);
			areai=1./area;
			for (int i=0; i<DIM; i++) bface->norm[i]=a[i]*areai;
			bface->area=area;
			bface=bface->next;
		}	while(bface!=root);
	}
}
void	Domain::createBoundaryCellList(DCellList *&root)
{
	if (dcell_root!=NULL)
	{	DCell	*cell=dcell_root;
		DCellList	*bp;
		deleteRing(root);
		do
		{	DNode **vert=cell->vert;
			for (int iv=0; iv<Nv; iv++)
			{	if (vert[iv]->state.boundary)
				{//Create a new cell-face pointer
					if (root==NULL)
					{	root=new DCellList;
						bp=root;
					}
					else
					{
						bp->next=new DCellList;
						bp=bp->next;
					}
					bp->cell=cell;
					break;
				}
			}
			cell=cell->next;
		}	while(cell!=dcell_root);
		bp->next=root;
	}
}
void	Domain::connectBoundaryCellsToFaces(BFaceList *&root)
{
	BFaceList	*bface,*prev;
	if (root==NULL)return;
	bface=root;
	do
	{	DCell	*cell=bface->cell;
		int	iface=bface->iface;
		ElementStatus	facetype=cell->facetype[iface];
		if(facetype>=boundary)//wall
		{	cell->neib[iface]=(BFaceList *)bface;
			cell->facetype[iface]=composite;
		}
		bface->type=cell->facetype[iface];
		prev=bface;
		bface=bface->next;
	}	while(bface!=root);
}
void	Domain::setBC()
{
int i=0;
	DNode *node_root=dnode_root;
	DCell *cell_root=dcell_root;	
	if(cell_root!=NULL)
	{	DCell	*cell=cell_root;
		do
		{	void	**neibs=cell->neib;
			DNode	**verts=cell->vert;
			ElementStatus	*facetype=cell->facetype;
			for(int iv=0;iv<Nv;iv++)
			if(neibs[iv]==NULL)///||facetype[iv]==composite)
			{	int	jv,bndflag;
				for(jv=0;jv<Nfv;jv++)
					if((bndflag=verts[(iv+jv+1)%Nv]->type)>=boundary)break;//belongs to the wall
//					if(jv==Nfv)//belongs to in/out-let
//					{	facetype[iv]=(ElementStatus)bndflag;
//						neibs[iv]=NULL;
//					}
//					else//belongs to the boundary
//					{	if(facetype[iv]==composite)
//						{	BFaceList	*face=(BFaceList*)neibs[iv];
//							face->type=(ElementStatus)bndflag;
//						}
//						else
						facetype[iv]=(ElementStatus)bndflag;//cylinder wall wall
//					}
			}
			cell=cell->next;
		}	while(cell!=cell_root);
	}
///		if(cell_root!=NULL)
///		{	DCell	*cell=cell_root;
///			do
///			{	DCell	**neibs=(DCell**)cell->neib;
///				DNode	**verts=cell->vert;
///				for(int iv=0;iv<Nv;iv++)
///				{	ElementStatus	*facetype=cell->facetype;
///					if(neibs[iv]==NULL)
///					{	ElementStatus	type=facetype[iv];
///						if(type>=boundary)//cylinder wall
///						for(int jv=0;jv<Nfv;jv++)
///							verts[(iv+jv+1)%Nv]->type=type;
///					}
///	//					if(facetype[iv]==composite)
///	//					{	BFaceList	*bface=(BFaceList*)neibs[iv];
///	//						ElementStatus	type=bface->type;
///	//						if(type>=boundary)//cylinder wall
///	//						for(int jv=0;jv<Nfv;jv++)
///	//							verts[(iv+jv+1)%Nv]->type=type;
///	//					}
///				}
///				cell=cell->next;
///			}	while(cell!=cell_root);
///		}
}
template	<class List>
void	Domain::deleteRing(List *&element)
{
	List	*current;
	if (element==NULL) return;
	current=element->next;
	while (current!=element)
	{	List	*kill=current;
		current=current->next;
		delete kill;
	}
	delete element;
	element=NULL;
}
///void	Domain::deleteRing(BFaceList *&element)
///{
///	BFaceList *current;
///	if (element==NULL) return;
///	current=element->next;
///	while (current!=element)
///	{	BFaceList	*kill=current;
///		current=current->next;
///		delete kill;
///	}
///	delete element;
///	element=NULL;
///}
void	Domain::setCellVolumes
(
	int	ivar	
)
{	const	double	onesixth=1./6.;
	int	loc=variable[ivar].loc;
	DCell	*root=dcell_root,
		*cell=root;
	if (root==NULL) return;
	cell=root;
	do 
	{	DNode	**verts=cell->vert;
		double 
			e[Nv][DIM],//edges
			*var=cell->var,
			*vol=var+loc,
			volume;
		//Edges
		EDGES(e,verts);
		volume=0.0;
		for (int ivert=0; ivert<Nv; ivert++)
		{	DNode	*vert=verts[ivert];
			double
				dvol,
				x[DIM],//face-center vector
				a[DIM];///=A[ivert],///=areas+DIM*ivert;
			AREA2(ivert,e,a);
			//Face-center 
			for (int i=0; i<DIM; i++)
			{	double	xcenter=0.0;
				for (int jvert=0; jvert<DIM; jvert++)
					xcenter+=verts[(ivert+jvert+1)%Nv]->x[i];
				x[i]=xcenter/(double)Nfv;
			}
			//Volume flux
			dvol=onesixth*SCLP(a,x);
			volume+=dvol;//Node-volume
		}
#ifdef DEBUG
		if (volume<SMALL)
		{	fprintf(stderr,"SetCellVolumes: volume=%g is too small\n",volume);
			exit(1);
		}
#endif
		*vol=volume;
		cell=cell->next;
	}	while(cell!=root);
	if(option.verbose)
	{	double	totvol=0.0;
		cell=root;
		do 
		{	totvol+=cell->var[loc];
			cell=cell->next;
		}	while(cell!=root);
		printf("Cell volume: %g\n",totvol);
	}
}
void	Domain::setNodeVolumes
	//Computes nodal volumes, using
	//previously computed cell-volumes
(	int	ivarvoln,
	int	ivarvolc
)
{
	int
		ivoln=variable[ivarvoln].loc,
		ivolc=variable[ivarvolc].loc;
	DNode	*node=dnode_root;
	DCell	*cell=dcell_root;
	if (dnode_root==NULL) return;
	do//Zero out source terms
	{	double
			*varn=node->var,
			*voln=varn+ivoln;
		*voln=0.0;
		node=node->next;
	}	while(node!=dnode_root);
	if (dcell_root==NULL)return;
	cell=dcell_root;
	do 
	{	DNode	**verts=cell->vert;
		double	volc=cell->var[ivolc];
		//Face Areas and Face Fluxes
		for (int ivert=0; ivert<Nv; ivert++)
			verts[ivert]->var[ivoln]+=volc;
		cell=cell->next;
	}	while(cell!=dcell_root);
	if(option.verbose)
	{	double	totvol=0.0;
		node=dnode_root;
		do 
		{	double	voln=node->var[ivoln];
			totvol+=voln;
			node=node->next;
		}	while(node!=dnode_root);
		printf("Nodal volume: %g\n",totvol/4);
	}
}
void	Domain::setNodeVolumes
(	//Computes nodal volumes
	int	ivarvol
)
{	const double	onesixth=1./6.;
	int	ivol=variable[ivarvol].loc;
	DNode	*node=dnode_root;
	DCell	*cell=dcell_root;
	if (dnode_root==NULL) return;
	do//Zero out source terms
	{	double
			*var=node->var,
			*vol=var+ivol;
		*vol=0.0;
		node=node->next;
	}	while(node!=dnode_root);
	if (dcell_root==NULL)return;
	cell=dcell_root;
	do //Compute the volumes
	{	DNode	**verts=cell->vert;
		void	**neibs=cell->neib;
		double	
			e[Nv][DIM];//tetrahedral edge-vectors					
		EDGES(e,verts);
		//Volumes
		for (int iv=0; iv<Nv; iv++)
		{	double
				*var=verts[iv]->var,
				dvol,
				a[DIM],//face area vector
				xf[DIM];//face center coordinates
			AREA2(iv,e,a);
			for (int i=0; i<DIM; i++)
			{ //Face center
				double	xcenter=0.0;
				for (int j=0; j<Nfv; j++)
					xcenter+=verts[(iv+j+1)%Nv]->x[i];
				xf[i]=xcenter/(double)Nfv;
			}
			dvol=onesixth*(SCLP(a,xf));
			var[ivol]+=dvol;
			if(neibs[iv]==NULL)//Fix boundary volumes
			for(int j=0;j<Nv1;j++)
			{	int	jv=(iv+j+1)%Nv;
				double	*v=verts[jv]->var;
					v[ivol]+=dvol;
			}
		}
		cell=cell->next;
	}	while(cell!=dcell_root);
	if(option.verbose)
	{	double	totvol=0.0;
		node=dnode_root;
		do 
		{	double	vol=node->var[ivol];
			totvol+=vol;
			node=node->next;
		}	while(node!=dnode_root);
		printf("Nodal volume (single-step scheme): %g\n",totvol/4);
	}
}
void	Domain::indexNeibNodes()
{//For each neib cell sets cell's index variable
	//to contain the local index of this cell
	//as indexed from the neighbor cell
	DCell	*cell=dcell_root;
	if (cell==NULL) return;
	do
	{	DNode	**verts=cell->vert;
		DCell	**neibs=(DCell**)cell->neib;
		cell->index=0;
		for (int iv=0; iv<Nv; iv++)
		{	DCell	*neib=neibs[iv];
			if (neib!=NULL&&cell->facetype[iv]!=composite)
			{	int	jv;
				DCell	**neibneibs=(DCell**)neib->neib;
				for (jv=0; jv<Nv&&neibneibs[jv]!=cell; jv++);
				cell->index+=jv<<(iv<<1);
				//The index is then retrieved as
				//jvert=(int)(mask[ivert]&cell->index)>>((uint)ivert<<1);
#ifdef DEBUG
				if (jv==Nv)
				{	ERROR("ERROR in indexNeibNodes: Internal grid inconsistency\n");
				}
#endif
			}
			///else cell->index=-1;
		}
		cell=cell->next;
	}	while(cell!=dcell_root);
}
void	Domain::cleanInternalFaces()
{
	DCell	*cell=dcell_root;
	if (cell==NULL)return;
	do
	{	///DCell	**neibs=(DCell**)cell->neib;
		ElementStatus	*types=cell->facetype;
		for (int iv=0; iv<Nv; iv++)
		if (cell->neib[iv]!=NULL)
		{	if(cell->facetype[iv]!=composite)
				types[iv]=internal;
		}
		else
		if (types[iv]==internal)
			types[iv]=boundary;
		cell=cell->next;
	}	while(cell!=dcell_root);
}
//-	void	Domain::allocCellVars
//-	(
//-		int	bufsize,//variable buffer size
//-		int	bndsize,//boundary variable size
//-		DCell	*cell_root
//-	)
//-	{	if (cell_root!=NULL) 
//-		{	DCell	*cell=cell_root;
//-			if (option.verbose)
//-				printf("Allocating memory for dynamic variables domain type %d ...",type);
//-			do
//-			{	int	size;
//-				ElementStatus	*facetype=cell->facetype;
//-				size=bufsize;
//-				for (int i=0; i<Nv; i++)
//-				if (facetype[i]!=internal)
//-					size+=bndsize;
//-				cell->var=new double[size];
//-				cell=cell->next;
//-			}	while(cell!=cell_root);
//-			if (option.verbose)
//-				printf(" done\n");
//-		}
//-	}
template <class Element>
void	Domain::allocElementVars
(
	int	bufsize,
	Element	*element_root
)
{	if (element_root!=NULL) 
	{	Element	*element=element_root;
		if (option.verbose)
			printf("Allocating memory for dynamic variables domain type %d ...",type);
		do
		{	if(element->var!=NULL) free(element->var);
			element->var=new double[bufsize];
			element=element->next;
		}	while(element!=element_root);
		if (option.verbose)
			printf(" done\n");
	}
}
template <class Element>
void	Domain::deleteGrid(Element *&root)
{
		if (root!=NULL) 
		{	for
			(	Element	*node=root->next; 
				;node=node->next
			)
			{	Element	*prev=node;	
				delete prev;
			}
			delete root;
		}
}
void	Domain::relax
(
	int	ixold,
 	int	ivoln,
 	int	ivolc,
	double	relaxation_factor
)
{	static	const	double
		onethird=1.0/3.0,
		onesixth=1.0/6.0;
	//Relax volume nodes
	if (dnode_root!=NULL)
	{	DNode	*node=dnode_root;
	do//Initialize
	{	double
			*x=node->x,
			*var=node->var,
			*vol=var+ivoln,
			*xold=var+ixold;
		//Initialize to zero
		if(node->state.boundary==0)
		{	*vol=0.0;
			for (int i=0; i<DIM; i++) xold[i]=0.0;
		}
		node=node->next;
	}	while (node!=dnode_root);
	}
	if (dcell_root!=NULL)
	{	DCell	*cell=dcell_root;
	do//Assemble
	{	double	*x,*y,vol,
			*varc=cell->var,
			*volc=varc+ivolc,
			e[Nv][DIM];//tetrahedral edge-vectors					
		DNode	**vert=cell->vert;
		//Edges of the tetrahedron:
		x=vert[Nv1]->x;
		for (int iv=0; iv<Nv; iv++)
		{ double
		    *y=vert[iv]->x;
			for (int k=0; k<DIM; k++)
				e[iv][k]=y[k]-x[k];
			x=y;
		}
		//Compute volumes and cell-centers
		//Loop over the faces
		vol=0.0;//volume of the cell
		for (int iv=0; iv<Nv; iv++)
		{	int
				iv1=(iv+2)%Nv,
				iv2=(iv+3)%Nv;
			DNode	*node=vert[iv];
			double
				*varn=node->var,//nodal variables
				*voln=varn+ivoln,//nodal volume
				dv,//volume increment
				xf[DIM],//face center coordinates
				a[DIM];///,area;//face area-vector
			y=node->var+ixold;
			//Coordinates of face vertexes and Three edges
			for (int i=0; i<DIM; i++)
			{ //Face center
				double	xcenter=0.0;
				for (int j=0; j<Nfv; j++)
					xcenter+=vert[(iv+j+1)%Nv]->x[i];
				xf[i]=xcenter/(double)Nfv;
			}
			//Face area vector	
			VECP(a,e[iv1],e[iv2]);
			if (SCLP(a,e[iv])>0.0)
				for (double *p=a; p-a<DIM; p++) *p*=-.5;
			else
				for (double *p=a; p-a<DIM; p++) *p*= .5;
			dv=(SCLP(xf,a));
			if (node->state.boundary==0)
			{	*voln+=dv;//nodal scheme
				for (int i=0; i<DIM; i++)
					y[i]+=xf[i]*dv;
			}
			vol+=dv;
		}
		*volc=vol;
		cell=cell->next;
	}	while(cell!=dcell_root);
	}
	//Move nodes to new positions
	if (dnode_root!=NULL)
	{	static const double
			//step=0.2*length_scale,
			sq_four_over_sq_three=2.0/pow(3.0,0.25),
			surface_relaxation=relaxation_factor,
			volume_relaxation=relaxation_factor;
		DNode *node=dnode_root;
		do 
		{	double
				//force[DIM],totforce,
				//factor,//normalization factor
				//*f=node->f,
				*varn=node->var,
				voln=varn[ivoln],
				*x=node->x,
				*y=varn+ixold;
			//for (int i=1; i<DIM; i++) force[i]=0.0;
			if (node->state.boundary==0)
			{	double	volni;
				if (voln<=SMALL)
				{	fprintf
					(	stderr,
						"Domain::relax:ERROR: volume=%g too small at %g,%g,%g\n",
						voln,x[0],x[1],x[2]
					);
					exit(1);
				}
				volni=0.75/voln;
				for (int i=0; i<DIM; i++)
					x[i]=volume_relaxation*x[i]+(1.0-volume_relaxation)*y[i]*volni;
			}
			node=node->next;
		}	while (node!=dnode_root);
	}
}
void	Domain::loadGeom(char *filename)
{
	FILE	*inp;
	OPENREAD(filename,inp);
	if (option.verbose) printf("File %s opened for reading\n",filename);
	if(strstr(filename,".dat")!=NULL)
		loadGeomBIN(inp);
	else
	if(strstr(filename,".asc")!=NULL)
		loadGeomASC(inp);
	else
		fprintf(stderr,"Unrecognized file type\n");
	fclose(inp);
	if (option.verbose)
	{	printf("File %s closed\n",filename);FLUSH;
	}
}
void	Domain::loadGeomBIN(FILE *inp)
{
	int
		mvar,
		inode,nn,
		icell,nc;//number of nodes and cells
	DNode	*current_node,*last_node,**node;
	DCell	*current_cell,*last_cell,**cell;
	fread(&type,sizeof(type),1,inp);
	if (option.verbose) printf("\tFile type: BINARY\n\tDomain type = %d\n",type);
	if (type!=dynamic)
	{	fprintf
		(	stderr,
			"ERROR: Can only read dynamic domains, current type=%d: read operation skipped\n",
			type
		);
		exit(0);
	}
	//Read node data
	mvar=getVarBufSize(nodes);
	if (dnode_root!=NULL)
		deleteRing(dnode_root);
	//Read the number of nodes
	fread(&nn,sizeof(nn),1,inp);
	if (option.verbose)
		printf("\tNumber of nodes = %d\n",nn);
	//Allocate the node list
	node=new DNode*[nn];
	dnode_root=new DNode(mvar);
	last_node=dnode_root;
	last_node->index=0;
	last_node->next=last_node->prev=last_node;
	node[0]=dnode_root;
	for (int i=1; i<nn; i++)
	{	last_node->next=new DNode(mvar);
		last_node->next->prev=last_node;
		last_node=last_node->next;
		last_node->index=i;
		node[i]=last_node;
	}
	last_node->next=dnode_root;
	dnode_root->prev=last_node;
	for (inode=0; inode<nn && !feof(inp); inode++)
	{	//input nodes
		DNode	*current=node[inode];
		double	*x=current->x;
		///int	itype;
		//Read type
		fread(&current->state,sizeof(current->state),1,inp);
		///current->type=(ElementStatus)itype;
		//Read coordinates
		fread(x,sizeof(x[0]),DIM,inp);
	}
	if (option.verbose) printf("%d Node data read\n",inode);
	if (inode<nn)
	{	fprintf(stderr,"ERROR: FILE TOO SHORT\n");
		exit(1);
	}
	//Read connectivities
	//Read cell data
	mvar=getVarBufSize(cells);
	if (dcell_root!=NULL)
		deleteRing(dcell_root);
	fread(&nc,sizeof(nc),1,inp);
	if (option.verbose) printf("\tNumber of cells = %d\n",nc);
	//Allocate cell-list
	dcell_root=new DCell(mvar);
	last_cell=dcell_root;
	last_cell->index=0;
	last_cell->next=last_cell->prev=last_cell;
	for (int i=1; i<nc; i++)
	{	last_cell->next=new DCell(mvar);
		last_cell->next->prev=last_cell;
		last_cell=last_cell->next;
		last_cell->index=i;
	}
	last_cell->next=dcell_root;
	dcell_root->prev=last_cell;
	//Read cell-node connectivity
	current_cell=dcell_root;
	for (icell=0; icell<nc; icell++)
	{	ElementStatus	*facetype=current_cell->facetype;
///		fread(&current_cell->facetype,sizeof(current_cell->facetype[0]),Nf,inp);
		fread(facetype,sizeof(facetype[0]),Nf,inp);
		//Read cell vertex-nodes
		for (int i=0; i<Nv; i++)
		{	int ivert;
			fread(&ivert,sizeof(int),1,inp);
			if (ivert>=0)
				current_cell->vert[i]=node[ivert];
#ifdef DEBUG
			else
			{	fprintf(stderr,"loadGeom: cell-node connectivity error (ivert=%d)\n",ivert);
				exit(1);
			}
			if (current_cell->vert[i]==NULL)
			{	fprintf
				(	stderr,
					"loadGeom: icell=%d, node[%d(%d)]=%x\n",
					icell,ivert,i,node[ivert]
				);
				exit(1);
			}
#endif
		}
		current_cell=current_cell->next;
	}
	delete node;
	//Allocate temporal cells array
	cell=new DCell*[nc];
	last_cell=dcell_root;
	for (int i=0; i<nc; i++)
	{
		cell[i]=last_cell;
		last_cell=last_cell->next;
	}
	//Read cell-cell connectivity
	last_cell=dcell_root;
	for (int ic=0; ic<nc; ic++)
	{	for (int i=0; i<Nf; i++)
		{	int	icell;
			fread(&icell,sizeof(icell),1,inp);
			if (icell>=0)
				last_cell->neib[i]=cell[icell];
			else
				last_cell->neib[i]=NULL;
		}
		last_cell=last_cell->next;
	}
	delete cell;
#ifdef DEBUG
	{//Check connectivity
		int	i=0;
		DNode	*node=dnode_root;
		DCell	*cell=dcell_root;
	do//Zero out indexes 
	{
		node->index=0; 
		node=node->next;
	}	while(node!=dnode_root);
	do 
	{	DNode	**verts=cell->vert;
		//Face Areas and Face Fluxes
		for (int ivert=0; ivert<Nv; ivert++)
			verts[ivert]->index=1;
		cell=cell->next;
	}	while(cell!=dcell_root);
	node=dnode_root;
	i=0;
	do//Check connectivity
	{
		if (node->index==0) 
		{	fprintf(stderr,"loadGeom:: grid inconsistency: loose node %d\n",i++);
			exit(1);
		}
		node=node->next;
	}	while(node!=dnode_root);
	}
#endif
}
namespace CellIndex
{
	struct Cells
	{	int	icell;//Index into cell[icell] array
		Cells *next;
	};
	struct Nodes
	{	DNode *node;
		Cells	*cells;
	};
	void insertIntoList(int icell, Cells *&cell)
	{//Insert into the list in accentind order
		if(cell==NULL)
		{	cell=new Cells;
			cell->next=NULL;
			cell->icell=icell;
			return;
		}
		if(icell==cell->icell)return;
		Cells *next=cell->next;
		if(next!=NULL)
		{	if(icell>next->icell) 
				insertIntoList(icell,next);
			else
			{	//insert newcell between cell and next
				Cells *newcell=new Cells;
				newcell->next=next;
				cell->next=newcell;
				if(icell>cell->icell)
				{//place icell between cell->icell and next->icell
					newcell->icell=icell;
				}
				else
				{//place icell before cell->cell
					newcell->icell=cell->icell;
					cell->icell=icell;
				}
			}
		}
		else
		{//next==NULL
			Cells *newcell=new Cells;
			newcell->next=NULL;
			cell->next=newcell;
			if(icell>cell->icell)
				newcell->icell=icell;
			else
			{//place icell before cell->cell
				newcell->icell=cell->icell;
				cell->icell=icell;
			}
		}
	}
}
void	Domain::loadGeomASC(FILE *inp)
{	//Node numbering starts from 1
	char buf[MAXLINLEN];
	int
		mvar,
		inode,nn,
		icell,nc;//number of nodes and cells
	DNode	*current_node,*last_node,**node;
	DCell	*current_cell,*last_cell,**cell;
	type=dynamic; // Only dynamic domains read in ASCII
//	fread(&type,sizeof(type),1,inp);
	if (option.verbose) printf("\tFile type: ASCII\n\tDomain type = %d\n",type);
	if (type!=dynamic)
	{	fprintf
		(	stderr,
			"ERROR: Can only read dynamic domains, current type=%d: read operation skipped\n",
			type
		);
		exit(0);
	}
	//Read node data
	mvar=getVarBufSize(nodes);
	if (dnode_root!=NULL)
		deleteRing(dnode_root);
	//Read the number of nodes
///	fread(&nn,sizeof(nn),1,inp);
	do
	{	fgets(buf,MAXLINLEN,inp);
	}	while(!feof(inp)&&buf[0]=='#');
	sscanf(buf,"%d %d",&nn,&nc);
	if (option.verbose)
		printf("\tNumber of nodes = %d\n\tNumber of cells = %d\n",nn,nc);
	//Allocate the node list
	node=new DNode*[nn];
	dnode_root=new DNode(mvar);
	last_node=dnode_root;
	last_node->index=0;
	last_node->next=last_node->prev=last_node;
	node[0]=dnode_root;
	for (int i=1; i<nn; i++)
	{	last_node->next=new DNode(mvar);
		last_node->next->prev=last_node;
		last_node=last_node->next;
		last_node->index=i;
		node[i]=last_node;
	}
	last_node->next=dnode_root;
	dnode_root->prev=last_node;
	for (inode=0; inode<nn && !feof(inp); inode++)
	{	//input nodes
		int	jnode;
		double y[DIM];
		///int	itype;
		//Read type
		///fread(&current->state,sizeof(current->state),1,inp);
		fgets(buf,MAXLINLEN,inp);
		sscanf(buf,"%d %g %g %g",&jnode,y,y+1,y+2);
		if(--jnode>=nn)
		{	fprintf(stderr,"Node index %d exeeds number of nodes %d\n",jnode,nn);
			exit(1);
		}
		DNode	*current=node[jnode];
		double	*x=current->x;
		///current->type=(ElementStatus)itype;
		//Read coordinates
		///fread(x,sizeof(x[0]),DIM,inp);
		for(int i=0;i<DIM;i++)x[i]=y[i];
	}
	if (option.verbose) printf("%d Node data read\n",inode);
	if (inode<nn)
	{	fprintf(stderr,"ERROR: FILE TOO SHORT\n");
		exit(1);
	}
	//Read cell data
	mvar=getVarBufSize(cells);
	if (dcell_root!=NULL)
		deleteRing(dcell_root);
///	fread(&nc,sizeof(nc),1,inp);
///	if (option.verbose) printf("\tNumber of cells = %d\n",nc);
	//Allocate cell-list
	dcell_root=new DCell(mvar);
	last_cell=dcell_root;
	last_cell->index=0;
	last_cell->next=last_cell->prev=last_cell;
	for (int i=1; i<nc; i++)
	{	last_cell->next=new DCell(mvar);
		last_cell->next->prev=last_cell;
		last_cell=last_cell->next;
		last_cell->index=i;
	}
	last_cell->next=dcell_root;
	dcell_root->prev=last_cell;
	//Read cell-node connectivity
	current_cell=dcell_root;
	for (icell=0; icell<nc; icell++)
	{	int jcell=0;
		char *p=buf;
		ElementStatus	*facetype=current_cell->facetype;
		fgets(buf,MAXLINLEN,inp);
		sscanf(buf,"%d %d",&jcell,facetype);
		while(!isspace(*p))p++;while(isspace(*p))p++;
		while(!isspace(*p))p++;while(isspace(*p))p++;
		///fread(facetype,sizeof(facetype[0]),Nf,inp);
		//Read cell vertex-nodes
		for (int i=0; i<Nv; i++)
		{	int ivert;
			///fread(&ivert,sizeof(int),1,inp);
			while(!isspace(*p))p++;while(isspace(*p))p++;
			sscanf(p,"%d",&ivert);
			if (--ivert>=0)
				current_cell->vert[i]=node[ivert];
#ifdef DEBUG
			else
			{	fprintf(stderr,"loadGeom: cell-node connectivity error (ivert=%d)\n",ivert);
				exit(1);
			}
			if (current_cell->vert[i]==NULL)
			{	fprintf
				(	stderr,
					"loadGeom: icell=%d, node[%d(%d)]=%x\n",
					icell,ivert,i,node[ivert]
				);
				exit(1);
			}
#endif
		}
		current_cell=current_cell->next;
	}
	delete node;
	//Allocate temporal cells array
	cell=new DCell*[nc];
	last_cell=dcell_root;
	for (int i=0; i<nc; i++)
	{
		cell[i]=last_cell;
		last_cell=last_cell->next;
	}
	// Cell-cell connectivity
#ifdef READ_CELL_CELL_CON
	//Read cell-cell connectivity
	last_cell=dcell_root;
	for (int ic=0; ic<nc; ic++)
	{	for (int i=0; i<Nf; i++)
		{	int	icell;
			fread(&icell,sizeof(icell),1,inp);
			if (icell>=0)
				last_cell->neib[i]=cell[icell];
			else
				last_cell->neib[i]=NULL;
		}
		last_cell=last_cell->next;
	}
#else	//Construct cell-cell connectivity
#ifdef CELLINDEX
	{	using namespace CellIndex;
		//Construct the list of cells for each node
		Nodes *nodes = new Nodes[nn];
		{//Index all nodes in nodes array
			int inode=0;
			DNode *current_node=dnode_root;
			do
			{	current_node->index=inode;
				nodes[inode].node=current_node;
				nodes[inode].cells=NULL;
				current_node=current_node->next;
				inode++;
			}	while(current_node!=dnode_root);
		}
		//Run through all cells and fill in nodes[]
		for(int ic=0;ic<nc;ic++)
		{	current_cell=cell[ic];
			current_cell->index=ic;//assign temporal cell-index
			DNode **verts=current_cell->vert;	
			for(int iv=0;iv<Nv;iv++)
			{	DNode *vert=verts[iv];
				int inode=vert->index;
				// Insert ic into cells in accending order
				insertIntoList(ic,nodes[inode].cells);
			}
		}
		//Run through the cells and find neighbors
		for(int ic=0;ic<nc;ic++)
		{	current_cell=cell[ic];
			DNode **verts=current_cell->vert;
			DCell **neibs=(DCell**)current_cell->neib;
			ElementStatus *facetype=current_cell->facetype;
			for(int iface=0;iface<Nv;iface++)
			{//Loop through the faces
				//If all three vertexes belong to 
				// the same cell which is not the current_cell
				// then that cell is the neighbor
				int 
					ivert0= iface      ,
					ivert1=(iface+1)%Nv,
					ivert2=(iface+2)%Nv;
				DNode 
					*vert0=verts[ivert0],
					*vert1=verts[ivert1],
					*vert2=verts[ivert2];
				int 
					inode0=vert0->index,
					inode1=vert1->index,
					inode2=vert2->index,
					hitboundary=0;
				Cells 
					*cell0=nodes[inode0].cells,
					*cell1=nodes[inode1].cells,
					*cell2=nodes[inode2].cells;
				// sharing the current vertex, vert. 
				// Loop through this list and merge the list
				// into neibs list
				while(hitboundary==0&&(cell0->icell==ic||cell0->icell!=cell1->icell||cell1->icell!=cell2->icell))
				{
					while(cell1!=NULL&&(cell1->icell<cell0->icell||cell1->icell==ic))
					{	if(cell1->next==NULL) 
						{	hitboundary=1;
							break;
						}
						cell1=cell1->next;
					}
					if(hitboundary)break;
					while(cell2!=NULL&&(cell2->icell<cell1->icell||cell2->icell==ic)) 
					{	if(cell2->next==NULL) 
						{	hitboundary=1;
							break;
						}
						cell2=cell2->next;
					}
					if(hitboundary)break;
					while(cell0!=NULL&&(cell0->icell<cell2->icell||cell0->icell==ic))
					{	if(cell0->next==NULL) 
						{	hitboundary=1;
							break;
						}
						cell0=cell0->next;
					}
					if(hitboundary)break;
				}
				if(hitboundary)
				{	neibs[iface]=NULL;
					facetype[iface]=boundary;
				}
				else
				if(cell0->icell==cell1->icell&&cell1->icell==cell2->icell)
				{	neibs[iface]=cell[cell0->icell];
					facetype[iface]=internal;
				}
				else
				{	fprintf(stderr,"Can't locate neighbors\n");
					exit(1);
				}
			}
		}
		//Delete nodes
		for(int in=0;in<nn;in++)
			deleteList(nodes[in].cells);
		delete nodes;
	}
#else //NO CELLINDEX
	//Head-on double looping to find neighbors
	last_cell=dcell_root;
	for (int ic=0; ic<nc; ic++)
	{
		if(option.verbose)
		{	printf
			(	"\rConstructing cell-cell connectivity: cell %d of %d (%5.2f%%)",
				ic,nc,100.0*(float)(ic+1)/(float)nc
			);fflush(stdout);
		}
		for (int iv=0; iv<Nf; iv++)
		{	///int	icell;
			//fread(&icell,sizeof(icell),1,inp);
			//Find cell sharing the same face
			DNode *verts[Nv];//three vertexes of the triangular face
			for (int i=0;i<Nv1;i++)	verts[i]=last_cell->vert[(iv+i)%Nv];
			//Scan all other cells and for each face see if it has
			// the same vertexes
			current_cell=last_cell->next;
			int found=0;
			while(!found&&current_cell!=last_cell)
			{	int match=0;
				for(int jv=0;jv<Nv;jv++)
				{	DNode *neibvert=current_cell->vert[jv];
					for(int i=0;i<Nv1;i++)
					if(verts[i]==neibvert) 
					if(++match==Nv1)break;
				}
				if(match==Nv1)found=1;
				current_cell=current_cell->next;
			}
			if (found)
				last_cell->neib[iv]=current_cell;
			else
				last_cell->neib[iv]=NULL;
		}
		last_cell=last_cell->next;
	}
	if(option.verbose)printf("\n");
#endif  //END NO-CELLINDEX
	delete cell;
#endif
#ifdef DEBUG
	{//Check connectivity
	if(option.verbose){printf("Checking connectivity\n");fflush(stdout);}
		int	i=0;
		DNode	*node=dnode_root;
		DCell	*cell=dcell_root;
	do//Zero out indexes 
	{
		node->index=0; 
		node=node->next;
	}	while(node!=dnode_root);
	do 
	{	DNode	**verts=cell->vert;
		//Face Areas and Face Fluxes
		for (int ivert=0; ivert<Nv; ivert++)
			verts[ivert]->index=1;
		cell=cell->next;
	}	while(cell!=dcell_root);
	node=dnode_root;
	i=0;
	do//Check connectivity
	{
		if (node->index==0) 
		{	fprintf(stderr,"loadGeom:: grid inconsistency: loose node %d\n",i++);
			exit(1);
		}
		node=node->next;
	}	while(node!=dnode_root);
	}
#endif
}
void	Domain::loadGeomBZIP(FILE *inp)
{
	int
		mvar,
		inode,nn,
		icell,nc;//number of nodes and cells
	DNode	*current_node,*last_node,**node;
	DCell	*current_cell,*last_cell,**cell;
//	BZFILE* b;
//	int     nBuf;
//	char    buf[MAXLINLEN];
//	int     bzerror;
//	int     nWritten;
//	b = bzReadOpen ( &bzerror, inp, 0, 0, NULL, 0 );
//	if (bzerror != BZ_OK) {
//	   bzReadClose ( &bzerror, b );
//	   /* handle error */
//		fprintf(stderr,"ERROR READING BZIPEED FILE\n");fflush(stderr);
//	bzerror = BZ_OK;
//	while (bzerror == BZ_OK && /* arbitrary other conditions */) {
//	   nBuf = bzRead ( &bzerror, b, buf, MAXLINLEN );
//	   if (bzerror == BZ_OK) {
//	      /* do something with buf[0 .. nBuf-1] */
//			printf("BUF=i'");
//			for(int i=0;i<MAXLINLEN;i++)
//				printf("%c",buf[i]);
//			printf("'\n");fflush(stdout);
//	   }
//	}
//	if (bzerror != BZ_STREAM_END) {
//	   bzReadClose ( &bzerror, b );
//	   /* handle error */
//	} else {
//	   bzReadClose ( &bzerror );
//	}
}
