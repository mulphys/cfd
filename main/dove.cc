#include <iostream>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "vecalg.h"
#include "io.h"
#include "main.h"
#include "geom.h"
#include "var.h"
#include "tool.h"
#include "domain.h"

#ifdef WITH_MPI
#include "mpi.h"
using namespace MPI;
#endif

#ifdef WITH_DOVE
int	Domain::countBoundaryNodes()
{
	int	n=0;
	n_nodes=0;
	DNode	*node_root=dnode_root,*node=node_root;
	if (node_root==NULL)return 0;
	do 
	{
		n_nodes++; //count the number of all nodes; added by zhang
		if(node->state.boundary==1)n++; //count boundary nodes
		node=node->next;
	}	while(node!=node_root);
printf("countBoundaryNodes:dom=%d='%s':nbv=%d\n",this-domain_root,name,n);FLUSH;
	return n;
}
int	Domain::indexBoundaryNodes()
{	int	n=countBoundaryNodes();
	if(bnode!=NULL)delete bnode;
	bnode=new DNode*[n];
	if(dnode_root!=NULL)
	{	n=0;
		DNode	*node_root=dnode_root,*node=node_root;
		do 
		{
			if(node->state.boundary==1)
				bnode[n++]=node;
			node=node->next;
		}	while(node!=node_root);
	}
	nbv=n;//set the nuber of boundary vertexes
	return n;
}
void	Domain::storeBoundaryCoordinates()
{	int	n=indexBoundaryNodes();
	if(Xb!=NULL)delete Xb;
	Xb=new double[DIM*n];
	for(int inode=0;inode<n;inode++)
	{	DNode	*node=bnode[inode];
		double
			*x=node->x,
			*xb=Xb+inode*DIM;
		for(int i=0;i<DIM;i++)xb[i]=x[i];
	}
}

/*added on june 12, base on email "Gasflow model fix" */
void	Domain::setBoundaryFaceTypes()
{
	DCell	*cell=dcell_root;
	do
	{ 
		DNode **verts=cell->vert;
		int	ivert;
		for (ivert=0;ivert<Nv;ivert++) 
		{	// Looping over each vertex of the cell
			DNode *vert=verts[ivert];
			if(vert->type<=(int)maxElementStatus) break;
		}
		if (ivert<Nv)
			cell->facetype[ivert]=dove;
			//cell->facetype[ivert]=(int)maxElementStatus-(int)(this-domain_root);
		cell=cell->next;
	}	while(cell!=dcell_root);

}
/*********************/

void	Domain::addDoveI
(
	int nbn,  //number of neighbor boundary nodes 
	double *X,  //coordinate of nodes
	Domain *dom//parallel domain
)
{	//Compute overlaps between this domain
	// and the set of points given by 
	// coordinates X
	//Loop through all the cells
	int	n;
	DCell	*cell_root=dcell_root,
		*cell=cell_root,**C;
	C=new DCell*[nbn];
	for(int i=0;i<nbn;i++)
		C[i]=NULL;
	n=0;
	do//Mark all the contained nodes
	{	for(int inode=0;inode<nbn;inode++)
		{	double	*x=X+inode*DIM;
			if(ihostcell(x,cell)==-1)//node is inside
			{	C[inode]=cell;
				n++;
			}
		}
		cell=cell->next;
	}	while(cell!=cell_root);
	//Use directed search later ...
	if(n>0)//non-zero overlap
	{	Dove	*dove = new Dove;
		dove->dom=dom;
		//dove->n=n;
		dove->I=new int[n];
		dove->C=new DCell*[n];
		n=0;
		for(int inode=0;inode<nbn;inode++)
		if(C[inode]!=NULL) 
		{
			dove->I[n]=inode;
			dove->C[n]=C[inode];
			n++;
		}
		dove->n=n; //moved here from above by zhang
		dove_ring->insertBefore(dove);
	}
	delete C;
}

void	Domain::addDoveX
(
	int nbn,  //number of nodes 
	double *X,  //coordinate of nodes
	Domain *dom
)
{	//Compute overlaps between this domain
	// and the set of points given by 
	// coordinates X
	//Loop through all the cells
	int	n;
	DCell	*cell_root=dcell_root,
		*cell=cell_root,**C;
	C=new DCell*[nbn];
	for(int i=0;i<nbn;i++)
		C[i]=NULL;
	//Mark all the contained nodes
	n=0;
	for(int inode=0;inode<nbn;inode++)
	{
        double        *x=X+inode*DIM;
        DCell        *cell_root=dcell_root,
	                *cell=cell_root;
    	do
        {
	        if(ihostcell(x,cell)==-1)//node is inside
			{
				C[inode]=cell;
				n++;
				break;
			}
			cell=cell->next;
        }while(cell!=cell_root);
	}
 
	if(n>0)//non-zero overlap
	{	Dove	*dove = new Dove;
		dove->dom=dom;
		//dove->n=n;
		dove->I=new int[n];
		dove->X=new double[DIM*n];
		dove->C=new DCell*[n];
		n=0;
		for(int inode=0;inode<nbn;inode++)
		if(C[inode]!=NULL) 
		{	double
				*y=dove->X+DIM*n,
				*x=X+DIM*inode;
			dove->I[n]=inode;
			dove->C[n]=C[inode];
			for(int i=0;i<DIM;i++)y[i]=x[i];
			n++;
		}
		dove->n = n; //moved here from the above by zhang
		dove_ring->insertBefore(dove);
	}
	delete C;
}

int	receiveInt(int iproc)
{	int	n=0;

	return n;
}
void Domain::sendVecArray(int n, double *X, Domain *dom)
{

}
void Domain::receiveArraySize(Domain *dom)
{

}
int Domain::receiveVecArraySize(Domain *dom)
{
	int n;//number of coordinates to be received;
	
	return 0;
}
void receiveArray(int n, double *A, int iproc)
{

}
void sendIArray(int n, int *I, int iproc)
{

}
void receiveIArray(int n, int *I, int iproc)
{

}
void Domain::getArray(int idom, int n, double *X)
{// gets array addressed to domain idom from this domain

}
void Domain::receiveVecArray(int n, double *X, Domain *dom)
{	int	jproc=dom->iproc;
	if(iproc!=jproc)
	{	receiveArray(n*DIM,X,jproc);

	}
	else
	{	int idom=this-domain_root;
		dom->getArray(idom,n*DIM,X);
	}
}

void	Domain::ConnectThisProcDomains(int iproc)
{//Determine overlaps between this domain
	//and all other domains on the same processor
	//Should be called after the storeBoundaryCoordinates
	int	idom=this-domain_root;
printf("ConnectThisProcDomains:dom=%d\n",idom);///DDD
	for(int jdom=1;jdom<ndomains;jdom++)
	{
	Domain
			*rdom=domain_root+(idom+ndomains-jdom)%ndomains;
		if(iproc==rdom->iproc)
		{	int	n=rdom->nbv;
printf("\tnbv=%d\n",n);///DDD
			if(n>0)
				addDoveI(n,rdom->Xb,rdom);
		}
	}
	if(!dove_ring->isEmpty())
	{	//Send the boundary overlap-status information to all 
		// the domains
		Dove	*dove=dove_ring->currentElement(), *dove0=dove;

		do
		{	Domain	*dom=dove->dom;
			int	n=dove->n,*I=dove->I;
			if(dom->iproc==iproc)
			{	DNode	**neib_bnode=dom->bnode;
printf("Sending boundary overlap of %d from %d to %d I[0]=%d I[n-1]=%d\n",n,idom,dom-domain_root,I[0],I[n-1]);///DDD

				for(int inode=0;inode<n;inode++)
					neib_bnode[I[inode]]->type=(int)maxElementStatus-(int)(dom-domain_root);
			}
			dove=dove_ring->nextElement();
		}while(dove!=dove0);
	}
}

#ifdef WITH_MPI
//send boundary coordinates from domi to domj
void	Domain::sendBoundary(Domain *domi, Domain *domj)
{
  std::cout<<domi->iproc<<" send to "<<domj->iproc<<" dom_nbv: "<<domi->nbv;
std::cout<<" Xb: "<<domi->Xb[0]<<" "<<domi->Xb[(domi->nbv)*DIM-1]<<std::endl;

	//send the number of boundary vertices first
	COMM_WORLD.Send(&domi->nbv,1,MPI_INT,domj->iproc,0);
	//send the coordinates
	if(domi->nbv>0)
	{
		COMM_WORLD.Send(domi->Xb, (domi->nbv)*DIM, MPI_DOUBLE, domj->iproc,0);
	}
}
//receive boundary coordinates from domj to domi and compute overlap
void	Domain::receiveBoundary(Domain *domi, Domain *domj)
{
	int domj_nbv;
	COMM_WORLD.Recv(&domj_nbv,1, MPI_INT, domj->iproc,0);
	if(domj_nbv>0)
	{
		double *domj_Xb = new double[domj_nbv*DIM];
		COMM_WORLD.Recv(domj_Xb, (domj_nbv)*DIM, MPI_DOUBLE, domj->iproc,0);
		std::cout<<domi->iproc<<" recv from "<<domj->iproc<<" dom_nbv: "<<domj_nbv<<" Xb[0]"<<domj_Xb[0]<<" "<<domj_Xb[domj_nbv*DIM-1]<<std::endl;

		addDoveX(domj_nbv,domj_Xb,domj);
		delete domj_Xb;
	}

}

void	Domain::sendOverlap()
{//Sends the overlaps to all the connected domains
	//	int	idom=this-domain_root;
	if(!dove_ring->isEmpty())
	{	//Send the boundary overlap-status information to all 
		// the domains
		Dove	*dove=dove_ring->currentElement(), *dove0=dove;
		do
		{	Domain	*dom=dove->dom;
			int	n=dove->n,*I=dove->I;
			//if(dom->iproc!=iproc) sendIArray(n,I,dom->iproc);
			if(dom->iproc!=iproc)
			{
				COMM_WORLD.Send(&n, 1, MPI_INT, dom->iproc, 0);
				COMM_WORLD.Send(I, n, MPI_INT, dom->iproc, 0);
std::cout<<iproc<<" send overlap to "<<dom->iproc<<" : n="<<n<<" I[0]="<<I[0]<<" I[n-1]="<<I[n-1]<<std::endl;
			}
			dove=dove_ring->nextElement();
		}	while(dove!=dove0);
	}
}

void	Domain::receiveOverlap()
{//Receives overlap info. from all the connected domains
	int	idom=this-domain_root;
	if(!dove_ring->isEmpty())
	{	Dove	*dove=dove_ring->currentElement(), *dove0=dove;
		//Get/receive the boundary status from the parallel domains 
		// and set the parallel-boundary flags
		do
		{	Domain	*dom=dove->dom;
			//int	n=dove->n,*I=dove->I,*J;
			int	n,*I=dove->I,*J;
			if(dom->iproc!=iproc)
			{	//receive n from dom
				COMM_WORLD.Recv(&n, 1, MPI_INT, dom->iproc, 0);
				//allocate J
				J=new int[n];
				//receiveIArray(n,J,dom->iproc);
				COMM_WORLD.Recv(J, n, MPI_INT, dom->iproc, 0);
std::cout<<iproc<<" recv overlap from "<<dom->iproc<<" : n= "<<n<<" J[0]"<<J[0]<<" J[n-1]"<<J[n-1]<<std::endl;

				//Set boundary node-flags
				for(int i=0;i<n;i++)
					bnode[J[i]]->type=(int)maxElementStatus-(int)(dom-domain_root);
				delete J;
			}
			dove=dove_ring->nextElement();
		}	while(dove!=dove0);
	}
}
#endif

void Domain::putVar(int ivar, Dove *dove)
{//update overlap region on the same-processor domain
	
}
void Domain::sendVar(int ivar, Dove *dove)
{//update overlap region on the other-processor domain
	
}
int	Domain::receiveVar(int ivar, Dove *dove)
{	int nov,//number of overlapped nodes received
		*Ion;//indexes of the overlapped nodes 
	double	*V;//array of values received
	return 0;
}

void Domain::updateDove(int ivar)
{//Exchange variable ivar in the overlap-regions in all connected domains
	// on this processor
	extern void getscl(int jvarloc,DCell *cell,double *x, double *scl);
	int
		dim=variable[ivar].dimension,
		varloc=variable[ivar].loc;
	if(dove_ring->isEmpty()) return;
	Dove	*dove0=dove_ring->currentElement(),*dove=dove0;
	do
	{	Domain	*dom=dove->dom;
		int
			n=dove->n,
			*I=dove->I;
		DCell	**C=dove->C;
		if(model!=dom->model)ERROR("Can't overlap different model domains yet");
		if(iproc==dom->iproc)
		{//this proc
			for(int i=0;i<n;i++)
			{	DNode	*node=dom->bnode[I[i]];//node to be updated
				double	*x=node->x;
				DCell	*cell=C[i];
				//Do interpolation
				for(int jvar=0;jvar<dim;jvar++)
				{	int	jvarloc=varloc+jvar;
					double	scl;	
					getscl(jvarloc,cell,x,&scl);
					node->var[jvarloc]=scl;
				}
			}
		}
		dove=dove_ring->nextElement();
	}	while(dove!=dove0);
}

void Domain::updateDove()
{//Exchange all variables in the overlap-regions in all connected domains
	// on this processor
	extern void getvarbuf(int varbufsize,DCell *cell,double *x, double *var);
	int	varbufsize=getVarBufSize(nodes);
	if(dove_ring->isEmpty()) return;

	Dove	*dove=dove_ring->currentElement(),*dove0=dove;

	do
	{	Domain	*dom=dove->dom;
		int
			n=dove->n,
			*I=dove->I;
		DCell	**C=dove->C;
		if(iproc!=dom->iproc)return;
		if(model!=dom->model)ERROR("Can't overlap different model domains yet");
		for(int i=0;i<n;i++)
		{	DNode	*node=dom->bnode[I[i]];//node to be updated
			double	*x=node->x,
					*var=node->var;
			DCell	*cell=C[i];
			//Do interpolation
			getvarbuf(varbufsize,cell,x,var);
		}
		dove=dove_ring->nextElement();
	} while(dove!=dove0);
}

void Domain::sendDove(int ivar)
{//Send overlap variables to connected domains on other processors
	extern void getscl(int jvarloc,DCell *cell,double *x, double *scl);
	int
		dim=variable[ivar].dimension,
		varloc=variable[ivar].loc;
	if(dove_ring->isEmpty()) return;
	Dove	*dove = dove_ring->currentElement(), *dove0=dove;

	// The list below runs trough all elements of dove-list.
	// Each element of the dove-list holds info about one overlapping domain.
	do //For each overlapping domain do
	{	Domain	*dom=dove->dom;
		int
			n=dove->n,
			*I=dove->I;
		DCell	**C=dove->C;
		if(model!=dom->model)ERROR("Can't overlap different model domains yet");
#ifdef WITH_MPI
		if(iproc!=dom->iproc)
		{	// The overlapping domain is run by a different process: use
			//message-passing
			//Allocate send buffer
			double
				*X=dove->X,
				*V=new double[n*dim];// send-buffer to store variables to send
				// loop over all the nodes of the neighboring domain, which
				// are inside this domain. Their number is stored as dove->n:
			for(int i=0;i<n;i++)
			{	double
							*x=X+DIM*i,  // coordinates of the node i
							*v=V+dim*i; //pointer to the variable at node i
				DCell	*cell=C[i];
				// Since a variable can be of any size: 
				// the loop below is needed:
				for(int j=0;j<dim;j++)
				{// Loop over all the variables:
					int	jvarloc=varloc+j; // get the variable index inside the
					                         // variable buffer
					double	scl;	// store interpolated variable value here
					//Do interpolation:
					// this is needed since the node of the neighbor domain is located
					//somewhere
					// inside the cell of the current domain, and it's position does not have
					// to coinside with any of the vertexes where variables are stored:

//                               o A
//                              / \
//                             /   \
//                            /     \
//                           /   X   \
//                          /         \
//                         /           \
//                        o-------------o
//                        B             C

					// Here: A,B,C are the nodes of the cell of this domain dove->C[i]
					// which host the neighbor domain's cell: X
					// To get the value of the variable at X we need to interpolated it
					// from A,B,C to X. 
					// This is accomplished with the function below:
					getscl(jvarloc,cell,x,&scl);
					//Put the interpolated value into the send-buffer
					v[j]=scl;
				}
			}
			//send I[i],i=0:n as an array of integers, and
			// send V[i],i=0:n*dim as an array of doubles 
			//... 

			COMM_WORLD.Send(&n, 1 , MPI_INT, dom->iproc, 0);
			COMM_WORLD.Send(I, n, MPI_INT, dom->iproc, 0);
			COMM_WORLD.Send(V, n*dim, MPI_DOUBLE, dom->iproc, 0);
			/*
			std::cout<<"in senddove iproc="<<dom->iproc<<" n="<<n;
			std::cout<<" I[0]="<<I[0]<<" I[n-1]"<<I[n-1];
			std::cout <<" V[0]="<<V[0]<<" V[n*dim-1]"<<V[n*dim-1]<<std::endl;
			*/
			delete V;
		}
#endif
		dove=dove_ring->nextElement();
	}	while(dove!=dove0);
}

void Domain::receiveDove(int ivar)
{//Receive overlap variables from connected domains on other processors
	extern void getscl(int jvarloc,DCell *cell,double *x, double *scl);
	int
		dim=variable[ivar].dimension,
		varloc=variable[ivar].loc;
	if(dove_ring->isEmpty()) return;
	Dove	*dove = dove_ring->currentElement(), *dove0=dove;

	do
	{	Domain	*dom=dove->dom;
		int
			n=dove->n,
			dim=variable[ivar].dimension,
			varloc=variable[ivar].loc,
			*I=dove->I;
		DCell	**C=dove->C;
#ifdef WITH_MPI
		if(iproc!=dom->iproc)
		{
			int n2,*I2;
			double *V2;
			COMM_WORLD.Recv(&n2, 1, MPI_INT, dom->iproc, 0);
			I2 = new int[n2];
			COMM_WORLD.Recv(I2, n2, MPI_INT, dom->iproc, 0);
			V2 = new double[n2*dim];
			COMM_WORLD.Recv(V2, n2*dim, MPI_DOUBLE, dom->iproc, 0);
			/*
			std::cout<<"in recvdove iproc="<<dom->iproc<<" n2="<<n2;
			std::cout<<" I2[0]="<<I2[0]<<" I2[n2-1]"<<I2[n2-1];
			std::cout<<" V2[0]="<<V2[0]<<" V2[n2*dim-1]"<<V2[n2*dim-1]<<std::endl;
			*/
			for (int i=0;i<n2;i++)
			{	double
					*v=V2+i*dim; // address of the variable at node i in the 
					               // received V-array.
				DNode *node=bnode[I2[i]];				

				for(int j=0;j<dim;j++)
				{
					int	jvarloc=varloc+j; // get the variable index inside the
					node->var[jvarloc]=v[j];					
				}
			}
			delete I2;
			delete V2;
		}
#endif
		dove = dove_ring->nextElement();
	}while(dove!=dove0);
}

#endif

