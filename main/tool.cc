/**************************************
 * TOOL MANIPULATION ROUTINES         *
***************************************/

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
//#include "particles.h"
#include "tool.h"
#include "domain.h"
/**************************************
 * CREATE NODE AND CELL POINTER LISTS *
***************************************
 Global algorithm: loops through the whole grid
 Used when the previous location of the tool is not known
***************************************/
int	Tool::createToolNodeList
(
	Domain	*D
)
{	int	ntn;//Number of inside tool nodes	
	DNode	*dnode_root=D->dnode_root;
	double	toolpos[DIM],toolradius;
	toolradius=getradius();
	getpos(toolpos);
	if (first_tool_dnode!=NULL)
	for 
	(	DNodeList *node=first_tool_dnode->next;
		; node=node->next
	)
	{//relaxation radius should be greater or equal to the tool radius
		if (node->node!=NULL)
		{
			node->node->state.insidetool=0;
		}
		if (node==first_tool_dnode) break;
	}
	deleteList(first_tool_dnode);
	//Create insidetool, relax and boundary-node pointer lists
	ntn=0;
	if (dnode_root!=NULL)
	for 
	(	DNode	*node=dnode_root->next;
		; node=node->next
	)
	{
		DNodeList
			*tp;//tool-pointer
		double
			*x=node->x,
			d=0.0;//distance between the node and tool center
		for (int i=0; i<DIM; i++)
		{	double	r=x[i]-toolpos[i];
			d+=r*r;
		}
		d=sqrt(d);
		if (d<=toolradius)
		{//Create insidetool node pointer
			if (first_tool_dnode==NULL)
			{	first_tool_dnode=new DNodeList;
				tp=first_tool_dnode;
			}
			else
			{	tp->next=new DNodeList;
				tp=tp->next;
			}
			tp->node=node;
			tp->next=first_tool_dnode;//make it a loop
			node->state.insidetool=1;
			ntn++;
		}
		else
			node->state.insidetool=0;
		if (node==dnode_root) break;
	}
	return ntn;
}
int	Tool::createLists
(
	Domain	*D
)
{	int
		toolinside=0,
		ntbn;//number of boundary nodes to be returned
	DNode	*dnode_root=D->dnode_root;
	DCell	*dcell_root=D->dcell_root;
	double	toolpos[DIM],toolradius,relaxradius;
	toolradius=getinsideradius();
	relaxradius=getrelaxradius();
	getpos(toolpos);
	deleteLists();
	//Create insidetool, relax and boundary-node pointer lists
	ntbn=0;
	if (dnode_root!=NULL)
	{	DNode	*node=dnode_root;
	do 
	{
		DNodeList
			*bp,//boundary pointer
			*tp,//tool-pointer
			*rp,//relaxation pointer
			*brp;//boundary relaxation pointer
		double
			*x=node->x,
			d=0.0;//distance between the node and tool center
		for (int i=0; i<DIM; i++)
		{	double	r=x[i]-toolpos[i];
			d+=r*r;
		}
		d=sqrt(d);
		if (d<=relaxradius)
		{ //Create insidetool node pointer
			if (first_relax_dnode==NULL)
			{	first_relax_dnode=new DNodeList;
				rp=first_relax_dnode;
			}
			else
			{	rp->next=new DNodeList;
				rp=rp->next;
			}
			rp->node=node;
			rp->next=first_relax_dnode;//make it a loop
			if (d<=toolradius)
			{	if (node->state.boundary)
				{//Create boundary node pointer
					if (first_tool_boundary_dnode==NULL)
					{	first_tool_boundary_dnode=new DNodeList;
						bp=first_tool_boundary_dnode;
					}
					else
					{
						bp->next=new DNodeList;
						bp=bp->next;
					}
					bp->node=node;
					bp->next=first_tool_boundary_dnode;//make it a loop				
					ntbn++;
				}
				else
				{//Create insidetool node pointer
					if (first_tool_dnode==NULL)
					{	first_tool_dnode=new DNodeList;
						tp=first_tool_dnode;
					}
					else
					{	tp->next=new DNodeList;
						tp=tp->next;
					}
					tp->node=node;
					tp->next=first_tool_dnode;//make it a loop
				}
				node->state.insidetool=1;
			}
			else
				node->state.insidetool=0;
			node->state.relax=1;
		}
		else//Outside of relaxation radius
			node->state.relax=0;
		node=node->next;
	}	while (node!=dnode_root);
	}
	//Create tool-cell pointer lists
//	deleteList(first_tool_boundary_dcell);
//	deleteList(first_relax_dcell);
	root_cell=NULL;
	if (dcell_root!=NULL)
	{	DCell	*c=dcell_root;
	do
	{	int	nrv,nbv;//number of tool and relaxation vertexes
		DCellList	*rc,*tc,*rbc,*tbc; //relaxation and tool cell-pointers
		DNode	**vert=c->vert;
		//Check if this is the tool_root cell;
		if ((root_cell=D->hostcell(toolpos,c))==c)
			toolinside=1;
		//Check the status of the vertexes
		nrv=nbv=0;
		for (int i=0; i<Nv; i++)
		{
			if (vert[i]->state.relax)
			{	nrv++;
				///if (vert[i]->state.insidetool&&vert[i]->type.boundary)
			//	if (vert[i]->state.insidetool)
			//		ntv++; 
			//	//else
			//	if (vert[i]->state.boundary)
			//		nbv++;
				if (vert[i]->state.insidetool)
				{ if (vert[i]->state.boundary)
					{	nbv++;break;
					}
				}
			}
		}
		if (nrv>0)//relaxation-cell
		{	if (first_relax_dcell==NULL)
			{	first_relax_dcell=new DCellList;
				rc=first_relax_dcell;
			}
			else
			{
				rc->next=new DCellList;
				rc=rc->next;
			}
			rc->cell=c;
			rc->next=first_relax_dcell;//make it a loop	
			//if (ntv>0&&nbv>0)
			if (nbv>0)
			{
				if (first_tool_boundary_dcell==NULL)
				{	first_tool_boundary_dcell=new DCellList;
					tbc=first_tool_boundary_dcell;
				}
				else
				{	tbc->next=new DCellList;
					tbc=tbc->next;
				}
				tbc->cell=c;
				tbc->next=first_tool_boundary_dcell;//make it a loop	
			}
		}
		c=c->next;
	}	while (c!=dcell_root);
	}
	if (option.verbose)
		if(toolinside==0&&state==add_cells)
			fprintf(stderr,"WARNING: Tool is out of domain bounds\n");
	return ntbn;
}
/***************************************
 *  UPDATE NODE AND CELL POINTER LISTS *
 ***************************************
 Local algorithm: 
 - Serch for the tool's root point
   using the previous root point position
	 and cell connectivity structure
 - Create tool nodes/cells pointer lists
   by means of recursive sweeping through the cells
	 starting from the critical point and using
	 cell-connectivity information
****************************************/
int	Tool::updateLists
(
	Domain	*D
)
{	int ntbn;//number of tool boundary nodes
	double	toolpos[DIM];
	if (root_cell==NULL) return 0;
	deleteLists();
	getpos(toolpos);//get new tool position
	//Find new root_cell, i.e. the cell containing toolpos point
	if ((root_cell=D->findcell(toolpos,root_cell))==NULL) 
		return createLists(D); 
	//Create node/cell pointer lists
	nbn=0;
	raid(root_cell);
	ntbn=nbn;
	return ntbn;
}
void	Tool::cleanSurface()
{
	if (first_relax_dnode!=NULL)
	for  //Reset toolsurface flags
	(	DNodeList *rnode=first_relax_dnode->next;
		; rnode=rnode->next
	)
	{//relaxation radius should be greater or equal to the tool radius
		DNode	*node=rnode->node;
		if (node==NULL) continue;
		node->state.toolsurface=0;
		if (rnode==first_relax_dnode) break;
	}
}
void	Tool::deleteLists()
{				
	//Clear flags
	if (first_relax_dnode!=NULL)
	for 
	(	DNodeList *rnode=first_relax_dnode->next;
		; rnode=rnode->next
	)
	{//relaxation radius should be greater or equal to the tool radius
		DNode	*node=rnode->node;
		if (node==NULL) goto endnodeloop;
		node->state.relax=0;
		node->state.insidetool=0;
		//node->state.toolsurface=0;
		endnodeloop: if (rnode==first_relax_dnode) break;
	}
	if (first_relax_dcell!=NULL)
	for 
	(	DCellList *rcell=first_relax_dcell->next;
		; rcell=rcell->next
	)
	{	DCell	*cell=rcell->cell;
		if (cell==NULL) goto endcelloop;
		cell->state.relax=0;
		endcelloop: if (rcell==first_relax_dcell) break;
	}
	deleteList(first_tool_dnode);
	deleteList(first_tool_boundary_dnode);
	deleteList(first_relax_dnode);
	deleteList(first_tool_boundary_dcell);
	deleteList(first_relax_dcell);
}
/************************************
 * RAID CELLS AROUND THE ROOT CELL  *
 * TO CREATE THE TOOL POINTER LISTS *
 ** *********************************/
void	Tool::raid(DCell *cell)
{	int
		nrv,//number of relaxation vertexes and nodes
		nbv;//number of boundary vertexes
	double
		relaxradius,
		toolradius;
	DNode	**vert;
	DCell	**neib;
	if (cell==NULL) return;
	if (cell->state.relax==1) return;
	vert=cell->vert;
	neib=(DCell**)cell->neib;
	relaxradius=getrelaxradius();
	toolradius=getinsideradius();
	//getpos(toolpos);
	//Set flags
	//Assign node pointers
	nrv=nbv=0;
	for (int iv=0; iv<Nv; iv++)
	{//Find distance from the tool-root
		DNode	*node=vert[iv];
		double	d,
			*x=vert[iv]->x;
			d=0.0;//distance between the node and tool center
		for (int i=0; i<DIM; i++)
		{	double	r=x[i]-rootpos[i];
			d+=r*r;
		}
		d=sqrt(d);
		if (d<=relaxradius)
		{ //Create insidetool node pointer
			nrv++;
			if (node->state.relax==0)
			{	DNodeList	*rp;
				if (first_relax_dnode==NULL)
				{	first_relax_dnode=new DNodeList;
					rp=first_relax_dnode;
				}
				else
				{
					rp=first_relax_dnode->next;
					first_relax_dnode->next=new DNodeList;
					first_relax_dnode=first_relax_dnode->next;
				}
				first_relax_dnode->node=node;
				first_relax_dnode->next=rp;//make it a loop
				node->state.relax=1;
				//nrn++;
			}
			if (d<=toolradius)
			{
				if (node->state.insidetool==0)
				{ node->state.insidetool=1;
					if (node->state.boundary)
					{//Create boundary node pointer
						DNodeList	*bp;
						if (first_tool_boundary_dnode==NULL)
						{	first_tool_boundary_dnode=new DNodeList;
							bp=first_tool_boundary_dnode;
						}
						else
						{	bp=first_tool_boundary_dnode->next;
							first_tool_boundary_dnode->next=new DNodeList;
							first_tool_boundary_dnode=first_tool_boundary_dnode->next;
						}
						first_tool_boundary_dnode->node=node;
						first_tool_boundary_dnode->next=bp;//make it a loop				
						nbn++;
					}
					else
					{//Create insidetool node pointer
						DNodeList	*tp;
						if (first_tool_dnode==NULL)
						{	first_tool_dnode=new DNodeList;
							tp=first_tool_dnode;
						}
						else
						{	tp=first_tool_dnode->next;
							first_tool_dnode->next=new DNodeList;
							first_tool_dnode=first_tool_dnode->next;
						}
						first_tool_dnode->node=node;
						first_tool_dnode->next=tp;//make it a loop
					}
				}
				if (node->state.boundary) nbv++;
			}
			else
				node->state.insidetool=0;
		}
		else//Outside of relaxation radius
			node->state.relax=0;
	}
	//Assign cell pointers
	if (nrv>0)//relaxation-cell
	{	DCellList	*rc;
		if (first_relax_dcell==NULL)
		{	first_relax_dcell=new DCellList;
			rc=first_relax_dcell;
		}
		else
		{
			rc=first_relax_dcell->next;
			first_relax_dcell->next=new DCellList;
			first_relax_dcell=first_relax_dcell->next;
		}
		first_relax_dcell->cell=cell;
		first_relax_dcell->next=rc;//make it a loop	
		cell->state.relax=1;	
		if (nbv>0)
		{	DCellList	*tbc;
			if (first_tool_boundary_dcell==NULL)
			{	first_tool_boundary_dcell=new DCellList;
				tbc=first_tool_boundary_dcell;
			}
			else
			{	tbc=first_tool_boundary_dcell->next;
				first_tool_boundary_dcell->next=new DCellList;
				first_tool_boundary_dcell=first_tool_boundary_dcell->next;
			}
			first_tool_boundary_dcell->cell=cell;
			first_tool_boundary_dcell->next=tbc;//make it a loop	
		}
	}
	else	return;
	//Raid neighbor cells
	for (int i=0; i<Nv; i++)
		raid(neib[i]);
}
void	Tool::setboundary()
{
	if (first_tool_boundary_dcell!=NULL)
	{	DCellList	*boundary_cell=first_tool_boundary_dcell;
		do
		{	DCell	
				*cell=boundary_cell->cell,
				**neib=(DCell**)cell->neib;
			ElementStatus	*ftype=cell->facetype;
			for (int i=0; i<Nf; i++)
				if (neib[i]==NULL)
					ftype[i]=etype;
				else
					ftype[i]=internal;
			boundary_cell=boundary_cell->next;
		}	while(boundary_cell!=first_tool_boundary_dcell);
	}
}
/***********************************
 *   GROW CELLS: CORAL MODEL       *
 ***********************************
 Cell refinment is done by 
 looping through the cells and
 splitting the longest edge exceeding
 the refinement scale.
 All the cells around this edge are
 then split successively
 ***********************************/
int	Tool::growCoral
(
	Domain	*domain	
)
{	int refine=0,
		tool_boundary_nodes,//number of tool-boundary nodes
		mnvar=domain->getVarBufSize(nodes),
		mcvar=domain->getVarBufSize(cells);
	double
		length_scale=getlengthscale(),//max distance between the nodes
	                  // above which the cell has to be refined
		tool_pos[DIM];
	getpos(tool_pos);
	do
	{	//REFINE-ADVANCE CYCLE
		refine=0;
		if (first_tool_boundary_dcell!=NULL)
		for 
		(	DCellList	*current_boundary_cell=first_tool_boundary_dcell->next;
			; current_boundary_cell=current_boundary_cell->next
		)
		{	//Loop through all vertexes and compute max edge length
			int
				ivert,jvert; //two vertexes which edge has to be refined
			DCell	*cell=current_boundary_cell->cell,
				*newclone,*oldclone,*first_clone;
			DNode	**vert=cell->vert;
			double	max_edge_length=0,*x,*y;
			//Find the longest edge
			x=vert[0]->x,*y;
			for (int iv=1; iv<=Nv; iv++)
			{	int	jv=iv%Nv;
				double	edge_length=0.0;
				y=vert[jv]->x;
				for (int k=0; k<DIM; k++)
				{	double	d=y[k]-x[k];
					edge_length+=d*d;
				}
				edge_length=sqrt(edge_length);
				if (edge_length>max_edge_length)
				{	max_edge_length=edge_length;
					ivert=iv-1;
					jvert=jv;
				}
				x=y;
			}
			//Examine two remaining vertexes
			for (int iv=0; iv<2; iv++)
			{	double	edge_length=0.0;
				x=vert[iv]->x;
				y=vert[iv+2]->x;
				for (int k=0; k<DIM; k++)
				{	double	d=y[k]-x[k];
					edge_length+=d*d;
				}
				edge_length=sqrt(edge_length);
				if (edge_length>max_edge_length)
				{	max_edge_length=edge_length;
					ivert=iv;
					jvert=iv+2;
				}
			}
			if (max_edge_length>length_scale)
			{//Refine the current cell along the 
				// edges between ivert and jvert
				int
					icell0,icell1,
					inext_neib_cell_old,
					iprev_neib_cell_old;
				DNode	*p,*pp,*newvert,
					*verti=vert[ivert],
				 	*vertj=vert[jvert],
					*common_vert,
					*vert0,*vert1;
				DCell	*cell0,*cell1,
					*next_neib_cell,
					*next_neib_cell_old,
					**neib=(DCell**)cell->neib;
				//Find two other vertexes beside ivert and jvert
				//and use them to loop around the edge
				vert0=vert1=NULL;
				for (int iv=0; iv<Nv; iv++)
				{	if (iv!=ivert&&iv!=jvert)
					{ if (vert0==NULL) 
						{	vert0=vert[iv];
							cell0=neib[iv];
							icell0=iv;
						}
						else
						{	vert1=vert[iv];
							cell1=neib[iv];
							icell1=iv;
						}
					}
				}
				//Create a new vertex
				newvert=new DNode(mnvar);
				//Initialize flags
				newvert->state.relax=0;
				newvert->state.insidetool=0;
				newvert->state.toolsurface=0;
				//Assign coordinates
				for (int k=0; k<DIM; k++)
					newvert->x[k]=0.5*(verti->x[k]+vertj->x[k]);
				//Insert the vertex into the DNode chain
				domain->insert(verti,newvert);
				//Split the current cell's (verti,vertj) edge
				first_clone=NULL;
				domain->addDCell
				(	mcvar,ivert,jvert,
					newvert,
					cell,
					-1,-1,
					first_clone,//the split neighbor cell's twin
					newclone //returned new cloned cell
				);
				first_clone=newclone;
				oldclone=first_clone;
				//First follow the direction of cell0 around the edge
				common_vert=vert1;
				next_neib_cell=cell0;
				inext_neib_cell_old=icell0;
				iprev_neib_cell_old=icell1;
				next_neib_cell_old=next_neib_cell;
				while(next_neib_cell!=NULL&&next_neib_cell!=cell)
				{ //Determine local indexes iivert,jjvert of vertexes 
					//verti and vertj inside neib_cell
					int	iivert=-1,jjvert=-1,
						inext_neib_cell,//local index to the next neighbor cell
						iprev_neib_cell;//local index of the vertex common between
						           //this cell and next neighbor cell 
											 //(which is not verti and vertj)
											 //This is also the index to the previous neighbor cell
											 //which is only true for a tetrahedron
					DNode	**next_vert=next_neib_cell->vert;
					DCell	**next_neib=(DCell**)next_neib_cell->neib;
					for (int iv=0; iv<Nv; iv++)
					{	DNode	*v=next_vert[iv];
						if (v==verti) iivert=iv;
						else
						if (v==vertj) jjvert=iv;
						else
						if (v==common_vert) inext_neib_cell=iv;
						else
							iprev_neib_cell=iv;
					}
					domain->addDCell
					(	mcvar,iivert,jjvert,newvert,
						next_neib_cell,
						iprev_neib_cell,
						inext_neib_cell_old,
						oldclone,newclone
					);
					//Goto the next neighbor cell
					common_vert=next_vert[iprev_neib_cell];
					next_neib_cell=next_neib[inext_neib_cell];
					inext_neib_cell_old=inext_neib_cell;
					oldclone=newclone;
					//ivert_old=iivert; jvert_old=jjvert;
				}
				if (next_neib_cell==NULL)//if the wall was hit then
				{//Go around in another direction 
								//do as before
					oldclone=first_clone;
					common_vert=vert0;
					next_neib_cell=cell1;
					inext_neib_cell_old=icell1;
					next_neib_cell_old=next_neib_cell;
					while(next_neib_cell!=NULL)
					{//	addDCell(neib_cell->cell,newvert);
						//Goto the next neighbor cell
						int	iivert=-1,jjvert=-1,
							inext_neib_cell,//local index to the next neighbor cell
							iprev_neib_cell;//local index of the vertex common between
						           //this cell and next neighbor cell 
											 //(which is not verti and vertj)
											 //This is also the index to the previous neighbor cell
											 //which is only true for a tetrahedron
						DNode	**next_vert=next_neib_cell->vert;
						DCell	**next_neib=(DCell**)next_neib_cell->neib;
						for (int iv=0; iv<Nv; iv++)
						{	DNode	*v=next_vert[iv];
							if (v==verti) iivert=iv;
							else
							if (v==vertj) jjvert=iv;
							else
							if (v==common_vert) inext_neib_cell=iv;
							else
								iprev_neib_cell=iv;
						}
						domain->addDCell
						(	mcvar,iivert,jjvert,newvert,
							next_neib_cell,
							iprev_neib_cell,
							inext_neib_cell_old,
							oldclone,newclone
						);
						//Goto the next neighbor cell
						common_vert=next_vert[iprev_neib_cell];
						next_neib_cell=next_neib[inext_neib_cell];
						inext_neib_cell_old=inext_neib_cell;
						oldclone=newclone;
					}
					//Set type
					newvert->state.boundary=1;
				}
				else// the wall was not hit
				{//Sew current_boundary_cell to the other end of the loop
					newclone->neib[inext_neib_cell_old]=first_clone;
					first_clone->neib[iprev_neib_cell_old]=newclone;
					//newvert->type.insidetool=1;
					//newvert->type.relax=1;
					newvert->state.boundary=0;
				}
				refine=1;
			}
			if (current_boundary_cell==first_tool_boundary_dcell) break;
		}//end refinement loop
		//Recount all boundary and near-boundary nodes
		// and reset the refinement flags
#ifdef GLOBAL_LIST_UPDATE
		tool_boundary_nodes=createLists(domain);
#else
		tool_boundary_nodes=updateLists(domain);
#endif
		// Relax all boundary and near-boundary nodes
		relax(domain,0.5);
		endloop: break;//DEBUG
	}	while (refine);
	advanceBoundaryNodes();
	return tool_boundary_nodes;
}
/***********************************
 *   GROW CELLS: CRYSTAL MODEL     *
 ***********************************
Select the biggest boundary face
Select where the new node should be placed so as to
        create a new cell (tetrahedron)
        with this face at the base.
See if this new node is close to already existing nodes
        and if so consider using the existing node
        instead of the new one.
Create the new cell and set the connectivities to
        existing cells.
 ***********************************/
int	Tool::growCrystal
(
	Domain	*domain	
)
{	int refine=0,
		tool_boundary_nodes,//number of tool-boundary nodes
		mnvar=domain->getVarBufSize(nodes),
		mcvar=domain->getVarBufSize(cells);
	double
		length_scale=getlengthscale(),//max distance between the nodes
	                  // above which the cell has to be refined
		tool_pos[DIM];
	getpos(tool_pos);
	do
	{	//REFINE-ADVANCE CYCLE
		refine=0;
		if (first_tool_boundary_dcell!=NULL)
		for 
		(	DCellList	*current_boundary_cell=first_tool_boundary_dcell->next;
			; current_boundary_cell=current_boundary_cell->next
		)
		{	//Loop through all vertexes and compute max edge length
			int
				ivert,jvert; //two vertexes which edge has to be refined
			DCell	*cell=current_boundary_cell->cell,
				*newclone,*oldclone,*first_clone;
			DNode	**vert=cell->vert;
			double	max_edge_length=0,*x,*y;
			//Find the longest edge
			x=vert[0]->x,*y;
			for (int iv=1; iv<=Nv; iv++)
			{	int	jv=iv%Nv;
				double	edge_length=0.0;
				y=vert[jv]->x;
				for (int k=0; k<DIM; k++)
				{	double	d=y[k]-x[k];
					edge_length+=d*d;
				}
				edge_length=sqrt(edge_length);
				if (edge_length>max_edge_length)
				{	max_edge_length=edge_length;
					ivert=iv-1;
					jvert=jv;
				}
				x=y;
			}
			//Examine two remaining vertexes
			for (int iv=0; iv<2; iv++)
			{	double	edge_length=0.0;
				x=vert[iv]->x;
				y=vert[iv+2]->x;
				for (int k=0; k<DIM; k++)
				{	double	d=y[k]-x[k];
					edge_length+=d*d;
				}
				edge_length=sqrt(edge_length);
				if (edge_length>max_edge_length)
				{	max_edge_length=edge_length;
					ivert=iv;
					jvert=iv+2;
				}
			}
			if (max_edge_length>length_scale)
			{//Refine the current cell along the 
				// edges between ivert and jvert
				int
					icell0,icell1,
					inext_neib_cell_old,
					iprev_neib_cell_old;
				DNode	*p,*pp,*newvert,
					*verti=vert[ivert],
				 	*vertj=vert[jvert],
					*common_vert,
					*vert0,*vert1;
				DCell	*cell0,*cell1,
					*next_neib_cell,
					*next_neib_cell_old,
					**neib=(DCell**)cell->neib;
				//Find two other vertexes beside ivert and jvert
				//and use them to loop around the edge
				vert0=vert1=NULL;
				for (int iv=0; iv<Nv; iv++)
				{	if (iv!=ivert&&iv!=jvert)
					{ if (vert0==NULL) 
						{	vert0=vert[iv];
							cell0=neib[iv];
							icell0=iv;
						}
						else
						{	vert1=vert[iv];
							cell1=neib[iv];
							icell1=iv;
						}
					}
				}
				//Create a new vertex
				newvert=new DNode(mnvar);
				//Initialize flags
				newvert->state.relax=0;
				newvert->state.insidetool=0;
				newvert->state.toolsurface=0;
				//Assign coordinates
				for (int k=0; k<DIM; k++)
					newvert->x[k]=0.5*(verti->x[k]+vertj->x[k]);
				//Insert the vertex into the DNode chain
				domain->insert(verti,newvert);
				//Split the current cell's (verti,vertj) edge
				first_clone=NULL;
				domain->addDCell
				(	mcvar,ivert,jvert,
					newvert,
					cell,
					-1,-1,
					first_clone,//the split neighbor cell's twin
					newclone //returned new cloned cell
				);
				first_clone=newclone;
				oldclone=first_clone;
				//First follow the direction of cell0 around the edge
				common_vert=vert1;
				next_neib_cell=cell0;
				inext_neib_cell_old=icell0;
				iprev_neib_cell_old=icell1;
				next_neib_cell_old=next_neib_cell;
				while(next_neib_cell!=NULL&&next_neib_cell!=cell)
				{ //Determine local indexes iivert,jjvert of vertexes 
					//verti and vertj inside neib_cell
					int	iivert=-1,jjvert=-1,
						inext_neib_cell,//local index to the next neighbor cell
						iprev_neib_cell;//local index of the vertex common between
						           //this cell and next neighbor cell 
											 //(which is not verti and vertj)
											 //This is also the index to the previous neighbor cell
											 //which is only true for a tetrahedron
					DNode	**next_vert=next_neib_cell->vert;
					DCell	**next_neib=(DCell**)next_neib_cell->neib;
					for (int iv=0; iv<Nv; iv++)
					{	DNode	*v=next_vert[iv];
						if (v==verti) iivert=iv;
						else
						if (v==vertj) jjvert=iv;
						else
						if (v==common_vert) inext_neib_cell=iv;
						else
							iprev_neib_cell=iv;
					}
					domain->addDCell
					(	mcvar,iivert,jjvert,newvert,
						next_neib_cell,
						iprev_neib_cell,
						inext_neib_cell_old,
						oldclone,newclone
					);
					//Goto the next neighbor cell
					common_vert=next_vert[iprev_neib_cell];
					next_neib_cell=next_neib[inext_neib_cell];
					inext_neib_cell_old=inext_neib_cell;
					oldclone=newclone;
					//ivert_old=iivert; jvert_old=jjvert;
				}
				if (next_neib_cell==NULL)//if the wall was hit then
				{//Go around in another direction 
								//do as before
					oldclone=first_clone;
					common_vert=vert0;
					next_neib_cell=cell1;
					inext_neib_cell_old=icell1;
					next_neib_cell_old=next_neib_cell;
					while(next_neib_cell!=NULL)
					{//	addDCell(neib_cell->cell,newvert);
						//Goto the next neighbor cell
						int	iivert=-1,jjvert=-1,
							inext_neib_cell,//local index to the next neighbor cell
							iprev_neib_cell;//local index of the vertex common between
						           //this cell and next neighbor cell 
											 //(which is not verti and vertj)
											 //This is also the index to the previous neighbor cell
											 //which is only true for a tetrahedron
						DNode	**next_vert=next_neib_cell->vert;
						DCell	**next_neib=(DCell**)next_neib_cell->neib;
						for (int iv=0; iv<Nv; iv++)
						{	DNode	*v=next_vert[iv];
							if (v==verti) iivert=iv;
							else
							if (v==vertj) jjvert=iv;
							else
							if (v==common_vert) inext_neib_cell=iv;
							else
								iprev_neib_cell=iv;
						}
						domain->addDCell
						(	mcvar,iivert,jjvert,newvert,
							next_neib_cell,
							iprev_neib_cell,
							inext_neib_cell_old,
							oldclone,newclone
						);
						//Goto the next neighbor cell
						common_vert=next_vert[iprev_neib_cell];
						next_neib_cell=next_neib[inext_neib_cell];
						inext_neib_cell_old=inext_neib_cell;
						oldclone=newclone;
					}
					//Set type
					newvert->state.boundary=1;
				}
				else// the wall was not hit
				{//Sew current_boundary_cell to the other end of the loop
					newclone->neib[inext_neib_cell_old]=first_clone;
					first_clone->neib[iprev_neib_cell_old]=newclone;
					//newvert->type.insidetool=1;
					//newvert->type.relax=1;
					newvert->state.boundary=0;
				}
				refine=1;
			}
			if (current_boundary_cell==first_tool_boundary_dcell) break;
		}//end refinement loop
		//Recount all boundary and near-boundary nodes
		// and reset the refinement flags
#ifdef GLOBAL_LIST_UPDATE
		tool_boundary_nodes=createLists(domain);
#else
		tool_boundary_nodes=updateLists(domain);
#endif
		// Relax all boundary and near-boundary nodes
		relax(domain,0.5);
		endloop: break;//DEBUG
	}	while (refine);
	advanceBoundaryNodes();
	return tool_boundary_nodes;
}
void	Tool::push()
{
	double
		*x0=current_framenode->x,
		*x1=current_framenode->next->x,
		df,//distance between two framenodes
		dfx[DIM];//its unit direction vector
	DNodeList	*bnode=first_tool_boundary_dnode;
	if(bnode==NULL)return;
	df=0.0;
	for(int i=0; i<DIM; i++)
	{	double	r=x1[i]-x0[i];
		dfx[i]=r;
		df+=r*r;
	}
	df=sqrt(df);
	for(int i=0; i<DIM; i++)dfx[i]/=df;
	do
	{	DNode	*node=bnode->node;
		double	*x=node->x,
			dx[DIM],
			d=0.0;
		for(int i=0; i<DIM; i++)
			dx[i]=x0[i]-x[i];
		d=SCLP(dfx,dx);
		if(d>0.0)
		for(int i=0; i<DIM; i++)
			x[i]+=d*dfx[i];
		bnode=bnode->next;
	}	while(bnode!=first_tool_boundary_dnode);
}
void	Tool::relax
(
 	Domain *domain,
	double	relaxation_factor
)
{	static	const	double
		onethird=1.0/3.0;
	int
		ivol=domain->variable[ivolume].loc,
		ixold_loc=domain->variable[ixold].loc,
		nodes_deleted,///node_deleted,
		cells_deleted;
	double
		toolpos[DIM],
		tool_radius=getradius(),
		length_scale=0.5*getlengthscale(),
		length_scalei=30.0/length_scale,
		length_scale3=pow(length_scale,3.0);
	getpos(toolpos);
	//RELAXATION
	//Relax volume nodes
	if (first_relax_dnode!=NULL&&first_relax_dcell!=NULL)
	do
	{	///node_deleted=0;
		DCellList	*c=first_relax_dcell;
		nodes_deleted=cells_deleted=0;
	//Initialize 
	for
	(	DNodeList	*node=first_relax_dnode->next;
		; node=node->next
	)
	{	DNode	*nd=node->node;
		double
			*xold=nd->var+ixold_loc;
		//Initialize to zero
		nd->var[ivol]=0.0;
		for (int i=0; i<DIM; i++) xold[i]=0.0;
		if (node==first_relax_dnode) break;
	}
	do
	{	int	kill,surekill;//delete cell flag
		double	vol,*x,*y,
			areas[Nv],
			e[Nv][DIM];//tetrahedral edge-vectors					
		DCell	*cell=c->cell;
		DNode	**vert=cell->vert;
		DCell	**neib=(DCell**)cell->neib;
		//Edges of the tetrahedron:
		x=vert[0]->x;
		for (int iv=0; iv<Nv; iv++)
		{ double
		    *y=vert[(iv+1)%Nv]->x;
			for (int k=0; k<DIM; k++)
				e[iv][k]=y[k]-x[k];
			x=y;
		}
		//Compute volumes and cell-centers
		//Loop over the faces
		kill=surekill=0;
		vol=0.0;//volume of the cell
		for (int iv=0; iv<Nv; iv++)
		{	int
				iv1=(iv+1)%Nv,
				iv2=(iv+2)%Nv;
			DNode	*node=vert[iv];
			double
				dv,//volume increment
				xf[DIM],//face center coordinates
				a[DIM],area;//face area-vector
			y=node->var+ixold_loc;
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
			if (SCLP(a,e[iv])<0.0)
				for (double *p=a; p-a<DIM; p++) *p*=-.5;
			else
				for (double *p=a; p-a<DIM; p++) *p*= .5;
			area=LENGTH(a);
			if (area<SMALL){kill=surekill=1;break;}
			areas[iv]=area;
			dv=(SCLP(xf,a));
///DDD			if (!vert[iv]->state.boundary)
			if (!node->state.boundary)///DDD
			{	node->var[ivol]+=dv;//nodal scheme
				for (int i=0; i<DIM; i++)
					y[i]+=xf[i]*dv;
			}
			vol+=dv;
			if (neib[iv]==NULL)
			{//Compute centers of mass
			 // of boundary node surface elements
				for (int jv=0; jv<Nfv; jv++)
				{	int	j=(iv+jv+1)%Nv;
					double	*yy=vert[j]->var+ixold_loc;
					for (int i=0; i<DIM; i++)
						yy[i]+=xf[i]*area;
					vert[j]->var[ivol]+=area;//store area as a volume
					//if actual volumes are needed use
					// vert[j]->volume+=dv;
				}
			}
		}
		if (vol>=0.01*length_scale3)
		{
//			double	volin=0.17909/vol;
//			double	voli=1.0/vol;
//			for (int iv=0; iv<Nv; iv++)
//			{
				//Compare the current volume with
				//that of a right tetrahedron
				//with face area area[iv]
				//rightvolume=1/3*area*height
				//area=1/2*edge*face_height
				//face_height=(edge^2-(1/2*edge)^2)^(1/2)=3^(1/2)/2*edge
				//area=1/2*(3/4)^(1/2)*edge^2=3^(1/2)/4*edge^2
				//height^2=face_height^2-(face_height-face_height_base)^2
				//face_height_base=(edge^2-height^2)^(1/2)
				//height^2=3/4*edge^2-(3/4*edge^2-3^0.5*edge*(edge^2-height^2)^0.5+edge^2-height^2)
				//3/4-3/4+3^0.5*(1-(height/edge)^2)^0.5=1
				//3*(1-(height/edge)^2)=1
				//height/edge=(2/3)^0.5
				//right_volume=1/3*3^0.5/4*edge^2*(2/3)^0.5*edge
				//=2^0.5/12*edge^3
				//=2^0.5/12*(4/3^0.5)^0.5*area^(3/2)
				//=(2/3^0.5)^0.5/6*area^1.5
				//=0.17909*area^1.5
//				double	aspect_ratio=volin*pow(areas[iv],1.5);
//				if (aspect_ratio>6.0)///||aspect_ratio<0.06)
//				if (voli*areas[iv]>length_scalei)///||aspect_ratio<0.06)
//				{	kill=1; break;}
//			}
				for (int iv=0; iv<Nv; iv++)
				{	if (vol<0.01*length_scale*areas[iv])
					{	kill=1; break;}
				}
		}
		else
			surekill=kill=1;
//		if(kill)//cell volume too small
		if(kill)//cell volume too small
		{	//If the cell is at the boundary it should be deleted
			int	ivneib=-1,ivbnd=-1,nbound=0;
			fprintf(stderr,"Cell volume %g too small at \n",vol);
			for (int iv=0; iv<Nv; iv++)
			{	double	*x=vert[iv]->x;
				fprintf
				(	stderr,
					"\t%g, %g, %g, state=%d\n",
					x[0],x[1],x[2],vert[iv]->state.boundary
				);
			}
			for (int iv=0; iv<Nv; iv++)
				if (neib[iv]==NULL) nbound++;
			if (nbound>0)
			{//Reset neighbor cell pointers
				for (int iv=0; iv<Nv; iv++)
				{	DCell	*neib_cell=neib[iv];
					if (neib_cell!=NULL)
					{	ivneib=iv;
						DCell	**neibneib=(DCell**)neib_cell->neib;
						for (int jv=0; jv<Nv; jv++)
						if (neibneib[jv]==cell)
						{	neib_cell->neib[jv]=NULL;
							neib_cell->state.boundary=1;
							neib_cell->facetype[jv]=etype;
							break;
						}
					}
					else
						ivbnd=iv;
				}
				if (nbound==1)
				{	//Set node ivbnd as a boundary node
					vert[ivbnd]->state.boundary=1;
				}
				if (nbound==Nv1)
				{	//Delete node ivneib
					int found=0;
					DNode	*vneib=vert[ivneib];
					{//Check if the vertex belongs to any other cell
						DCell *dcell=cell->next;
						do
						{	DNode	**verts=dcell->vert;
							for (int i=0; i<Nv; i++)
								if (verts[i]==vneib)
								{	found=1; break;}
							dcell=dcell->next;
						}	while(!found&&dcell!=cell);
					}
					if (!found)
					{
						unlinkNodeListPointer(vneib,first_relax_dnode);
///					unlinkNodeListPointer(vneib,first_tool_boundary_dnode);
///					unlinkNodeListPointer(vneib,first_tool_dnode);
						domain->deleteElement(vert[ivneib],domain->dnode_root);
						nodes_deleted++;
						fprintf
						(	stderr,
							"Boundary vertex deleted, deleted vertexes = %d\n",
							nodes_deleted
						);FLUSH;
					}
				}
				//Delete cell
				domain->deleteElement(c->cell,domain->dcell_root);
///				unlinkCellListPointer(cell,first_tool_boundary_dcell);
				c->cell=NULL;
				cells_deleted++;
				fprintf
				(	stderr,
					"Boundary cell deleted, deleted cells = %d\n",
					cells_deleted
				);FLUSH;
			}
			else
			if (surekill)
			{	fprintf(stderr,"Can't delete internal cell: Aborting\n");
					exit(1);
			}
			else
			{	fprintf(stderr,"WARNING: Badaspect ratio of internal cell\n");
			}
		}
		if (nodes_deleted||cells_deleted)break;
		c=c->next;
	}	while(c!=first_relax_dcell);
		if (cells_deleted||nodes_deleted)
			updateLists(domain);
	}	while(nodes_deleted||cells_deleted);
	//Move nodes to new positions
	if (first_relax_dnode!=NULL)
	{	static const double
			//step=0.2*length_scale,
			sq_four_over_sq_three=2.0/pow(3.0,0.25),
			surface_relaxation=relaxation_factor,
			volume_relaxation=relaxation_factor;
		DNodeList *n=first_relax_dnode;
		do 
		{	DNode	*node=n->node;
			double
				//force[DIM],totforce,
				//factor,//normalization factor
				//*f=node->f,
				volume=node->var[ivol],
				*x=node->x,
				*y=node->var+ixold_loc;
			//for (int i=0; i<DIM; i++) force[i]=0.0;
			if (volume<=SMALL)
			{	fprintf
				(	stderr,
					"Tool::relax:ERROR: volume=%g too small at %g,%g,%g\n",
					volume,x[0],x[1],x[2]
				);
				exit(1);//goto endloop;
			}
			if (!node->state.boundary)
			{	double	volumei=0.75/volume;
				for (int i=0; i<DIM; i++)
					x[i]=volume_relaxation*x[i]+(1.0-volume_relaxation)*y[i]*volumei;
						//force[i]=volume_relaxation*(voli*y[i]-x[i]);
			}
			else//Boundary node
			if (node->state.toolsurface)
			{//On the tool surface
				double	dd,d[DIM],
					areai=1.0/volume;
				for (int i=0; i<DIM; i++)
				x[i]=surface_relaxation*x[i]+(1.0-surface_relaxation)*y[i]*areai;
				//put back to the tool surface
				dd=0.0;
				for (int i=0; i<DIM; i++)
				{	double	b=x[i]-toolpos[i];
					dd+=b*b;
					d[i]=b;
				}
				dd=tool_radius/sqrt(dd);
				for (int i=0; i<DIM; i++)
					x[i]=toolpos[i]+dd*d[i];
			}
			///	else//Boundary inside the tool
			///	{	//Apply the force that moves every
			///		// boundary node to the vertexes of 
			///		// the right triangle with the same area
			///		double	&area=volume,
			///			edge=sq_four_over_sq_three*sqrt(area);
			///		//if the distance between the nodes is
			///		// bigger than edge pool them closer
			///		// otherwise pool them apart
			///		///...
			///	}
			endloop: n=n->next;
		}	while (n!=first_relax_dnode);
	}
}
void	Tool::advanceBoundaryNodes
(
)
{	double
		toolpos[DIM],
		tool_insideradius=getinsideradius(),
		tool_radius=getradius(),
		length_scale=getlengthscale(),//max distance between the nodes
		optimal_separation=0.95*length_scale;
	getpos(toolpos);
	//Compute forces from the neighbor nodes
	if (first_tool_boundary_dnode!=NULL)
	for 
	(	DNodeList	*node=first_tool_boundary_dnode->next;
		; node=node->next
	)
	{	DNode	*dnode=node->node;
		double	f,
			*x=dnode->x,//coordinates of the node
			d[DIM],dd, //distance vector and its length
			fneib[DIM],//Force from neighbor nodes
			ftool[DIM];//Force from the tool
		for (int i=0; i<DIM; i++)fneib[i]=ftool[i]=0.0;
		//FORCE OF THE NEIGHBORS
		//Contribution from the boundary nodes
		for
		(	DNodeList	*neib=node->next;
			neib!=node;
			neib=neib->next
		)
		{	double
				neib_force,
				*y=neib->node->x;//coordinates of the neighbor
			dd=0.0;
			for (int i=0; i<DIM; i++)
			{	double	r=x[i]-y[i];
				dd+=r*r;
				d[i]=r;
			}
			dd=sqrt(dd);
			neib_force=force_neib(dd,0.3*length_scale);
			if (dd>SMALL)
				neib_force/=dd;
			for (int i=0; i<DIM; i++)
				fneib[i]+=neib_force*d[i];
		}
		//Contribution from the tool nodes
		if (first_tool_dnode!=NULL)
		for
		(	DNodeList	*neib=first_tool_dnode->next;
			;	neib=neib->next
		)
		{	double
				neib_force,
				*y=neib->node->x;//coordinates of the neighbor
			dd=0.0;
			for (int i=0; i<DIM; i++)
			{	double	r=x[i]-y[i];
				dd+=r*r;
				d[i]=r;
			}
			dd=sqrt(dd);
			neib_force=force_neib(dd,0.5*length_scale)/dd;
			if (dd>SMALL)
				neib_force/=dd;
			for (int i=0; i<DIM; i++)
				fneib[i]+=neib_force*d[i];
			if (neib==first_tool_dnode) break;
		}
		//FORCE OF THE TOOL
		//Distance from the tool center
		dd=0.0;
		for (int i=0; i<DIM; i++)
		{	double	r=x[i]-toolpos[i];
			dd+=r*r;
			d[i]=r;
		}
		dd=sqrt(dd);
		if (dd<SMALL) dd=SMALL;
		//MOVE THE NODE
		f=force_tool(dd,tool_radius)/dd;
		dd=0.0;
		for (int i=0; i<DIM; i++)
		{	double	r,dx=f*d[i]+fneib[i];
			x[i]+=dx;
			r=d[i]+dx;
			dd+=r*r;
			d[i]=r;
		}
		dd=sqrt(dd);
		if (dd>=tool_insideradius)
		{//Keep the node inside the tool
			double	r=1.0-tool_radius/dd;
			for (int i=0; i<DIM; i++)
				x[i]-=r*d[i];
			dnode->state.toolsurface=1;
		}
		if (node==first_tool_boundary_dnode) break;
	}
}
void	Tool::setpos(double x, double y, double z)
{
	center[0]=x;center[1]=y;center[2]=z;
}
void	Tool::setdefaults()
{
	ivolume=0;
	ixold=1;
	nbn=0;
	surface_roughness=SMALL;
	refinement=0.5;//0.7;
	type=coral;
	etype=boundary;
	mode.active=1;
	///mode.construct=1;
	///mode.destroy=0;
	state=add_cells;
	mode.manual=0;
	mode.setboundary=0;
	mode.setcells=0;
	mode.setnodes=0;
	rootpos=center;
	first_relax_dnode=NULL;
	first_tool_boundary_dnode=NULL;
	first_tool_dnode=NULL;
	root_cell=NULL;
	first_tool_boundary_dcell=NULL;
	first_relax_dcell=NULL;
	first_framenode=NULL;
	last_framenode=NULL;
	current_framenode=NULL;
	next=NULL;
	if (option.verbose)
	{	printf("Tool created, roughness=%g\n",surface_roughness);fflush(stdout);
	}
}
Tool::Tool()
{
	setdefaults();
}
Tool::Tool(char *filename)
{
	setdefaults();
	configfile=strdup(filename);
	load();
}
Tool::~Tool()
{
	deleteLists();
	deleteList(first_framenode);
	delete configfile;
}
void	Tool::load()
{	int	level=0;
	char	word[MAXLINLEN];
	double	roughness;
	FILE	*inp;
	if ((inp=fopen(configfile,"r"))==NULL)
		ERROR1("CAN'T OPEN ",configfile);
	first_tool_boundary_dnode=NULL;
	first_tool_boundary_dcell=NULL;
	first_tool_dnode=NULL;
	setradius(LARGE);
	setpos(0.0,0.0,0.0);
	while (!feof(inp))
	{	GETWORD(word,inp);
		if(strcmp(word,"shape")==0)
		{	GETWORD(word,inp);
			if(strcmp(word,"sphere")==0)
			{	shape=inside_sphere;
				for (int i=0; i<DIM; i++)
				{	double	coord;
					GETWORD(word,inp);
					sscanf(word,"%lg",&coord);
					center[i]=coord;
				}
				GETWORD(word,inp);
				sscanf(word,"%lg",&roughness);
				setroughness(roughness);
				if(option.verbose|option.debug)
				{	INDENT(1);
					printf
					(	"Tool: sphere, origin=%g,%g,%g, radius=%g\n",
						center[0],
						center[1],
						center[2],
						sphere.radius
					);fflush(stdout);
				}
			}
			else
			{	fprintf(stderr,"ERROR in file %s: Tool %s is unknown\n",configfile,word);
				exit(1);
			}
		}
		else
		if(strcmp(word,"state")==0)
		{	
			GETWORD(word,inp);
			if(strcmp(word,"build")==0)
				state=add_cells;
			else
			if(strcmp(word,"destroy")==0)
				state=remove_cells;
			else
			if(strcmp(word,"push")==0)
				state=push_cells;
			else
				state=passive;
			if (option.verbose)printf("\tTool state = %d\n",state);
		}
		else
		if(strcmp(word,"mode.setboundary")==0)
		{	GETWORD(word,inp);
			if(*word=='1'||*word=='Y'||*word=='y')
				mode.setboundary=1;
			else
				mode.setboundary=0;	
			if (option.verbose)printf("\tmode.build=%d\n",mode.setboundary);
		}
		else	break;
	}
	if(strcmp(word,FRAME_KEYWORD)==0)
	{	int 
			gotcircle=0,
///			shape_type=-1,
			build_mode=-1;
		GETWORD(word,inp);
		if (*word!='{')
		{	fprintf
			(	stderr,
				"ERROR: Keyword %s should be followed by }\n",
				FRAME_KEYWORD
			);
			exit(1);
		}
		level++;
		while(!feof(inp))
		{
			int	tool_type;
			double
				tool_radius,
				refinement_factor;
			Frame	*fp;
			GETWORD(word,inp);
			while (*word=='!') {SKIPLINE(inp);GETWORD(word,inp);}
			if (*word=='{') ERROR("Can't have { in frame description");
			if (*word=='}') {level--;break;}
			if (strcmp(word,"node")==0)
			{	///int	ip;
				double
					tool_radius,
					refinement_factor,
					x[DIM];
				GETWORD(word,inp);//TYPE OR BOUNDARY CONDITIONS
				sscanf(word,"%d",&tool_type);
				GETWORD(word,inp);//TOOL_RADIUS
				sscanf(word,"%lg",&tool_radius);
				GETWORD(word,inp);//CELL/TOOL SIZE RATIO
				sscanf(word,"%lg",&refinement_factor);
				GETWORD(word,inp);//TOOL COORDINATE 0
				sscanf(word,"%lg",x+0);
				GETWORD(word,inp);//TOOL COORDINATE 1
				sscanf(word,"%lg",x+1);
				GETWORD(word,inp);//TOOL COORDINATE 3
				sscanf(word,"%lg",x+2);
				addnode(x,tool_type,tool_radius,refinement_factor);
			}
			else
			if (strcmp(word,"segment")==0)
			{	int	n;
				double
					tool_radius,
					refinement_factor,
					d,dx,
					x0[DIM],x1[DIM];
				GETWORD(word,inp);//TYPE OR BOUNDARY CONDITIONS
				sscanf(word,"%d",&tool_type);
				GETWORD(word,inp);
				sscanf(word,"%lg",&tool_radius);
				GETWORD(word,inp);
				sscanf(word,"%lg",&refinement_factor);
				GETWORD(word,inp);
				sscanf(word,"%lg",x0+0);
				GETWORD(word,inp);
				sscanf(word,"%lg",x0+1);
				GETWORD(word,inp);
				sscanf(word,"%lg",x0+2);
				GETWORD(word,inp);
				sscanf(word,"%lg",x1+0);
				GETWORD(word,inp);
				sscanf(word,"%lg",x1+1);
				GETWORD(word,inp);
				sscanf(word,"%lg",x1+2);
				d=0.0;
				for (int i=0; i<DIM; i++)
				{	double	dd=x1[i]-x0[i];
					d+=dd*dd;
				}
				d=sqrt(d);
				dx=0.5*refinement_factor*tool_radius;
				n=(int)(d/dx)+1;
				dx=d/(double)n;
				for (int i=0; i<n; i++)
				{	double	x[DIM];
					for (int j=0; j<DIM; j++)
						x[j]=x0[j]+i*dx*(x1[j]-x0[j])/d;
					addnode(x,tool_type,tool_radius,refinement_factor);
				}
				if(option.verbose)
					printf("Number of tool frame nodes = %d\n",n);
			}
			else
			if (strcmp(word,"circle")==0)
			{	int	n;
				double
					tool_radius,
					refinement_factor,
					r,a,d,da,dx,
					pi=4.*atan(1.),
					x0[DIM],//center of the circle
					norm[DIM],//normal to the circle plane
					e0[DIM],//random unit vector
					e1[DIM];//orthonormal vectors in the plane of the circle
				GETWORD(word,inp);//TYPE OR BOUNDARY CONDITIONS
				sscanf(word,"%d",&tool_type);
				GETWORD(word,inp);
				sscanf(word,"%lg",&tool_radius);
				GETWORD(word,inp);
				sscanf(word,"%lg",&refinement_factor);
				GETWORD(word,inp);
				sscanf(word,"%lg",x0+0);
				GETWORD(word,inp);
				sscanf(word,"%lg",x0+1);
				GETWORD(word,inp);
				sscanf(word,"%lg",x0+2);
				GETWORD(word,inp);
				sscanf(word,"%lg",norm+0);
				GETWORD(word,inp);
				sscanf(word,"%lg",norm+1);
				GETWORD(word,inp);
				sscanf(word,"%lg",norm+2);
				GETWORD(word,inp);
				sscanf(word,"%lg",&r);
				d=2.*pi*r;
				dx=0.55*refinement_factor*tool_radius;///BUG: dx=0.6 endless loop?
				n=(int)ceil(d/dx);//number of points in the circle
				//dx+=d/(double)n-dx;
				da=2.0*pi/(double)n;
				if (gotcircle%2==0)da=-da;
				do//select random vector
				{	randvec(e0);
				}	while (fabs(SCLP(norm,e0))>0.9);
				VECP(e1,e0,norm);
				d=LENGTH(e1);
				for (int i=0; i<DIM; i++) e1[i]/=d;
				VECP(e0,e1,norm);
				d=LENGTH(e0);
				for (int i=0; i<DIM; i++) e0[i]/=d;
				//do the circle
				if (option.debug|option.verbose)
				{	printf
					(	"Tool: rad=%g; Circle: center=%g, %g, %g; R=%g\n",
						tool_radius,x0[0],x0[1],x0[2],r
					);
				}
				a=0.0;
				for (int i=0; i<n; i++)
				{	double	x[DIM],
						a=i*da,
						sa=sin(a),
						ca=cos(a);
					for (int j=0; j<DIM; j++)
						x[j]=x0[j]+r*(e0[j]*sa+e1[j]*ca);
					addnode(x,tool_type,tool_radius,refinement_factor);
				}
				gotcircle++;
			}
			else
			if (strcmp(word,"sector")==0)
			{	int	n,m;
				double
					tool_radius,
					refinement_factor,
					offset,direction,
					r,a,d,da,dx,
					pi=4.*atan(1.),
					x0[DIM],//center of the circle
					norm[DIM],//normal to the circle plane
					e0[DIM],//random unit vector
					e1[DIM];//orthonormal vectors in the plane of the circle
				GETWORD(word,inp);//TYPE OR BOUNDARY CONDITIONS
				sscanf(word,"%d",&tool_type);
				GETWORD(word,inp);
				sscanf(word,"%lg",&tool_radius);
				GETWORD(word,inp);
				sscanf(word,"%lg",&refinement_factor);
				GETWORD(word,inp);
				sscanf(word,"%lg",x0+0);
				GETWORD(word,inp);
				sscanf(word,"%lg",x0+1);
				GETWORD(word,inp);
				sscanf(word,"%lg",x0+2);
				GETWORD(word,inp);
				sscanf(word,"%lg",norm+0);
				GETWORD(word,inp);
				sscanf(word,"%lg",norm+1);
				GETWORD(word,inp);
				sscanf(word,"%lg",norm+2);
				GETWORD(word,inp);
				sscanf(word,"%lg",&r);
				GETWORD(word,inp);
				sscanf(word,"%lg",&offset);
				GETWORD(word,inp);
				sscanf(word,"%lg",&direction);
				d=2.*pi*r;
				dx=0.55*refinement_factor*tool_radius;///BUG: dx=0.6 endless loop?
				n=(int)ceil(d/dx);//number of points in the circle
				m=(int)(n*(offset/360.0));
				//dx+=d/(double)n-dx;
				da=2.0*pi/(double)n;
				if (direction==1)da=-da;
				do//select random vector
				{	randvec(e0);
				}	while (fabs(SCLP(norm,e0))>0.9);
				VECP(e1,e0,norm);
				d=LENGTH(e1);
				for (int i=0; i<DIM; i++) e1[i]/=d;
				VECP(e0,e1,norm);
				d=LENGTH(e0);
				for (int i=0; i<DIM; i++) e0[i]/=d;
				//do the circle
				if (option.debug|option.verbose)
				{	printf
					(	"Tool: rad=%g; Circle: center=%g, %g, %g; R=%g\n",
						tool_radius,x0[0],x0[1],x0[2],r
					);
				}
				a=0.0;
				for (int i=0; i<n; i++)
				{	double	x[DIM],
						a=da*((i+m)%n),
						sa=sin(a),
						ca=cos(a);
					for (int j=0; j<DIM; j++)
						x[j]=x0[j]+r*(e0[j]*sa+e1[j]*ca);
					addnode(x,tool_type,tool_radius,refinement_factor);
				}
			}
			else
			{	ERROR1("Unknown keyword '%s'\n",word);}
		}
		if (feof(inp)){ERROR1("FILE TOO SHORT",configfile);}
	}
	else
	{	fprintf
		(	stderr,
			"ERROR in %s: Unknown keyword '%s', expecting '%s'\n",
			configfile,word,FRAME_KEYWORD
		);
		exit(1);
	}
	fclose(inp);
	if (level!=0)
	{	fprintf
		(	stderr,
			"ERROR in %s: unmatched braces '{}'\n",
			configfile
		);
		exit(1);
	}
}
void	Tool::addnode
(
	double	*x,//coordinates
	int	tool_type,
	double	radius,
	double	refinement
)
{
	if (first_framenode==NULL)
	{	first_framenode=new Frame;
		last_framenode=first_framenode;
		first_framenode->next=last_framenode->next=first_framenode;
	}
	else
	{
		last_framenode->next=new Frame;
		last_framenode=last_framenode->next;
		last_framenode->next=first_framenode;
	}
	last_framenode->type=(ElementStatus)tool_type;
	last_framenode->radius=radius;
	last_framenode->refinement=refinement;
	for (int i=0; i<DIM; i++)
		last_framenode->x[i]=x[i];
}
void	Tool::move(double dx[])
{	for (int i=0; i<DIM; i++)
		center[i]+=dx[i];
}
void	Tool::pos(double x[])
{	for (int i=0; i<DIM; i++)
		rootpos[i]=x[i];
}
void	Tool::getpos(double x[])
{	for (int i=0; i<DIM; i++)
		x[i]=rootpos[i];
}
void	Tool::setype(ElementStatus new_type)
{
	etype=new_type;
}
void	Tool::setradius(double radius)
{
	sphere.radius=radius;
}
void	Tool::setroughness(double roughness)
{	const double	smallestroughness=0.0001;
	if (roughness<=0.0)roughness=smallestroughness;
	else
	if (roughness>=1.0)
	{	printf
		(	"WARNING: roughness parameter %g is greater than 1: set to 0\n",
			roughness
		);
		roughness=smallestroughness;
	}
	surface_roughness=roughness;
}
double	Tool::getradius()
{
	return sphere.radius;
}
double	Tool::getinsideradius()
{
	return (1.0-surface_roughness)*sphere.radius;
}
double	Tool::getrelaxradius()
{//Relaxation radius must be
 // greater than or equal to 
 // the tool radius
	//return 1.5*shape.radius;
	return 1.2*sphere.radius;
	//return 0.99*shape.radius;
}
double	Tool::getlengthscale()
{
	//return 0.7*shape.radius;
	return refinement*sphere.radius;
}
void	Tool::setrefinement(double	r)
{
	//return 0.7*shape.radius;
	if (r<0.5)
	{	fprintf(stderr,"WARNING: refinment factor %g is too small, increased to 0.5\n",r); 
		r=0.5;
	}
	refinement=r;
}
void	Tool::setVolInd(int ivol)
{
	ivolume=ivol;
}
void	Tool::setXoldInd(int jxold)
{
	ixold=jxold;
}
double	Tool::force_tool
(
	double	distance,
	double	interaction_radius
)
{
	//static const double	force=0.005;
	static const double	force=0.01*sphere.radius;
	if (distance<interaction_radius)
		return force*interaction_radius;
	return 0.0;
}
double	Tool::force_neib
(
	double	distance,
	double	interaction_radius
)
{
	//static const double	force=0.002;
	static const double	force=0.002*sphere.radius;
	if (distance<interaction_radius)
		return force;
	return 0.0;
}
void	Tool::insert(DCell	*newcell)
{
	DCellList	*p=new DCellList;
	p->next=first_tool_boundary_dcell->next;
	first_tool_boundary_dcell->next=p;
	p->cell=newcell;
}
template	<class List>
void	Tool::deleteList(List *&element)
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
void	Domain::addDCell
(	int	varbufsize,
 	int	ivert,
	int	jvert,
 	DNode	*&newvert,
//	DCellList	*current_boundary_cell,
	DCell	*&cell,
	int	iprev_neib_cell,
	int	ioldclone_neib_next,//local index of the next neighbor of the old clone
	DCell	*&oldclone,//neighbor cell's twin cell
	DCell	*&newclone
	//ElementStatus	new_type
)
{	DNode	**vert=cell->vert;
	DCell	**neib=(DCell**)cell->neib;
	DCellList	*p;
	//Create a new cell
	newclone=new DCell(cell,varbufsize);//clone cell to newcell
	//Insert cell into the DCell chain
	insert(cell,newclone);///-,dcell_root);
	vert[ivert]=newvert;
	newclone->vert[jvert]=newvert;
	neib[jvert]=newclone;
	newclone->neib[ivert]=cell;
	if (newclone->neib[jvert]!=NULL)
	{	int	i;
		DCell
			*cloneneib=(DCell*)newclone->neib[jvert],
			**nb=(DCell**)cloneneib->neib;
		for (i=0; i<Nv&&nb[i]!=cell; i++);
#ifdef DEBUG
		if (i==Nv)
		{	fprintf(stderr,"ERROR: addDCell: connectivity problem\n");
			fflush(stderr);
			exit(1);
		}
#endif
		nb[i]=newclone;
		//newclone->neib[jvert]->type[i]=internal;
	}
	//else
		//newclone->type[jvert]=new_type;
	if (oldclone!=NULL)
	{
		newclone->neib[iprev_neib_cell]=oldclone;
///newclone->facetype[iprev_neib_cell]=internal;///DDD
		oldclone->neib[ioldclone_neib_next]=newclone;
///oldclone->facetype[ioldclone_neib_next]=internal;///DDD
	}
	//else///BUG: produces a volume=0 error
		//newclone->type[iprev_neib_cell]=new_type;
	//Create a new boundary_cell
	//a separate loop later will reset
	//the actual status of this cell: 
	//(boundary, nearboundary, internal)
	//Insert the new reference to the new cell
	//into boundary cell pointer list 
	//before the first.
//		p=current_boundary_cell->next;
//		current_boundary_cell->next=new DCellList;
//		current_boundary_cell=current_boundary_cell->next;
//		current_boundary_cell->next=p;
//		current_boundary_cell->cell=newclone;
	tool->insert(newclone);
}
void	Tool::destroy
(
	Domain	*dom
)
{
	DNodeList	*tnode;
	if (first_tool_dnode==NULL) return;
	tnode=first_tool_dnode;
	do
	{	dom->deleteDNode(tnode->node);
		tnode=tnode->next;
	}	while (tnode!=first_tool_dnode);
}
void	Domain::deleteDNode(DNode *&dnode)
{
	Bond	*neib=dnode->bond;
	if (neib!=NULL)
	{	//Unlink the node
		for
		(	Bond	*neib_next=neib->next; 
			; neib_next=neib_next->next
		)
		{	DNode	*neib_node=neib_next->node;
///		int	nneib=dnode->nneib;
///		for (int	i=0; i<nneib; i++)
///		{	DNode *neib=dnode->neib[i];
///			if (neib!=NULL)//Remove neighbor pointers
///			{
			Bond	*neib_neib=neib_node->bond;
			for
			(	Bond	*neib_neib_next=neib_neib->next;
				; neib_neib_next=neib_neib_next->next
			)
			{	DNode	*neib_neib_node=neib_neib_next->node;
///				for (int j=0; j<neib->nneib; j++)
///				if (neib->neib[j]==dnode)
///				{	neib->neib[j]=NULL; 
				if (neib_neib_node==dnode)
				{	//delete neib_neib from the list
					if (neib_neib_next==neib_neib)
					{	//only one neighbor pointer
						delete neib_node->bond;
						neib_node->bond=NULL;
					}
					else
					{
						Bond	*tmp=neib_neib_next->next;
						delete neib_neib->next;
						neib_neib->next=tmp;
					}
///					neib_neib->node=NULL;
					break;
				}
				neib_neib=neib_neib_next;
///				dnode->neib[i]=NULL;
			}
			if (neib_next==neib)
			{
				delete dnode->bond;
				dnode->bond=NULL;
			}
			else
			{
				Bond	*tmp=neib_next->next;
				delete	neib->next;
				neib->next=tmp;
			}
			neib=neib_next;
		}
	}
	deleteElement(dnode,dnode_root);
}
template <class Element>
void	Domain::deleteElement(Element *&element, Element *&root)
{
	Element *prev,*next;
	if (element==NULL) return;
	prev=element->prev;
	next=element->next;
	prev->next=next;
	next->prev=prev;
	if (element==root) root=next;
	delete element;
	element=NULL;
}
///	void	Domain::deleteNode(DNode *&element)///DDD
///	{
///		DNode *prev,*next,
///			*root=dnode_root;
///		if (element==NULL) return;
///		prev=element->prev;
///		next=element->next;
///		if (element==root) root=next;
///		prev->next=next;
///		next->prev=prev;
///	///	delete element;
///	///	element=NULL;
///	}
///	void	Domain::deleteCell(DCell *&element)
///	{
///		DCell *prev,*next,
///			*root=dcell_root;
///		if (element==NULL) return;
///		prev=element->prev;
///		next=element->next;
///		if (element==root) root=next;
///		prev->next=next;
///		next->prev=prev;
///	///	delete element->var;
///	///	element->var=NULL;
///		delete element;///DDD
///		element=NULL;
///	}
///template <class Element>
///void	Domain::deleteElement(Element *&element)
///{
///	Element *prev,*next;
///	if (element==NULL) return;
///	prev=element->prev,
///	next=element->next;
///	prev->next=next;
///	next->prev=prev;
///	delete element;
///	element=NULL;
///}
void	Tool::unlinkNodeListPointer(DNode *kill, DNodeList *root)
{
	DNodeList	*node=root;
	if (root!=NULL)
	do
	{	if (node->node==kill)
			node->node=NULL;
		node=node->next;
	}	while(node!=root);
}
void	Tool::unlinkCellListPointer(DCell *kill, DCellList *root)
{
	DCellList	*cell=root;
	if (root!=NULL)
	do
	{	if (cell->cell==kill)
			cell->cell=NULL;
		cell=cell->next;
	}	while(cell!=root);
}

#include "templates.cc"

