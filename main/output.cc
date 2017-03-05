#include	<stdlib.h>
#include	<ctype.h>
#include	<string.h>
#include	<stdio.h>
#include	<math.h>

#include "main.h"
#include "io.h"
#include "geom.h"
#include "var.h"
#include "tool.h"
#include "domain.h"
#include "output.h"
namespace	Output
{
	int	nxwdump=0;//Window dump counter
	struct OutputStatus	status;
	char	ioconfigfile[MAXLINLEN];
	enum	OutputTypes	outype=default_output;
	char	*outypename[maxoutypes]=
	{	(char*)"Default",
		(char*)"User",
		(char*)"Tecplot",
		(char*)"Ensight"
	};
	void	initOutput()
	{
		status.active=0;
	}
	void	setOutput()
	{
		SCOPY(ioconfigfile,configfile);
	}
	void	outputBoundary()
	{
	}
	void	toggle()
	{
		status.active=status.active==1?0:1;
	}
	void	initxwd(int i)
	{
		nxwdump=i;
	}
	void	xwd(char *windowname)
	{
		int	sysreturned;
		char
			*s,
			dumpfilename[MAXLINLEN],
			command[MAXLINLEN];
		for (s=windowname; !isalpha(*s); s++);
		sprintf(dumpfilename,"%s%03d.xwd",s,nxwdump++);
		sprintf(command,"xwd -name %s -out %s",windowname,dumpfilename);
		if (option.verbose)
		{	printf("Dumping window to '%s' ... ",dumpfilename);fflush(stdout);
		}
		sysreturned=system(command);
		if (option.verbose) 
			printf("done\n");
		if (sysreturned!=0)
		{	fprintf(stderr,"Window dump to %s failed with status %d\n",sysreturned);
			fflush(stderr);
		}
	}
}
using namespace Output;
void	Domain::saveGeom(char *filename)
{
	printf("Saving geometry of domain %s in format: ",filename);
	switch(outype)
	{
		case default_output:
		printf("internal\n");
		saveGeomDefault(filename);
		break;
		case user_output:
		printf("user\n");
		saveGeomUser(filename);
		break;
		case tecplot_output:
		printf("tecplot\n");
		saveGeomTecplot(filename);
		break;
		case ensight_output:
		printf("ensight\n");
		saveGeomEnsight(filename);
		break;
		default:
		WARNING("WRONG OUTPUT TYPE SPECIFIED: FILE-OUTPUT SKIPPED");
	}
}
void	Domain::saveGeomEnsight(char *filename)
{

}
void	Domain::saveGeomTecplot(char *filename)
{

}
void	Domain::saveGeomUser(char *filename)
{

}
void	Domain::saveGeomDefault(char *filename)
{	char outgeofilename[MAXLINLEN];
	FILE	*out;
	if (type!=dynamic)
	{	fprintf
		(	stderr,
			"ERROR: Only dynamic domain output is possible: write operation skipped\n"
		);
		return;
	}
#ifdef DEBUG
	{//Check connectivity
		int	i=0;
		DNode	*node=dnode_root;
		DCell	*cell=dcell_root;
	do//Zero out source terms
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
		{	fprintf(stderr,"SaveGeom:: grid inconsistency: loose node %d\n",i++);
			exit(1);
		}
		node=node->next;
	}	while(node!=dnode_root);
	}
#endif
	cleanInternalFaces();
	sprintf(outgeofilename,"%sgeo.dat",filename);
	OPENWRITE(outgeofilename,out);
	if (option.verbose)
	{	printf("File %s opened for writing\n",outgeofilename);
		printf("\tDomain type = %d\n",type);FLUSH;
	}
	fwrite(&type,sizeof(type),1,out);
	if (dnode_root!=NULL)
	{	//Count the number of nodes and set indexes
		int	nn=0;
		DNode	*node=dnode_root;
		do
		{	node->index=nn++;
			node=node->next;
		}	while (node!=dnode_root);
		//Output number of nodes
		fwrite(&nn,sizeof(nn),1,out);
		//Output nodes
		node=dnode_root;
		do
		{	//Write type
			//int	itype=(int)node->type;
			struct NodeState	state=node->state;
			state.insidetool=0;
			state.toolsurface=0;
			state.relax=0;
			fwrite(&state,sizeof(node->state),1,out);
			//fwrite(&node->nneib,sizeof(node->nneib),1,out);
			//Coordinates
			fwrite(node->x,sizeof(node->x[0]),DIM,out);
			node=node->next;
		}	while (node!=dnode_root);
	}
	if (dcell_root!=NULL)
	{	//Store cell data
		const	int	minus_one=-1;
		int nc=0;
		DCell	*cell=dcell_root;
		//Count number of cells
		do
		{	cell->index=nc++;
			cell=cell->next;
		}	while (cell!=dcell_root);
		//Write the number of cells
		fwrite(&nc,sizeof(nc),1,out);
		//Write cell-node connectivity
		cell=dcell_root;
		do
		{	//write cell-type
			fwrite(&cell->facetype,sizeof(cell->facetype[0]),Nf,out);
			//write cell vertex-nodes
			for (int i=0; i<Nv; i++)
			if (cell->vert[i]!=NULL)
				fwrite(&cell->vert[i]->index,sizeof(int),1,out);
#ifdef DEBUG
			else
			{	fprintf(stderr,"saveGeom: cell-node connectivity error\n");
				exit(1);
			}
#endif
			cell=cell->next;
		}	while (cell!=dcell_root);
		//Write cell-cell connectivity
		cell=dcell_root;
		do
		{	for (int i=0; i<Nf; i++)
			{	DCell	*cellneib=(DCell*)cell->neib[i];
				if (cellneib!=NULL&&cell->facetype[i]!=composite)
					fwrite(&cellneib->index,sizeof(int),1,out);
				else
					fwrite(&minus_one,sizeof(minus_one),1,out);
			}
			cell=cell->next;
		}	while (cell!=dcell_root);
	}
	fclose(out);
	if (option.verbose)
	{	printf("Writing %s done\n",filename);FLUSH;
	}	
	indexNeibNodes();//reset indexes
}
//SET OUTPUT
void	Domain::setOutputDefault(){outype=default_output;}
void	Domain::setOutputEnsight(){outype=ensight_output;}
void	Domain::setOutputTecplot(){outype=tecplot_output;}
void	Domain::setOutputUser(){outype=user_output;}
// OUTPUT ALL DATA
void	Domain::outputData(char *filename)
{
	if(option.verbose)printf("Writing data to %s in %s format ...",filename,outypename[outype]);
	switch(outype)
	{
		case default_output:
		outputDefault(filename);
		break;
		case user_output:
		outputUser(filename);
		break;
		case tecplot_output:
		outputTecplot(filename);
		break;
		case ensight_output:
		outputEnsight(filename);
		break;
		default:
		WARNING("WRONG OUTPUT TYPE SPECIFIED: FILE-OUTPUT SKIPPED");
	}
	if(option.verbose)printf(" done\n");
};
void	Domain::saveDispVarDefault(int ivar, char *filename)
{
	char tecfilename[MAXLINLEN];
	DNode	*node_root=dnode_root,
		*node=node_root;
	DCell	*cell_root=dcell_root,
		*cell=cell_root;
	int	n=0;
	if(ivar<0)
	{	fprintf(stderr,"No active variable defined: write operation skipped\n");
		return;
	}
	printf("Saving variable %s: to %s in default format\n",variable[ivar].name,filename);///DDD
	do//Count the nodes
	{	n++;
		node=node->next;
	}	while(node!=node_root);
	sprintf(tecfilename,"%s%dtec.dat",filename,ivar);
	FILE	*out=fopen(tecfilename,"w");
	printf
	(	"Writing %d variable on %d nodes to %s at time %g in Tecplot format ",
		ivar,n,tecfilename,runtime.current
	);FLUSH;
	if(node_root==NULL) return;
	//Output header
	fprintf(out,"Time= %11.4f\n",runtime.current);
	fprintf(out,"VAR=");
	//for(int i=0; i<nvar; i++)
	{	int	i=ivar,
			varloc=variable[i].loc,
			varlen=variable[i].dimension;
		char *name=variable[i].name;
		if(varlen>1)
		for(int j=0;j<varlen;j++)
			fprintf(out,"%s%d ",name,j);
		else
			fprintf(out,"%s ",name);
		fprintf(out,"\n");
	}
	//Output data
	n=0;
	do//Output nodal variables
	{	double
			*x=node->x,
			*var=node->var;
		fprintf(out,"%d\t%g\t%g\t%g",n++,x[0],x[1],x[2]);
		//for(int i=0; i<nvar; i++)
		{	int	i=ivar,
				varloc=variable[i].loc,
				lenvar=variable[i].dimension;
			double	*v=var+varloc;
			for(int j=0; j<lenvar; j++)
			{
				fprintf(out,"\t%g",v[j]);
			}
			fprintf(out,"\n");
		}
		node=node->next;
	}	while(node!=node_root);
	fclose(out);
}
void	Domain::outputDefault(char *filename)
{
printf("outputDefault: %s",filename);///DDD
}
void	Domain::outputEnsight(char *filename)
{
printf("outputEnsight",filename);///DDD
}

void	Domain::outputTecplot(char *filename)
{
	char tecfilename[MAXLINLEN];
	DNode	*node_root=dnode_root,
		*node=node_root;
	DCell	*cell_root=dcell_root,
		*cell=cell_root;
	int	n=0;
	do//Count the nodes
	{	n++;
		node=node->next;
	}	while(node!=node_root);
	sprintf(tecfilename,"%stec.dat",filename);
	FILE	*out=fopen(tecfilename,"w");
	printf
	(	"Writing %d variables on %d nodes to %s at time %g in Tecplot format ",
		nvar,n,tecfilename,runtime.current
	);FLUSH;
	if(node_root==NULL) return;
	//Output header
	fprintf(out,"VAR=");
	for(int i=0; i<nvar; i++)
	{	int
			varloc=variable[i].loc,
			varlen=variable[i].dimension;
		char *name=variable[i].name;
		if(varlen>1)
		for(int j=0;j<varlen;j++)
			fprintf(out,"%s%d ",name,j);
		else
			fprintf(out,"%s ",name);
		fprintf(out,"\n");
	}
	//Output data
	n=0;
	do//Output nodal variables
	{	double
			*x=node->x,
			*var=node->var;
		fprintf(out,"%d\t%g\t%g\t%g",n++,x[0],x[1],x[2]);
		for(int i=0; i<nvar; i++)
		{	int
				varloc=variable[i].loc,
				lenvar=variable[i].dimension;
			double	*v=var+varloc;
			for(int j=0; j<lenvar; j++)
				fprintf(out,"\t%g",v[j]);
			fprintf(out,"\n");
		}
		node=node->next;
	}	while(node!=node_root);
	fclose(out);
}
void	Domain::outputUser(char *filename)
{
printf("outputUser",filename);///DDD
}

