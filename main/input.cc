#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include "vecalg.h"
#include "main.h"
#include "io.h"
#include "geom.h"
#include "var.h"
#include "tool.h"
#include "domain.h"
#include "box.h"
#include "cylinder.h"
#include "input.h"
#include "output.h"
#include "templates.cc"


void	Domain::setGeom
(
	char inpfilename[]
)
{
	if(strstr(inpfilename,".xml\0")!=NULL)
		setGeomXML(inpfilename);
	else
		setGeomCFG(inpfilename);
}

void	Domain::setGeomCFG
(
	char inpfilename[]
)
{	int	
///		type_specified=0,
		gotword=0,
		idomain=this-domain_root;//,ivar;
	char	*p,word[MAXLINLEN];
	FILE	*inp;
	//Corrdinates of basic element types
	coordinates=new Variable[maxelements];
	volumes=new Variable[maxelements];
	for (int itype=0; itype<maxelements; itype++)
	{	coordinates[itype].type=itype;
		coordinates[itype].dimension=DIM;
		volumes[itype].type=itype;
		volumes[itype].dimension=1;
		Ne[itype]=0;
	}
	if (option.verbose|option.debug)
	{	printf
		(	"Reading geometry for domain no. %d: from %s\n",
			idomain+1,inpfilename
		);FLUSH;
	}
	//Scan 1: Setting up geometry and allocating variables
	//Locate the domain
	OPENREAD(inpfilename,inp);
	//Skip first idomain domains
	level=0; 
	while (!feof(inp)&&idomain>=0&&level>=0)
	{	GETWORD(word,inp);
		if (*word=='!') {SKIPLINE(inp);continue;}
		if (*word=='{') level++;
		if (*word=='}') level--;
		if (level==0&&strcmp(word,(char*)DOMAIN_KEYWORD)==0)
		{	idomain--;
			GETWORD(word,inp);
			SCOPY(name,word);//domain name as found in configfile
			//if (strcmp(word,name)==0) break;
		}
	}
	idomain=this-domain_root;
	if (feof(inp)) 
	{	fprintf
		(	stderr,
			"Domain name %s not found in %s\n",
			name,inpfilename
		);exit(1);
	}
	//GETWORD(word,inp);//Get domain name
	//for (p=word; isspace(*p); p++);
	//while (!isspace(*p))p++;*p='\0';
	//SCOPY(name,word);FLUSH;
	if (option.verbose|option.debug)
	{	INDENT(0);printf("%s name='%s'\n",(char*)DOMAIN_KEYWORD,name);FLUSH;
	}
	while (!feof(inp)&&*word!='{')
		GETWORD(word,inp);
	level++;
	//Scan domain description
	while (!feof(inp)&&level>=1)
	{	if (!gotword) GETWORD(word,inp);
		gotword=0;
		if (*word=='!') {SKIPLINE(inp);continue;}
		if (*word=='{') level++;
		if (*word=='}') level--; 
		if (level==0) break;
		if (level==1)
		{
			if (strcmp(word,(char*)MODEL_KEYWORD)==0)
			{	GETWORD(word,inp);
				model=atoi(word);
				if(option.verbose)
					printf("Model %d\n",model);
			}
		/*added by zhang, missing in original, I guess*/
#ifdef WITH_DOVE
			else
			if (strcmp(word,"proc")==0)
			{	GETWORD(word,inp);
				iproc=atoi(word);
				if(option.verbose)
					printf("iproc: %d\n",iproc);
			}
#endif
		/***********/
			else
			if (strcmp(word,(char*)DOMAIN_TYPE_KEYWORD)==0)
			{///	type_specified=1;
				GETWORD(word,inp);
//?				if (strcmp(word,"box")==0)
//?				{	int	
//?						gotdimensions=0,
//?						gotgridsize=0,
//?						gotboundary=0;
//?					Box	b;
//?					type=box;
//?					if (option.verbose|option.debug)
//?					{	INDENT(0);printf("Type = box\n");
//?					}
//?					for (int i=0; i<3 && !feof(inp); i++)
//?					{//READ DIMENSIONS, GRID AND BOUNDARY
//?						GETWORD(word,inp);
//?						if(strcmp(word,DIMENSIONS_TAG)==0)
//?						{	if(option.verbose|option.debug)
//?							{	INDENT(1);printf("%s",DIMENSIONS_TAG);
//?							}
//?							for (int i=0; i<DIM; i++)
//?							{
//?								GETWORD(word,inp);
//?								sscanf(word,"%lg",xmin+i);
//?								GETWORD(word,inp);
//?								sscanf(word,"%lg",xmax+i);
//?								if(option.verbose|option.debug)
//?								{	printf("\t%lg : %lg",xmin[i],xmax[i]);FLUSH;
//?								}
//?							}
//?							if(option.verbose|option.debug)printf("\n");
//?							gotdimensions=1;
//?							fixedlimits=1;
//?						}
//?						if (strcmp(word,GRIDSIZE_KEYWORD)==0)
//?						{	if(option.verbose|option.debug)
//?							{	INDENT(1);printf("%s:",GRIDSIZE_KEYWORD);
//?							}
//?							for (int i=0; i<DIM; i++)
//?							{	char	*p;
//?								GETWORD(word,inp);
//?								sscanf(word,"%ld",b.nx+i);
//?								if(option.verbose|option.debug)printf("\t%ld ",b.nx[i]);FLUSH;
//?							}
//?							if(option.verbose|option.debug)printf("\n");
//?							gotgridsize=-1;
//?						}
//?						if (strcmp(word,BOUNDARY_CONDITIONS_TAG)==0)
//?						{	if(option.verbose|option.debug)
//?							{	INDENT(1);printf("%s:",BOUNDARY_CONDITIONS_TAG);
//?							}
//?							for (int i=0;i<6;i++)
//?							{	GETWORD(word,inp);
//?								sscanf(word,"%d",b.bc+i);
//?								if(option.verbose|option.debug)
//?									printf(" %d",b.bc[i]);
//?								b.bc[i]=-b.bc[i];//stored as negative numbers
//?							}
//?							if(option.verbose|option.debug)printf("\n");
//?							gotboundary=1;
//?						}
//?					}
//?					if 
//?					(    !gotdimensions
//?						|| !gotgridsize
//?						|| !gotboundary
//?					)	
//?					{	fprintf
//?						(	stderr,
//?							"One of the keywords not found in the description of domain %d: %s=%d or %s=%d or %s=%d\n",
//?							name,
//?							DIMENSIONS_TAG,gotdimensions,
//?							GRIDSIZE_KEYWORD,gotgridsize,
//?							BOUNDARY_CONDITIONS_TAG,gotboundary
//?						); exit(1);
//?					}
//?					Ne[nodes]=1;
//?					for (int i=0; i<DIM; i++) Ne[nodes]*=b.nx[i]+1;
//?					Ne[edges]=0;
//?					Ne[cells]=1;
//?					for (int i=0; i<DIM; i++) Ne[cells]*=b.nx[i];
//?					Ne[cells]*=5;
//?					nbv=2*((b.nx[0]+1)*(b.nx[1]+1)+(b.nx[0]+b.nx[1])*(b.nx[2]-1));
//?					Ne[boundary_nodes]=nbv;
//?					nbf=0;
//?					for (int i=0; i<DIM; i++) nbf+=4*b.nx[(i+1)%DIM]*b.nx[(i+2)%DIM];
//?					{	//Determine the number of faces
//?						int nif=4*Ne[cells]-nbf;//number_of_internal_faces
//?						if (nif%2==1)ERROR("Grid inconsistency in the number of faces");
//?						nif/=2;
//?						Ne[boundary_faces]=nbf;
//?						Ne[faces]=nif+nbf;
//?					}
//?					//Ne[points]=0;
//?					if (option.verbose|option.debug)
//?					{
//?						for (int itype=0; itype<maxelements; itype++)
//?						{	INDENT(1);
//?							printf("Number of %s: %d\n",elementype[itype],Ne[itype]);
//?						}
//?					}
//?					alloc(nodes,node);
//?					//alloc(edges,edge);
//?					alloc(faces,face);
//?					alloc(cells,cell);
//?					allocgeovar(coordinates);
//?					allocgeovar(volumes);
//?					b.create(xmin,xmax,coordinates[nodes].val,node);
//?					//nbv=b.setBoundaryVertexes(node);
//?					//if (option.verbose|option.debug)
//?					//{	printf("Number of boundary vertexes: nbv=%d\n",nbv);FLUSH;
//?					//}
//?					b.connect(cell);
//?				}
//?				else
//?				if (strcmp(word,"cylinder")==0)
//?				{	int	
//?						gotdimensions=0,
//?						gotgridsize=0,
//?						gotboundary=0,
//?						nx[2],//number of rings (nx[0]), and slices (nx[1])
//?						bc[3];//boundary conditions on cylinder walls:
//?						      //bc[0]: wall, bc[1]: bottom, bc[2]: top
//?					double	dd[2];//dd[0]: cylinder radius; dd[1]: cylinder length
//?					Cylinder	cyl;
//?					type=cylinder;
//?					if (option.verbose|option.debug)
//?					{	INDENT(0);printf("Type = cylinder\n");
//?					}
//?					for (int item=0; item<3 && !feof(inp); item++)
//?					{	GETWORD(word,inp);
//?						if(*word=='!')
//?						{	SKIPLINE(inp);item--;continue;
//?						}
//?						if(strcmp(word,DIMENSIONS_TAG)==0)
//?						{	if(option.verbose|option.debug)
//?							{	INDENT(1);printf("%s",DIMENSIONS_TAG);
//?							}
//?							GETWORD(word,inp);//length
//?							sscanf(word,"%lg",&dd[0]);
//?							GETWORD(word,inp);//radius
//?							sscanf(word,"%lg",&dd[1]);
//?							if(option.verbose|option.debug)
//?							{	printf("\tlenth=%lg, radius=%lg\n",dd[0],dd[1]);FLUSH;
//?							}
//?							//Compute limits
//?							//	for (int i=0; i<2; i++)
//?							//	{	char	*p;
//?							//		GETWORD(word,inp);
//?							//		sscanf(word,"%lg",xmin+i);
//?							//		GETWORD(word,inp);
//?							//		sscanf(word,"%lg",xmax+i);
//?							//		if(option.verbose|option.debug)
//?							//		{	printf("\t%lg : %lg",xmin[i],xmax[i]);FLUSH;
//?							//		}
//?							//	}
//?							//	if(option.verbose|option.debug)printf("\n");
//?							gotdimensions=1;
//?							fixedlimits=1;
//?							continue;
//?						}
//?						if (strcmp(word,GRIDSIZE_KEYWORD)==0)
//?						{	if(option.verbose|option.debug)
//?							{	INDENT(1);printf("%s:",GRIDSIZE_KEYWORD);
//?							}
//?							GETWORD(word,inp);//get no. rings
//?							sscanf(word,"%ld",&nx[0]);
//?							if(option.verbose|option.debug)printf(" %ld rings, ",nx[0]);FLUSH;
//?							GETWORD(word,inp);//get no. slices
//?							sscanf(word,"%ld",&nx[1]);
//?							if(option.verbose|option.debug)printf("%ld slices\n",nx[1]);FLUSH;
//?							///	for (int i=0; i<DIM; i++)
//?							///	{	char	*p;
//?							///		GETWORD(word,inp);
//?							///		sscanf(word,"%ld",b.nx+i);
//?							///		if(option.verbose|option.debug)printf("\t%ld ",b.nx[i]);FLUSH;
//?							///	}
//?							///	if(option.verbose|option.debug)printf("\n");
//?							cyl.setgrid(nx);
//?							gotgridsize=1;
//?							continue;
//?						}
//?						if (strcmp(word,BOUNDARY_CONDITIONS_TAG)==0)
//?						{	if(option.verbose|option.debug)
//?							{	INDENT(1);printf("%s (wall, bottom, top):",BOUNDARY_CONDITIONS_TAG);
//?							}
//?							for (int i=0;i<3;i++)
//?							{	GETWORD(word,inp);
//?								sscanf(word,"%d",bc+i);
//?								if(option.verbose|option.debug)
//?									printf(" %d",bc[i]);
//?								bc[i]=-bc[i];//stored as negative numbers
//?								cyl.setbc(bc);
//?							}
//?							if(option.verbose|option.debug)printf("\n");
//?							gotboundary=1;
//?							continue;
//?						}
//?					}
//?					if 
//?					(    !gotdimensions
//?						|| !gotgridsize
//?						|| !gotboundary
//?					)	
//?					{	fprintf
//?						(	stderr,
//?							"One of the keywords not found in the description of domain %s: %s=%d or %s=%d or %s=%d\n",
//?							name,
//?							DIMENSIONS_TAG,gotdimensions,
//?							GRIDSIZE_KEYWORD,gotgridsize,
//?							BOUNDARY_CONDITIONS_TAG,gotboundary
//?						); exit(1);
//?					}
//?					cyl.init(nx,dd);
//?					Ne[nodes]=cyl.getnodes();
//?					///for (int i=0; i<DIM; i++) Ne[nodes]*=b.nx[i]+1;
//?					Ne[edges]=0;
//?					Ne[cells]=cyl.getcells();
//?					///for (int i=0; i<DIM; i++) Ne[cells]*=b.nx[i];
//?					///Ne[cells]*=5;
//?					//Count boundary faces ...
//?					//nbf=0;
//?					//for (int i=0; i<DIM; i++) nbf+=4*b.nx[(i+1)%DIM]*b.nx[(i+2)%DIM];
//?					{//Determine the number of boundary faces
//?						int	
//?							nrings=nx[0],
//?							nslices=nx[1],
//?							nringfaces,nsectorfaces,
//?							nwallfaces,nwallsectorfaces;
//?						nsectorfaces=nrings*nrings;
//?						nringfaces=6*nsectorfaces;
//?						nwallsectorfaces=2*nrings;
//?						nwallfaces=6*nwallsectorfaces*nslices;
//?						nbf=2*nringfaces+nwallfaces;
//?					}
//?					{	//Determine the number of faces
//?						int nif=4*Ne[cells]-nbf;//number_of_internal_faces
//?						if (nif%2==1)ERROR("Grid inconsistency in the number of faces");
//?						nif/=2;
//?						Ne[boundary_faces]=nbf;
//?						Ne[faces]=nif+nbf;
//?					}
//?					if (option.verbose|option.debug)
//?					{
//?						for (int itype=0; itype<maxelements; itype++)
//?						{	INDENT(1);
//?							printf("Number of %s: %d\n",elementype[itype],Ne[itype]);
//?						}
//?					}
//?					alloc(nodes,node);
//?					//alloc(edges,edge);
//?					alloc(faces,face);
//?					alloc(cells,cell);
//?					allocgeovar(coordinates);
//?					allocgeovar(volumes);
//?					cyl.create(coordinates[nodes].val,node,cell);
//?					//nbv=b.setBoundaryVertexes(node);
//?					//if (option.verbose|option.debug)
//?					//{	printf("Number of boundary vertexes: nbv=%d\n",nbv);FLUSH;
//?					//}
//?					//cyl.connect(cell);
//?				}
//?				else
				if (strcmp(word,"dynamic")==0)
				{	char	geofile[MAXLINLEN];
					enum GridType	gridtype=unstructured_tetrahedral_grid;
					int
						///gotdimensions=0,
						///gotmaxnodes=0,
						///gotmaxcells=0,
						gotnodes=0,
						gotcells=0,
						gotgrid=0;
					type=dynamic;
					//first_framenode=NULL;
					if(option.verbose|option.debug)
					{
						INDENT(0);printf("\tType = dynamic\n");
					}
					GETWORD(word,inp);
					if(strcmp(word,DIMENSIONS_TAG)==0)
					{	if(option.verbose|option.debug)
						{	INDENT(1);printf("%s",DIMENSIONS_TAG);
						}
						for (int i=0; i<DIM; i++)
						{
							GETWORD(word,inp);
							sscanf(word,"%lg",xmin+i);
							GETWORD(word,inp);
							sscanf(word,"%lg",xmax+i);
							if(option.verbose|option.debug)
							{	printf("\t%lg : %lg",xmin[i],xmax[i]);FLUSH;
							}
						}
						if(option.verbose|option.debug)printf("\n");
						fixedlimits=1;
					}
					else ERROR1((char*)"KEYWORD NOT FOUND: ",DIMENSIONS_TAG);
					GETWORD(word,inp);
					if(strcmp(word,(char*)GRIDTYPE_KEYWORD)==0)
					{	//Get grid-type
						GETWORD(word,inp);
						if(strcmp(word,"box")==0)///gridtype=box;
						{	int	
								gotdimensions=0,
								gotgridsize=0,
								gotboundary=0;
							Box	b;
							if(option.debug|option.verbose)
							{	printf("\t%s=box",(char*)GRIDTYPE_KEYWORD);FLUSH;
							}
							for (int i=0; i<3 && !feof(inp); i++)
							{//READ BOX DIMENSIONS, GRID AND BOUNDARY
								GETWORD(word,inp);
								if(strcmp(word,GRIDDIMENSIONS_TAG)==0)
								{	if(option.verbose|option.debug)
									{	INDENT(1);printf("%s",GRIDDIMENSIONS_TAG);
									}
									for (int i=0; i<DIM; i++)
									{
										GETWORD(word,inp);
										sscanf(word,"%lg",xmin+i);
										GETWORD(word,inp);
										sscanf(word,"%lg",xmax+i);
										if(option.verbose|option.debug)
										{	printf("\t%lg : %lg",xmin[i],xmax[i]);FLUSH;
										}
									}
									if(option.verbose|option.debug)printf("\n");
									gotdimensions=1;
									fixedlimits=1;
								}
								if (strcmp(word,(char*)GRIDSIZE_KEYWORD)==0)
								{	if(option.verbose|option.debug)
									{	INDENT(1);printf("%s:",GRIDSIZE_KEYWORD);
									}
									for (int i=0; i<DIM; i++)
									{	char	*p;
										GETWORD(word,inp);
										sscanf(word,"%ld",b.nx+i);
										if(option.verbose|option.debug)printf("\t%ld ",b.nx[i]);FLUSH;
									}
									if(option.verbose|option.debug)printf("\n");
									gotgridsize=-1;
								}
								if (strcmp(word,BOUNDARY_CONDITIONS_TAG)==0)
								{	if(option.verbose|option.debug)
									{	INDENT(1);printf("%s:",BOUNDARY_CONDITIONS_TAG);
									}
									for (int i=0;i<6;i++)
									{	GETWORD(word,inp);
										sscanf(word,"%d",b.bc+i);
										if(option.verbose|option.debug)
											printf(" %d",b.bc[i]);
										b.bc[i]=-b.bc[i];//stored as negative numbers
									}
									if(option.verbose|option.debug)printf("\n");
									gotboundary=1;
								}
							}
							if 
							(    !gotdimensions
								|| !gotgridsize
								|| !gotboundary
							)	
							{	fprintf
								(	stderr,
									"One of the keywords not found in the description of domain %d: %s=%d or %s=%d or %s=%d\n",
									name,
									GRIDDIMENSIONS_TAG,gotdimensions,
									GRIDSIZE_KEYWORD,gotgridsize,
									BOUNDARY_CONDITIONS_TAG,gotboundary
								); exit(1);
							}
							Ne[nodes]=1;
							for (int i=0; i<DIM; i++) Ne[nodes]*=b.nx[i]+1;
							Ne[edges]=0;
							Ne[cells]=1;
							for (int i=0; i<DIM; i++) Ne[cells]*=b.nx[i];
							Ne[cells]*=5;
							nbv=2*((b.nx[0]+1)*(b.nx[1]+1)+(b.nx[0]+b.nx[1])*(b.nx[2]-1));
							Ne[boundary_nodes]=nbv;
							nbf=0;
							for (int i=0; i<DIM; i++) nbf+=4*b.nx[(i+1)%DIM]*b.nx[(i+2)%DIM];
							{	//Determine the number of faces
								int nif=4*Ne[cells]-nbf;//number_of_internal_faces
								if (nif%2==1)ERROR("Grid inconsistency in the number of faces");
								nif/=2;
								Ne[boundary_faces]=nbf;
								Ne[faces]=nif+nbf;
							}
							if (option.verbose|option.debug)
							{
								for (int itype=0; itype<maxelements; itype++)
								{	INDENT(1);
									printf("Number of %s: %d\n",elementype[itype],Ne[itype]);
								}
							}
							if(option.verbose){printf("Allocating %d brick nodes\n",Ne[nodes]);FLUSH;}
							alloc(Ne[nodes],nodes,dnode_root);
							if(option.verbose){printf("Allocating %d brick cells\n",Ne[cells]);FLUSH;}
							alloc(Ne[cells],cells,dcell_root);
							b.create(xmin,xmax,dnode_root);
							b.connect(dnode_root,dcell_root);
						}
						else
						if(strcmp(word,"cylinder")==0)///gridtype=cylinder;
//void	read_cylinder
//(
//	char	*&word,
//	FILE	*inp
//)
				{	int	
						Ne[maxelements],
						gotdimensions=0,
						gotgridsize=0,
						gotboundary=0,
						nx[2],//number of rings (nx[0]), and slices (nx[1])
						bc[3];//boundary conditions on cylinder walls:
						      //bc[0]: wall, bc[1]: bottom, bc[2]: top
					double	dd[2];//dd[0]: cylinder radius; dd[1]: cylinder length
					Cylinder	cyl;
					if (option.verbose|option.debug)
					{	INDENT(0);printf("\t%s=cylinder\n",GRIDTYPE_KEYWORD);
					}
					for (int i=0; i<3 && !feof(inp); i++)
					{	GETWORD(word,inp);
						if(*word=='!')
						{	SKIPLINE(inp);i--;continue;
						}
						if(strcmp(word,GRIDDIMENSIONS_TAG)==0)
						{	if(option.verbose|option.debug)
							{	INDENT(1);printf("%s",GRIDDIMENSIONS_TAG);
							}
							GETWORD(word,inp);//length
							sscanf(word,"%lg",&dd[0]);
							GETWORD(word,inp);//radius
							sscanf(word,"%lg",&dd[1]);
							if(option.verbose|option.debug)
							{	printf("\tlenth=%lg, radius=%lg\n",dd[0],dd[1]);FLUSH;
							}
							gotdimensions=1;
							fixedlimits=1;
						}
						if (strcmp(word,(char*)GRIDSIZE_KEYWORD)==0)
						{	if(option.verbose|option.debug)
							{	INDENT(1);printf("%s:",GRIDSIZE_KEYWORD);
							}
							GETWORD(word,inp);//get no. rings
							sscanf(word,"%ld",&nx[0]);
							if(option.verbose|option.debug)printf(" %ld rings, ",nx[0]);FLUSH;
							GETWORD(word,inp);//get no. slices
							sscanf(word,"%ld",&nx[1]);
							if(option.verbose|option.debug)printf("%ld slices\n",nx[1]);FLUSH;
							cyl.setgrid(nx);
							gotgridsize=1;
						}
						if (strcmp(word,BOUNDARY_CONDITIONS_TAG)==0)
						{	if(option.verbose|option.debug)
							{	INDENT(1);printf("%s (wall, bottom, top):",BOUNDARY_CONDITIONS_TAG);
							}
							for (int i=0;i<3;i++)
							{	GETWORD(word,inp);
								sscanf(word,"%d",bc+i);
								if(option.verbose|option.debug)
									printf(" %d",bc[i]);
								bc[i]=-bc[i];//stored as negative numbers
								cyl.setbc(bc);
							}
							if(option.verbose|option.debug)printf("\n");
							gotboundary=1;
						}
					}
					if 
					(    !gotdimensions
						|| !gotgridsize
						|| !gotboundary
					)	
					{	fprintf
						(	stderr,
							"One of the keywords not found in the description of domain %s: %s=%d or %s=%d or %s=%d\n",
							name,
							GRIDDIMENSIONS_TAG,gotdimensions,
							GRIDSIZE_KEYWORD,gotgridsize,
							BOUNDARY_CONDITIONS_TAG,gotboundary
						); exit(1);
					}
					cyl.init(nx,dd);
					Ne[nodes]=cyl.getnodes();
					Ne[edges]=0;
					Ne[cells]=cyl.getcells();
					if(option.verbose)
					{	printf("No.of nodes= %d\n",Ne[nodes]);FLUSH;
						printf("No.of cells= %d\n",Ne[cells]);FLUSH;
					}
					{//Determine the number of boundary faces
						int	
							nrings=nx[0],
							nslices=nx[1],
							nringfaces,nsectorfaces,
							nwallfaces,nwallsectorfaces;
						nsectorfaces=nrings*nrings;
						nringfaces=6*nsectorfaces;
						nwallsectorfaces=2*nrings;
						nwallfaces=6*nwallsectorfaces*nslices;
						nbf=2*nringfaces+nwallfaces;
						if(option.verbose)
						{	printf("No.of boundary faces= %d\n",nbf);FLUSH;
						}
					}
					{	//Determine the number of faces
						int nif=4*Ne[cells]-nbf;//number_of_internal_faces
						if (nif%2==1)ERROR("Grid inconsistency in the number of faces");
						nif/=2;
						Ne[boundary_faces]=nbf;
						Ne[faces]=nif+nbf;
						if(option.verbose)
						{	printf("No.of faces= %d\n",Ne[faces]);FLUSH;
						}
					}
					//Ne[points]=0;
					//	if (option.verbose|option.debug)
					//	{
					//		for (int itype=0; itype<maxelements; itype++)
					//		{	INDENT(1);
					//			printf("Number of %s: %d\n",elementype[itype],Ne[itype]);
					//		}
					//	}
					if(option.verbose){printf("Allocating %d cylinder nodes\n",Ne[nodes]);FLUSH;}
					alloc(Ne[nodes],nodes,dnode_root);
					if(option.verbose){printf("Allocating %d cylinder cells\n",Ne[cells]);FLUSH;}
					alloc(Ne[cells],cells,dcell_root);
					cyl.create(dnode_root,dcell_root);
					setBC();
				}
///}//END read_cylinder
						else
						if(strcmp(word,"tetra")==0)
						{//gridtype=unstructured_tetrahedral_grid
							while 
							(	!feof(inp)&&
								( gotnodes==0
								||gotcells==0
								)
							)
							{	GETWORD(word,inp);
								if (*word=='!') SKIPLINE(inp);
								if (*word=='{') level++;
								if (*word=='}') level--; 
								if (level==0) break;
								if(!gotnodes&&strcmp(word,NODES_TAG)==0)
								{
									int	nn,//number of nodes
										mvar=getVarBufSize(nodes);
									DNode	*current,*last;//initiate the pool
																			 //of dynamic nodes
									GETWORD(word,inp);
									sscanf(word,"%d",&nn);
									if(option.verbose|option.debug)
									{	INDENT(1);printf("Loading %d %s ...",nn,NODES_TAG);
									}
									for (int i=0; i<nn; i++)
									{//-	int	type;
										double	x[DIM];
										if (i>0)
										{	DNode	*old=current;
											current->next=new DNode(mvar);
											current=current->next;
											current->prev=old;
											current->next=dnode_root;
											dnode_root->prev=current;
										}
										else
										{	dnode_root=new DNode(mvar);
											current=dnode_root;
											current->next=current->prev=current;
										}
										for (int j=0; j<DIM; j++)
										{
											GETWORD(word,inp);
											sscanf(word,"%lg",x+j);
											current->x[j]=x[j];
										}
										current->state.boundary=1;
									}
									last=current;
									last->next=dnode_root;
									dnode_root->prev=last;
									gotnodes=1;
									if(option.verbose|option.debug)printf(" done\n");
								}
								else
								if(!gotcells&&strcmp(word,CELLS_TAG)==0)
								{	int	nc,//number of cells
										mvar=getVarBufSize(cells);
									DCell	*current,*last;//initialize the pool  of dynamic cells
									if (dnode_root==NULL)
									{	fprintf(stderr,"CAN'T LOAD CELLS BEFORE NODES ARE DEFINED\n");
										exit(1);
									}
									GETWORD(word,inp);
									sscanf(word,"%d",&nc);
									if(option.verbose|option.debug)
									{	INDENT(1);printf("Loading %d %s ...",nc,CELLS_TAG);FLUSH;
									}
									current=dcell_root;
									for (int i=0; i<nc; i++)
									{	if (i>0)
										{	DCell	*old=current;
											current->next=new DCell(mvar);
											current=current->next;
											current->prev=old;
											current->next=dcell_root;
											dcell_root->prev=current;
											//last_dcell=current;
										}
										else
										{	dcell_root=new DCell(mvar);
											current=dcell_root;
											current->next=current->prev=current;
										}
										for (int j=0; j<Nv; j++)
										{	int	vert;
											DNode	*p=dnode_root;
											GETWORD(word,inp);
											sscanf(word,"%d",&vert);vert--;
											for (int k=0; k<vert; k++)
												p=p->next;
											current->vert[j]=p;
											current->neib[j]=NULL;
										}
										if(option.debug)
											printf
											(	"\n\tvert:%x,%x,%x,%x",
												current->vert[0],
												current->vert[1],
												current->vert[2],
												current->vert[3]
											);
									}
									last=current;
									if (current!=NULL)
									{	current->next=dcell_root;
										dcell_root->prev=last;
									}
									//Compute connectivity
									//...
									//Compute boundary and near-boundary cells
									current=dcell_root;
									nbc=0;//number of boundary cells
									if(option.verbose|option.debug){printf(" done\n");FLUSH;}
									gotcells=1;
								}
							}
						}
						else
						{	fprintf
							(	stderr,
								"WRONG GRID-TYPE: %s\n",word
							);exit(1);
						}
					}
				}
				else
				if (strcmp(word,"user")==0)
				{	char	geofile[MAXLINLEN];
					type=user;
					int
						gotdimensions=0,
						gotfile=0;
//-						gottool=0;
					if(option.verbose|option.debug)
					{	INDENT(0);printf("\tType = %d (user)\n",type);
					}
					do
					{	GETWORD(word,inp);
						if (*word=='!') SKIPLINE(inp);
						if (*word=='{') level++;
						if (*word=='}') level--;
						if (level==0) break;
						if (!gotfile&&strcmp(word,FILENAME_TAG)==0)
						{	GETWORD(word,inp);
							strcpy(geofile,word);
							loadGeom(geofile);
							gotfile=1;
						}
						if(!gotdimensions&&strcmp(word,DIMENSIONS_TAG)==0)
						{	if(option.verbose|option.debug)
							{	INDENT(1);printf("%s",DIMENSIONS_TAG);
							}
							for (int i=0; i<DIM; i++)
							{	GETWORD(word,inp);
								sscanf(word,"%lg",xmin+i);
								GETWORD(word,inp);
								sscanf(word,"%lg",xmax+i);
								if(option.verbose|option.debug)
								{	printf("\t%lg : %lg",xmin[i],xmax[i]);FLUSH;
								}
							}
							if(option.verbose|option.debug)printf("\n");
							gotdimensions=1;
							fixedlimits=1;
						}
//-							if(!gottool&&strcmp(word,TOOL_TAG)==0)
//-							{	GETWORD(word,inp);//tool file name
//-								tool=new Tool(word);
//-								gottool=1;
//-							}
					}	while(!gotfile||!gotdimensions);
					if (gotfile!=1)
						ERROR("Keyword \"file\" not found for domain \"user\"");
					gotword=1;
				}
				else
				{	fprintf(stderr,"\nUNKNOWN DOMAIN TYPE: '%s'\n",word);
					exit(1);
				}
			}
			else
			if (strcmp(word,(char*)OUTPUT_TYPE_KEYWORD)==0)
			{	using namespace Output;
				GETWORD(word,inp);
				SCOPY(output_type_name,word);
				if (strcmp(word,"default")==0)
					outype=default_output;
				else
				if (strcmp(word,"ensight")==0)
					outype=ensight_output;
				else
				if (strcmp(word,"tecplot")==0)
					outype=tecplot_output;
				else
				if (strcmp(word,"user")==0)
					outype=user_output;
				else
				{	WARNING1("Bad output type",word);
					SCOPY(output_type_name,"default");
				}
				if(option.verbose) 
				{	printf("Ouput type is %d: %s\n",(int)outype,output_type_name);FLUSH;
				}
			}
			if (memcmp(word,TOOL_TAG,strlen(TOOL_TAG))==0)
			{	GETWORD(word,inp);//tool file name
				tool=new Tool(word);
			}
			else
			if (memcmp(word,(char*)POINTS_KEYWORD,strlen(POINTS_KEYWORD))==0)
			{//Get the number of points (particles) in the domain
				GETWORD(word,inp);
				sscanf(word,"%d",&mp);
				Ne[points]=mp;
				np=0;
				if (option.verbose)
					printf("\tNumber of points: %d\n",mp);
				alloc(points,origin);
				allocvar(points,coordinates);
				initp();//set active particles to 0
			}
			else
			if (memcmp(word,ROTATE_TAG,strlen(ROTATE_TAG))==0)
			{	char axis=word[strlen(ROTATE_TAG)];
				int	nn=Ne[nodes],
					iaxis,jaxis,kaxis;
				double angle,*X=coordinates[nodes].val;
				GETWORD(word,inp);
				sscanf(word,"%lg ",&angle);
				if (option.verbose)
					printf("\tRotating domain around %c-axis by %g degrees\n",axis,angle);
				angle*=PI/180.0;
				switch (axis)
				{
					case 'x': iaxis=0; break;
					case 'y': iaxis=1; break;
					case 'z': iaxis=2; break;
					default:
						fprintf
						(	stderr,
							"UNRECOGNIZED ROTATION AXIS: %c, SHOULD BE: x, y or z\n",
							axis
						);
						exit(1);
				}
				jaxis=(iaxis+1)%DIM;
				kaxis=(iaxis+2)%DIM;
				if(type==dynamic)
				{	DNode	*node=dnode_root;
					if(dnode_root!=NULL)
					do
					{	double
							tmp,r,r2=0.0,
							*x=node->x;
						tmp     = x[jaxis]*cos(angle)+x[kaxis]*sin(angle);
						x[kaxis]=-x[jaxis]*sin(angle)+x[kaxis]*cos(angle);
						x[jaxis]=tmp;
						node=node->next;
					}	while(node!=dnode_root);
				}
				else
				for (int i=0; i<nn; i++)
				{	double	tmp,
						r,r2=0.0,
						*x=X+DIM*i;
					//VALID FOR DIM=3 ONLY
					tmp     = x[jaxis]*cos(angle)+x[kaxis]*sin(angle);
					x[kaxis]=-x[jaxis]*sin(angle)+x[kaxis]*cos(angle);
					x[jaxis]=tmp;
				}
			}
			else
			if (strcmp(word,TRANSLATE_TAG)==0)
			{	int	nn=Ne[nodes];
				double *X=coordinates[nodes].val,
					dx[DIM];
				for (int i=0; i<DIM; i++)
				{	GETWORD(word,inp);
					sscanf(word,"%lg ",dx+i);
				}
				if(option.verbose|option.debug)
				{
					printf("\tTranslating domain %d by (%g, %g, %g)",idomain,dx[0],dx[1],dx[2]);
				}
				if(type==dynamic)
				{	DNode	*node=dnode_root;
					if(dnode_root!=NULL)
					do
					{	double	*x=node->x;
						for (int j=0; j<DIM; j++)
							x[j]+=dx[j];
						node=node->next;
					}	while(node!=dnode_root);
				}
				else
				for (int i=0; i<nn; i++)
				{	double	*x=X+DIM*i;
					for (int j=0; j<DIM; j++)
						x[j]+=dx[j];
				}
				if(option.verbose|option.debug)printf(".\n");
			}
			//if (ivar<nvar&&level==1&&strcmp(word,VARIABLE_KEYWORD)==0)
			//{	char	*p;
			//	Variable *var=variable+ivar++;
			//	GETWORD(word,inp);
			//	var->setvar(word,&level,type,Ne,inp,inpfilename);
			//}
			//Count number of local variables, procedures and functions 
			// for the procedure
			// Procedures are not callable from a procedure at this time
			//-	if (iproc<nproc&&level==1&&strcmp(word,PROCEDURE_KEYWORD)==0)
			//-	{	int	localvar=0;
			//-		Procedure *proc=procedure+iproc++;
			//-		proc->setzero();
			//-		GETWORD(word,inp);
			//-		while (*word=='!'&&!feof(inp)) {SKIPLINE(inp);GETWORD(word,inp);}
			//-		if (*word!='{') 
			//-		{	SCOPY(proc->name,word);
			//-			//printf(": '%s'\n",proc->name);
			//-			GETWORD(word,inp);
			//-			while (*word=='!'&&!feof(inp)) {SKIPLINE(inp);GETWORD(word,inp);}
			//-			if (*word!='{')
			//-			{	fprintf(stderr,"%s: %d: %s",PROCEDURE_KEYWORD,iproc,proc->name);
			//-				ERROR("Sybol '{' missing after "PROCEDURE_KEYWORD);
			//-			}
			//-		}
			//-		//else printf("\n");
			//-		level++;
			//-		//Scan procedure description
			//-		localvar=0;
			//-		while (!feof(inp)&&level>=2)
			//-		{	GETFWORD(word,inp);
			//-			if (*word=='!') {SKIPLINE(inp);continue;}
			//-			if (*word=='{') level++;
			//-			if (*word=='}') level--;
			//-			if (level<2)break;
			//-			if (strcmp(word,PROCEDURE_TYPE_KEYWORD)==0)
			//-			{	GETWORD(word,inp);
			//-				if (isdigit(*word))
			//-					sscanf(word,"%d",&proc->type);
			//-				else
			//-				{	for (int itype=0; itype<DIM+1; itype++)
			//-						if (strcmp(word,elementype[itype])==0) proc->type=itype;
			//-				}
			//-				continue;
			//-			}
			//-			if (strcmp(word,PROCEDURE_NUMBER_KEYWORD)==0)
			//-			{	GETWORD(word,inp);
			//-				sscanf(word,"%d",&proc->niter);
			//-				continue;
			//-			}
			//-			if (level==2)
			//-			{	if (strcmp(word,VARIABLE_KEYWORD)==0) proc->nvar++;
			//-				if (strchr(word,'\(')!=NULL) proc->nfunc++;
			//-			}
			//-		}
			//-		if(option.verbose|option.debug)
			//-		{	INDENT(1);
			//-			printf
			//-			(	"Allocating %d*%d=%d bytes for variables\n",
			//-				proc->nvar,sizeof(struct Variable),
			//-				proc->nvar*sizeof(struct Variable)
			//-			);
			//-		}
			//-		proc->var=new Variable[proc->nvar];
			//-		if (proc->nfunc==0)
			//-		{	FLUSH;
			//-			fprintf
			//-			(	stderr,
			//-				"WARNING: No functions specified for %s '%s'\n",
			//-				PROCEDURE_KEYWORD,proc->name
			//-			);fflush(stderr);
			//-		}
			//-		else
			//-		{
			//-			if(option.verbose|option.debug)
			//-			{	INDENT(1);
			//-				printf
			//-				(	"Allocating %d*%d=%d bytes for functions\n",
			//-					proc->nfunc,sizeof(struct Function),
			//-					proc->nfunc*sizeof(struct Function)
			//-				);
			//-			}
			//-			proc->function=new Function[proc->nfunc];
			//-		}
			//-	}
			//-	if (strchr(word,'\(')!=NULL) nfunc++;
		}//if level==1
	}//while(!feof)
///	if(type_specified==0) ERROR("Domain type not specified");
	if (level!=0)
		ERROR("Unmatched braces");
	if (model==-1)
	{	fprintf
		(	stderr,
			"No model specified for domain %d in %s\n",
			idomain+1,inpfilename
		);	exit(1);
	}
	fclose(inp);
	if(option.debug){printf("Model=%d\nFile %s closed\n",model,inpfilename);FLUSH;}
//-		if (nfunc>0)//Aloocate domain functions
//-		{	function=new Function[nfunc];
//-			if(option.verbose|option.debug)
//-			{	INDENT(0);
//-				printf
//-				(	"Allocating %d*%d=%d bytes for level=%d functions\n",
//-					nfunc,sizeof(struct Function),
//-					nfunc*sizeof(struct Function),
//-					level
//-				);
//-			}
//-		}
	//PASS 2: Seting up procedures
	//Locate the domain
//-	 OPENREAD(inpfilename,inp);
//-	 printf("File %s: PASS 2: Setting up procedures.\n",inpfilename);
//-	 level=0; idomain=ndomain;
//-	 printf("Domain %d: '%s'\n",ndomain+1,name);FLUSH;
//-	 while (!feof(inp)&&idomain>=0&&level>=0)
//-	 {	GETWORD(word,inp);
//-		 if (*word=='!') {SKIPLINE(inp);continue;}
//-		 if (*word=='{') level++;
//-		 if (*word=='}') level--;
//-		 if (level==0&&strcmp(word,DOMAIN_KEYWORD)==0)idomain--;
//-	 }
//-	 GETWORD(word,inp);//Get domain name
//-	 while (!feof(inp)&&*word!='{')
//-		 GETWORD(word,inp);
//-	 level++;
//-	 //Scan domain description
//-	 ivar=0;iproc=0;ifunc=0;
//-	 operation=NULL;
//-	 while (!feof(inp)&&level>=1)
//-	 {	struct Operation	*current;
//-		 GETWORD(word,inp);
//-		 if (*word=='!') {SKIPLINE(inp);continue;}
//-		 if (*word=='{') level++;
//-		 if (*word=='}') level--; 
//-		 if (level==0) break;
//-		 //Scan procedure description
//-		 if (iproc<nproc&&level==1&&strcmp(word,PROCEDURE_KEYWORD)==0)
//-		 {	//Procedure header
//-			 int localvar,localfunc;
//-			 Procedure *proc=procedure+iproc++;
//-			 if(option.verbose|option.debug)
//-			 {
//-				 INDENT(0);printf("%s %d",PROCEDURE_KEYWORD,iproc);FLUSH;
//-				 if (proc->name!=NULL) printf(": '%s'\n",proc->name);
//-				 else printf("\n");
//-				 INDENT(1);printf("%s=%s\n",PROCEDURE_TYPE_KEYWORD,elementype[proc->type]);
//-				 INDENT(1);printf("%s=%d\n",PROCEDURE_NUMBER_KEYWORD,proc->niter);
//-				 INDENT(1);printf("local variables = %d\n",proc->nvar);
//-			 }
//-			 GETWORD(word,inp);
//-			 while (*word=='!'&&!feof(inp)) {SKIPLINE(inp);GETWORD(word,inp);}
//-			 if (*word!='{') 
//-			 {	GETWORD(word,inp);
//-				 while (*word=='!'&&!feof(inp)) {SKIPLINE(inp);GETWORD(word,inp);}
//-				 if (*word!='{') 
//-					 ERROR("Symbol '{' missing after "PROCEDURE_KEYWORD);
//-			 }
//-			 level++;
//-			 //Procedure body
//-			 localvar=0;localfunc=0;
//-			 while (!feof(inp)&&level>1)
//-			 {	int	rank;//,ifoundvar[3];
//-				 char	*s,*q,*arg[3];
//-				 Function	*func;
//-				 //Variable	*var[3];
//-				 //	struct	Function
//-				 //	{	int	rank,type,narg,*ind;
//-				 //		char	name[MAXLINLEN];
//-				 //	}	f;
//-//				for (int i=0; i<3; i++)ifoundvar[i]=0;
//-				 GETWORD(word,inp);
//-				 if (*word=='!') {SKIPLINE(inp);continue;}
//-				 if (*word=='{') {level++;continue;}
//-				 if (*word=='}') {if(--level<2)break;else continue;}
//-				 if 
//-				 (	  strcmp(word,PROCEDURE_TYPE_KEYWORD)==0
//-					 ||strcmp(word,PROCEDURE_NUMBER_KEYWORD)==0
//-				 )
//-				 {	//TYPE: skip - was read in PASS 1
//-					 GETWORD(word,inp);
//-					 continue;
//-				 }
//-				 //Local variable inside procedure
//-				 if (level==2&&strcmp(word,VARIABLE_KEYWORD)==0)
//-				 {	char	*p;
//-					 Variable *v;
//-					 if (localvar==proc->nvar) ERROR("Local varibles mismatch");
//-					 v=proc->var+localvar++;
//-					 GETWORD(word,inp);
//-					 //v->rank=getrank(word);
//-					 //setind(word,v->rank,&ind[0]);
//-					 v->setvar(word,&level,proc->type,inp,inpfilename);
//-					 continue;
//-				 }
//-				 //Function description
//-				 if (localfunc>=proc->nfunc) 
//-				 {	fprintf
//-					 (	stderr,
//-						 "WARNING: The total number of arguments is too big for the total number of functions: %d\n",
//-						 proc->nfunc
//-					 );
//-					 break;
//-				 }
//-				 func=proc->function+localfunc++;
//-				 //get the first function argument
//-				 setrank(word,&rank,&func->vind[0]);
//-//			findvar(word,rank,localvar,&ifoundvar[0],proc->var,variable,&(func->var[0]));
//-				 if 
//-				 (	   findvar(word,rank,localvar,proc->var,&func->var[0])==0
//-					 && findvar(word,rank,nvar,variable,&func->var[0])==0
//-				 )	ERROR1("Word not found:",word);
//-				 //Function name
//-				 GETWORD(word,inp);
//-				 setrank(word,&func->rank,&func->ind);
//-				 //get function type
//-				 if ((p=strchr(word,'_'))==NULL)
//-					 func->type==0; //default type
//-				 else
//-				 {	int	i; *p++='\0';
//-					 for (i=0; i<maxfunctypes; i++) 
//-						 if (strcmp(p,functype[i])==0)
//-						 {	func->type=i; break; }
//-					 if (i==maxfunctypes)
//-						 ERROR1("Unknown argument type",p);
//-				 }
//-				 SCOPY(func->name,word);
//-				 //get function arguments 
//-				 GETWORD(word,inp);
//-				 //count function arguments
//-				 func->narg=2;
//-				 for (p=word;*p!='\0';p++) 
//-					 if (*p==',') func->narg++;
//-					 else
//-					 if (*p=='[')
//-						 while (*p!='\0'&&*p!=']')p++;
//-				 if (func->narg<2) ERROR1("Too few arguments in function",func->name);
//-				 if (func->narg>3) ERROR1("Too many arguments in function",func->name);
//-				 if (option.verbose|option.debug)
//-				 {
//-					 INDENT(0);printf("Function %s\n",func->name);FLUSH;
//-					 INDENT(1);printf("type= %s\n",functype[func->type]);FLUSH;
//-					 INDENT(1);printf("narg= %d\n",func->narg);FLUSH;
//-					 INDENT(1);printf("rank= %d, ind=",func->rank);FLUSH;
//-					 for (int i=0; i<func->rank; i++) {printf(" %c",func->ind[i]);FLUSH;}
//-					 printf("\n");
//-					 INDENT(1);
//-					 printf
//-					 (	"target= %s, type=%d, rank=%d, ind=",
//-						 func->var[0]->name,func->var[0]->type,func->var[0]->rank
//-					 );
//-					 for (int i=0; i<func->var[0]->rank; i++) printf(" %c",func->vind[0][i]);
//-					 printf("\n");FLUSH;
//-				 }
//-				 //assign argument's names to arg[i]
//-				 p=word;//beginning of the arguments
//-				 for (int i=1; i<func->narg; i++)
//-				 for (arg[i]=p;*p!='\0';p++)
//-					 if (*p==',') 
//-					 {	*p++='\0'; break;}
//-					 else
//-					 if (*p=='[')
//-					 {	//*p++='\0';
//-						 while (*p!='\0'&&*p!=']')p++;
//-					 }
//-				 //Determine the source arguments
//-				 for (int i=1; i<func->narg; i++)
//-				 {
//-					 setrank(arg[i],&rank,&func->vind[i]);
//-//					findvar(arg[i],rank,localvar,&ifoundvar[i],proc->var,variable,&func->var[i]);
//-//					findvar(arg[i],rank,localvar,proc->var,&func->var[i]);
//-					 if 
//-					 (	   findvar(arg[i],rank,localvar,proc->var,&func->var[i])==0
//-						 && findvar(arg[i],rank,nvar,variable,&func->var[i])==0
//-					 )	ERROR1("Word not found:",word);
//-					 if(option.verbose|option.debug)
//-					 {	INDENT(1);
//-						 printf
//-						 (	"source[%d]='%s', rank=%d, ind=",
//-							 i,func->var[i]->name,func->var[i]->rank
//-						 );
//-						 for (int j=0; j<func->var[i]->rank; j++) 
//-							 printf(" %c",func->vind[i][j]);printf("\n");
//-					 }
//-				 }
//-				 if //check for recursion
//-				 (	func->narg==2
//-					 &&strcmp(func->var[0]->name,func->var[1]->name)==0
//-				 )	ERROR1("Recursive invocation not allowed for function",func->name);
//-				 func->setind(func->var);//Setting up function indexes
//-				 func->setfun();//Setting up  function pointers
//-			 }
//-			 if (localfunc<proc->nfunc)
//-				 ERROR("Some function definitions were incomplete");
//-		 }
//-		 if (level==1&&strcmp(word,DOMAIN_OPERATION_KEYWORD)==0)
//-		 {
//-			 if (operation==NULL)
//-			 {	operation=new Operation;
//-			 }
//-		 }
//-		 if (level==1&&strchr(word,'\(')!=NULL)
//-		 {	struct Function	*func;
//-			 if (ifunc>=nfunc) ERROR1("Function numbers don't agree at",word);
//-			 func=function+ifunc++;
//-		 }
//-	 }
//-	 fclose(inp);
}
//	int	Domain::findvar
//	(
//		char	*word,
//		int	rank,
//		int	nvar,//localvar,
//	//	int	*ifoundvar,
//	//	Variable *procvar,
//		Variable	*variable,
//		Variable	**var
//	)
//	{	int	ivar;
//		//determine if it's an earlier defined local variable
//	//		if (!*ifoundvar)
//	//		{	for (int ivar=0; ivar<localvar; ivar++)
//	//			if (strcmp(word,procvar[ivar].name)==0)
//	//			{	if (rank!=procvar[ivar].rank)
//	//					ERROR1("Declaration/Invocation ranks do not match for",word);
//	//				*var=procvar+ivar;*ifoundvar=1;
//	//				break;
//	//			}
//	//		}
//		//determine if it's a global variable 
//	//	if (!*ifoundvar)
//	//	{
//			for (ivar=0; ivar<nvar; ivar++)
//			if (strcmp(word,variable[ivar].name)==0)
//			{	if (rank!=variable[ivar].rank)
//					ERROR1("Declaration/Invocation ranks do not match for",word);
//				*var=variable+ivar;//ifoundvar=1;
//				break;
//			}
//	//	}
//	//	if (ifoundvar==0)
//		if (ivar==nvar) return 0;
//		else return 1;
//	}
void	Domain::setGeomXML
(
	char inpfilename[]
)
{	
	using namespace Input;
	const char *doctype=DOCTYPE;
	int	
		gotword=0,
		idomain=this-domain_root;//,ivar;
	char	*p,word[MAXLINLEN];
	FILE	*inp;
	if(option.verbose)printf("Domain::setGeomXML:inpfilename=%s\n",inpfilename);
	xmlDocPtr doc;
	xmlNodePtr root;
	//Corrdinates of basic element types
	coordinates=new Variable[maxelements];
	volumes=new Variable[maxelements];
	for (int itype=0; itype<maxelements; itype++)
	{	coordinates[itype].type=itype;
		coordinates[itype].dimension=DIM;
		volumes[itype].type=itype;
		volumes[itype].dimension=1;
		Ne[itype]=0;
	}
	if (option.verbose|option.debug)
	{	printf
		(	"Reading geometry for domain no. %d: from %s\n",
			idomain+1,inpfilename
		);FLUSH;
	}
	doc = xmlParseFile(inpfilename);
	if (doc == NULL ) 
	{	fprintf(stderr,"XML parser failed in %s\n",inpfilename);
		return;
	}
	root = xmlDocGetRootElement(doc);
	if (root == NULL) 
	{	fprintf(stderr,"Empty document: %s\n",inpfilename);
		xmlFreeDoc(doc);
		return;
	}
	if (xmlStrcmp(root->name, (const xmlChar *) doctype)) 
	{	fprintf
		(	stderr,
			"document of the wrong type, root node != %s in %s\n",
			doctype,inpfilename
		);
		xmlFreeDoc(doc);
		return;
	}
	for 
	(	xmlNodePtr cur = root->xmlChildrenNode; 
		cur != NULL; cur = cur->next
	) 
	{	//SCAN DESCRIPTION AND FIND THIS DOMAIN
		if ((!xmlStrcmp(cur->name, (const xmlChar *)DOMAIN_KEYWORD)))
		{	if(idomain>0){idomain--;continue;}	
			char modelid[MAXLINLEN];
			if(getCharAttr(cur,(char*)DOMAIN_NAME_TAG,name)==0) 
			{	fprintf
				(	stderr,"CAN'T GET %s FOR %s IN %s\nAborting\n",
					DOMAIN_NAME_TAG,DOMAIN_KEYWORD,inpfilename
				); exit(1);
			}
			if(option.verbose) printf("Domain: %s\n",name);
//?			idomain=this-domain_root;
			if(getCharAttr(cur,(char*)MODEL_KEYWORD,modelid)==0) exit(1);
			if(option.verbose) printf("Model ID: %s\n",modelid);
			if
			(	(	model=getIntAttr
					(	doc,(char*)MODEL_KEYWORD,(char*)MODEL_ID_TAG,modelid,(char*)"type")
				)==-1
			) 	exit(1);
			if(option.verbose){ printf("Model no.: %d\n",model);fflush(stdout);}
#ifdef WITH_DOVE
			if((iproc=parseInt(doc,cur,DOMAIN_PROC_TAG))<0) exit(1);
			if(option.verbose){ printf("Domain processor: %d\n",iproc);fflush(stdout);}
#endif
			if(parseWord(doc,cur,(char*)DOMAIN_TYPE_KEYWORD,word)==0) exit(1);
			if(option.verbose){ printf("Domain type: %s\n",word);fflush(stdout);}
			if (strcmp(word,"dynamic")==0)
			{	char	geofile[MAXLINLEN];
				enum GridType	gridtype=unstructured_tetrahedral_grid;
				type=dynamic;
				//first_framenode=NULL;
				if(option.verbose|option.debug)
				{	INDENT(0);printf("\tType dynamic is valid\n");
				}
				for 
				(	xmlNodePtr next=cur->xmlChildrenNode;
					next != NULL;
					next=next->next
				) 
				{	if ((!xmlStrcmp(next->name, (const xmlChar *)GEOMETRY_TAG)))
					{	if(parseWord(doc,next,(char*)BOUNDS_TAG,word)!=0) 
						{	sscanf
							(	word,"%lg %lg %lg %lg %lg %lg",
								xmin,xmax,xmin+1,xmax+1,xmin+2,xmax+2
							);
							fixedlimits=1;
							if(option.verbose)
							printf
							(	"\tbounds: %g:%g x %g:%g x %g:%g\n",
								xmin[0],xmax[0],xmin[1],xmax[1],xmin[2],xmax[2]
							);
						}
						for 
						(	xmlNodePtr geo = next->xmlChildrenNode;
							geo != NULL; geo=geo->next
						) 
						{	
							if(!xmlStrcmp(geo->name, (const xmlChar *)SHAPE_TAG))
							{	if(option.verbose)
								{	printf("Geometry type: %s\n",SHAPE_TAG);fflush(stdout);
								}
								if(parseWord(doc,geo,(char*)TYPE_TAG,word)==0) exit(1);
								if(strcmp(word,"box")==0)///gridtype=box;
								{	Box	b; //BOX
									if(option.verbose)printf("\tGrid type=box\n");
									if(parseWord(doc,geo,(char*)DIMENSIONS_TAG,word)==0) exit(1);
									sscanf
									(	word,"%d %d %d",
										b.nx,b.nx+1,b.nx+2
									);
									if(option.verbose)
									printf
									(	"\tDimensions: %d x %d x %d\n",
										b.nx[0],b.nx[1],b.nx[2]
									);
									if(parseWord(doc,geo,(char*)BOUNDARY_CONDITIONS_TAG,word)==0) exit(1);
									sscanf
									(	word,"%d %d %d %d %d %d",
										b.bc,b.bc+1,b.bc+2,b.bc+3,b.bc+4,b.bc+5
									);
									for(int i=0;i<6;i++)b.bc[i]=-b.bc[i];//stored as negative numbers
									if(option.verbose)
									printf
									(	"\tBoundary: %d:%d; %d:%d %d:%d\n",
										b.bc[0],b.bc[2],b.bc[3],b.bc[4],b.bc[5],b.bc[6]
									);
									Ne[nodes]=1;
									for (int i=0; i<DIM; i++) Ne[nodes]*=b.nx[i]+1;
									Ne[edges]=0;
									Ne[cells]=1;
									for (int i=0; i<DIM; i++) Ne[cells]*=b.nx[i];
									Ne[cells]*=5;
									nbv=2*((b.nx[0]+1)*(b.nx[1]+1)+(b.nx[0]+b.nx[1])*(b.nx[2]-1));
									Ne[boundary_nodes]=nbv;
									nbf=0;
									for (int i=0; i<DIM; i++) nbf+=4*b.nx[(i+1)%DIM]*b.nx[(i+2)%DIM];
									{	//Determine the number of faces
										int nif=4*Ne[cells]-nbf;//number_of_internal_faces
										if (nif%2==1)ERROR("Grid inconsistency in the number of faces");
										nif/=2;
										Ne[boundary_faces]=nbf;
										Ne[faces]=nif+nbf;
									}
									if (option.verbose|option.debug)
									{
										for (int itype=0; itype<maxelements; itype++)
										{	INDENT(1);
											printf("Number of %s: %d\n",elementype[itype],Ne[itype]);
										}
									}
									if(option.verbose)
									{printf("Allocating %d brick nodes\n",Ne[nodes]);FLUSH;}
									alloc(Ne[nodes],nodes,dnode_root);
									if(option.verbose)
									{printf("Allocating %d brick cells\n",Ne[cells]);FLUSH;}
									alloc(Ne[cells],cells,dcell_root);
									b.create(xmin,xmax,dnode_root);
									b.connect(dnode_root,dcell_root);
								}//END BOX
								else
								if(strcmp(word,"cylinder")==0)///gridtype=cylinder;
								{	int	//CYLINDER
										Ne[maxelements],
										nx[2],//number of rings (nx[0]), and slices (nx[1])
										bc[3];//boundary conditions on cylinder walls:
										      //bc[0]: wall, bc[1]: bottom, bc[2]: top
									double	dd[2];//dd[0]: cylinder radius; dd[1]: cylinder length
									Cylinder	cyl;
									if (option.verbose|option.debug)
									{	INDENT(0);printf("\t%s=cylinder\n",(char*)GRIDTYPE_KEYWORD);
									}
									if(parseWord(doc,geo,(char*)DIMENSIONS_TAG,word)==0) exit(1);
									sscanf
									(	word,"%ld %ld",
										nx,nx+1
									);
									cyl.setgrid(nx);
									if(option.verbose|option.debug)
									{	printf("\t%d rings, %d slices\n",nx[0],nx[1]);FLUSH; 
									}
									if(parseWord(doc,geo,(char*)BOUNDS_TAG,word)==0) exit(1);
									sscanf
									(	word,"%lg %lg",
										dd,dd+1
									);
									if(option.verbose)
									{	printf("\tlenth=%lg, radius=%lg\n",dd[0],dd[1]);FLUSH;
									}
									if(parseWord(doc,geo,(char*)BOUNDARY_CONDITIONS_TAG,word)==0) exit(1);
									sscanf
									(	word,"%d %d %d",
										bc,bc+1,bc+2
									);
									for(int i=0;i<3;i++)bc[i]=-bc[i];//stored as negative numbers
									cyl.setbc(bc);
									cyl.init(nx,dd);
									Ne[nodes]=cyl.getnodes();
									Ne[edges]=0;
									Ne[cells]=cyl.getcells();
									if(option.verbose)
									{	printf("No.of nodes= %d\n",Ne[nodes]);FLUSH;
										printf("No.of cells= %d\n",Ne[cells]);FLUSH;
									}
									{//Determine the number of boundary faces
										int	
											nrings=nx[0],
											nslices=nx[1],
											nringfaces,nsectorfaces,
											nwallfaces,nwallsectorfaces;
										nsectorfaces=nrings*nrings;
										nringfaces=6*nsectorfaces;
										nwallsectorfaces=2*nrings;
										nwallfaces=6*nwallsectorfaces*nslices;
										nbf=2*nringfaces+nwallfaces;
										if(option.verbose)
										{	printf("No.of boundary faces= %d\n",nbf);FLUSH;
										}
									}
									{	//Determine the number of faces
										int nif=4*Ne[cells]-nbf;//number_of_internal_faces
										if (nif%2==1)ERROR("Grid inconsistency in the number of faces");
										nif/=2;
										Ne[boundary_faces]=nbf;
										Ne[faces]=nif+nbf;
										if(option.verbose)
										{	printf("No.of faces= %d\n",Ne[faces]);FLUSH;
										}
									}
									if(option.verbose){printf("Allocating %d cylinder nodes\n",Ne[nodes]);FLUSH;}
									alloc(Ne[nodes],nodes,dnode_root);
									if(option.verbose){printf("Allocating %d cylinder cells\n",Ne[cells]);FLUSH;}
									alloc(Ne[cells],cells,dcell_root);
									cyl.create(dnode_root,dcell_root);
									setBC();
								}//END CYLINDER
							}
							else
							if(!xmlStrcmp(geo->name, (const xmlChar *)GRID_TAG))
							{	if(option.verbose)
								{	printf("Geometry type: %s\n",GRID_TAG);fflush(stdout);
								}
								printf("Structured grids are not implemented yet.\nAborting\n");
								exit(1);
							}
							else
							if(!xmlStrcmp(geo->name, (const xmlChar *)MESH_TAG))
							{//gridtype=unstructured_tetrahedral_grid
								const char *mesh_type[]=
								{	"node_cell",
									"data_file"
								};
								int itype=-1,ntypes=sizeof(mesh_type)/sizeof(mesh_type[0]);
								if((parseWord(doc,geo,(char*)TYPE_TAG,word))==0) 
								{	fprintf(stderr,"No type specified for the mesh.\nAborting\n");
									exit(1);
								}
								for(itype=0;itype<ntypes;itype++)
									if(strcmp(mesh_type[itype],word)==0) break;
								if(itype>=ntypes)
								{	fprintf(stderr,"Wrong mesh type '%s'\n",word);
									fprintf(stderr,"Valid types:\n");
									for(itype=0;itype<ntypes;itype++)
										fprintf(stderr,"\t%s\n",mesh_type[itype]);
									fprintf(stderr,"Aborting\n");
									exit(1);
								}
								if(option.verbose|option.debug)
								{	INDENT(1);printf("Mesh type: %s",mesh_type[itype]);
								}
								switch(itype)
								{
								case 0: //node_cell mesh
								{
								int nn=0,nc=0; //number of nodes and cells
								for 
								(	xmlNodePtr mesh = geo->xmlChildrenNode;
									mesh != NULL; mesh=mesh->next
								) 
								if(!xmlStrcmp(mesh->name, (const xmlChar *)NODES_TAG))
								{	//LOAD NODES
									//	mvar=getVarBufSize(nodes);
									if((nn=parseInt(doc,mesh,(char*)"number"))<=0) 
									{	fprintf(stderr,"No nodes found for the mesh.\nAborting\n");
										exit(1);
									}
									if(option.verbose|option.debug)
									{	INDENT(1);printf("Loading %d %s ...",nn,NODES_TAG);
									}
									int ncoords=0;
									for 
									(	xmlNodePtr coords = mesh->xmlChildrenNode;
										coords != NULL;coords = coords->next
									) 
									{	if 
										(	(!xmlStrcmp
											(	coords->name, (const xmlChar *)COORDINATES_TAG)
											)
										) 
										{	xmlChar *data=xmlNodeListGetString
											(doc,coords->xmlChildrenNode,1);
											ncoords=parseStringNodes(nn,(char*)data);
											xmlFree(data);
										}
									}
									if(nn<=0||ncoords!=nn)
									{	fprintf
										(	stderr,
											"Mismatch in the number of nodes %d and coordinates: %d\n",
											nn,ncoords
										);
										exit(1);
									}
									///loadList(doc,mesh->xmlChildrenNode,COORDINATES_TAG);
								}
								else
								if ((!xmlStrcmp(mesh->name, (const xmlChar *)CELLS_TAG)))
								{	if(nn<=0) 
									{	fprintf(stderr,"Number of nodes: %d is too small\n",nn);
										continue;
									}
									if (dnode_root==NULL)
									{	fprintf(stderr,"CAN'T LOAD CELLS BEFORE NODES ARE DEFINED\n");
										exit(1);
									}
									nc=parseInt(doc,mesh,(char*)"number");
									if(option.verbose|option.debug)
									{	INDENT(1);printf("Loading %d %s\n",nc,CELLS_TAG);
									}
									if(nc<=0) 
									{
									//	exit(1);
										DCell	*current,*last;
										current=dcell_root;
										last=current;
										if (current!=NULL)
										{	current->next=dcell_root;
											dcell_root->prev=last;
										}
										//Compute connectivity
										//...
										//Compute boundary and near-boundary cells
										current=dcell_root;
										nbc=0;//number of boundary cells
										if(option.verbose|option.debug){printf(" done\n");FLUSH;}
									}
									else
									{
									int nnodes=0,ncells=0;
									for 
									(	xmlNodePtr submesh = mesh->xmlChildrenNode;
										submesh != NULL; submesh = submesh->next
									) 
									{	
										if 
										(	(!xmlStrcmp
											(	submesh->name, (const xmlChar *)NODES_TAG)
											)
										) 
										{	xmlChar *data=xmlNodeListGetString
											(doc,submesh->xmlChildrenNode,1);
											ncells=parseStringCells(nn,nc,(char*)data);
											if(ncells!=nc)
											{	fprintf
												(	stderr,
													"Mismatch in the number of cells declared, %d, and actually found: %d\n",
													nc,ncells
												);
												exit(1);
											}
											xmlFree(data);
										}
									}//END for
									}//END if
								}//END if
								if(nn<=0)
								{	fprintf(stderr,"Empty data set in domain %s\n",name);
								}
								}//END case node_cell
								break;
								case 1: //data_file
									if((parseWord(doc,geo,(char*)FILENAME_TAG,word))==0) 
									{	fprintf(stderr,"No filename specified for the mesh.\nAborting\n");
										exit(1);
									}
									strcpy(geofile,word);
									loadGeom(geofile);
									break; //END data_file
								}//END switch itype
							}
							else
							if((!xmlStrcmp(geo->name, (const xmlChar *)POINTS_TAG)))
							{	
								if((mp=parseInt(doc,geo,(char*)NUMBER_TAG))==-1) 
								{	printf("Cant load points: mp=%d\n",mp);
									exit(1);
								}
								Ne[points]=mp;
								np=0;
								if (option.verbose) 
									printf("\tNumber of %s: %d\n",POINTS_TAG,mp);
								alloc(points,origin);
								allocvar(points,coordinates);
								initp();//set active particles to 0
							}//END POINTS_TAG	
							else
							if((!xmlStrcmp(geo->name, (const xmlChar *)TRANSFORM_TAG)))
							{	
								int	nn=Ne[nodes];
								double angle[DIM],
									*X=coordinates[nodes].val;
								for (int i=0;i<DIM;i++) angle[i]=0.0;
								parseFloat(doc,geo,(char*)ROTATEX_TAG,angle[0]);
								parseFloat(doc,geo,(char*)ROTATEY_TAG,angle[1]);
								parseFloat(doc,geo,(char*)ROTATEZ_TAG,angle[2]);
								for (int iaxis=0;iaxis<DIM;iaxis++)
								{	int 
										jaxis=(iaxis+1)%DIM,
										kaxis=(iaxis+2)%DIM;
									double a=PI*angle[iaxis]/180.0;
									if(type==dynamic)
									{	DNode	*node=dnode_root;
										if(dnode_root!=NULL)
										do
										{	double
												tmp,
												*x=node->x;
											tmp     = x[jaxis]*cos(a)+x[kaxis]*sin(a);
											x[kaxis]=-x[jaxis]*sin(a)+x[kaxis]*cos(a);
											x[jaxis]=tmp;
											node=node->next;
										}	while(node!=dnode_root);
									}
									else
									for (int i=0; i<nn; i++)
									{	double	tmp,
											*x=X+DIM*i;
										//VALID FOR DIM=3 ONLY
										tmp     = x[jaxis]*cos(a)+x[kaxis]*sin(a);
										x[kaxis]=-x[jaxis]*sin(a)+x[kaxis]*cos(a);
										x[jaxis]=tmp;
									}
								}
								if(option.verbose)printf("\trotate: %g %g %g\n",angle[0],angle[1],angle[2]);
								if(parseWord(doc,geo,(char*)TRANSLATE_TAG,word)!=0)
								{	double dx[DIM];
									for (int i=0;i<DIM;i++) dx[i]=0.0;
									sscanf(word,"%lg %lg %lg",dx,dx+1,dx+2);
									if(option.verbose)
									{	printf
										(	"\t%s domain %d by %g %g %g\n",TRANSLATE_TAG,name,dx[0],dx[1],dx[2]
										);fflush(stdout);
									}
									if(type==dynamic)
									{	DNode	*node=dnode_root;
										if(dnode_root!=NULL)
										do
										{	double	*x=node->x;
											for (int j=0; j<DIM; j++)
												x[j]+=dx[j];
											node=node->next;
										}	while(node!=dnode_root);
									}
									else
									for (int i=0; i<nn; i++)
									{	double	*x=X+DIM*i;
										for (int j=0; j<DIM; j++)
											x[j]+=dx[j];
									}
								}
							}//END TRANSFORM_TAG	
						}//END geo
					}//END GEOMETRY_TAG
					else
					if ((!xmlStrcmp(next->name, (const xmlChar *)OUTPUT_TAG)))
					{	using namespace Output;
						if(parseWord(doc,next,(char*)TYPE_TAG,word)==0) 
						{	fprintf
							(	stderr,"Tag %s not found in %s section\n",
								TYPE_TAG,OUTPUT_TAG
							);exit(1);
						}
						SCOPY(output_type_name,word);
						if (strcmp(word,"default")==0)
							outype=default_output;
						else
						if (strcmp(word,"ensight")==0)
							outype=ensight_output;
						else
						if (strcmp(word,"tecplot")==0)
							outype=tecplot_output;
						else
						if (strcmp(word,"user")==0)
							outype=user_output;
						else
						{	WARNING1("Bad output type",word);
							SCOPY(output_type_name,"default");
						}
						if(option.verbose) 
						{	printf("Ouput type is %d: %s\n",(int)outype,output_type_name);FLUSH;
						}
					}//END OUTPUT_TAG
					else
					if ((!xmlStrcmp(next->name, (const xmlChar *)TOOL_TAG)))
					{	char filename[MAXLINLEN];
						if(parseWord(doc,next,(char*)FILENAME_TAG,filename)==0)
						{	fprintf(stderr,"CAN'T GET %s FOR %s\n",FILENAME_TAG,TOOL_TAG);
							exit(1);
						}
						if(option.verbose)
						{	printf("%s %s: %s\n",TOOL_TAG,FILENAME_TAG,filename);fflush(stdout);
						}
						tool=new Tool(filename);
					}//END TOOL_TAG
				}//END next
			}
			break;
		}
	}
	xmlFreeDoc(doc);
	if(option.debug){printf("Model=%d\nFile %s closed\n",model,inpfilename);FLUSH;}
}
int	Domain::parseStringNodes(int nn, char *data)
{
	int mvar, inode,
		datasize=strlen(data);
	char *inp=data; //current data-pointer
	DNode	*current,*last;
	//Read node data
	mvar=getVarBufSize(nodes);
	if (dnode_root!=NULL)
		deleteRing(dnode_root);
	//Read the number of nodes
//	fread(&nn,sizeof(nn),1,inp);
	if (option.verbose)
		printf("\tNumber of nodes = %d\n",nn);
	//Allocate the node list
	for (inode=0; inode<nn && inp-data<datasize; inode++)
	{	if(inode>0)
		{	DNode	*old=current;
			current->next=new DNode(mvar);
			current=current->next;
			current->prev=old;
			current->next=dnode_root;
			dnode_root->prev=current;
		}
		else
		{	dnode_root=new DNode(mvar);
			current=dnode_root;
			current->next=current->prev=current;
		}
		for (int i=0;i<DIM;i++)
		{
			while( isspace(*inp)&&inp-data<datasize)inp++;
			if(inp-data==datasize)break;
			sscanf(inp,"%lg",current->x+i);
			while(!isspace(*inp)&&inp-data<datasize)inp++;
			if(inp-data==datasize)break;
		}
		current->state.boundary=1;
	}
	last=current;
	last->next=dnode_root;
	dnode_root->prev=last;
	if (option.verbose) printf("%d Node data read\n",inode);
	if (inode<nn)
	{	fprintf(stderr,"ERROR IN %s: Number of nodes specified %d does not match actually found: %d\n",nn,inode);
		exit(1);
	}
	return inode;
}
int	Domain::parseStringCells(int nn, int nc, char *data)
{
	//Read connectivities
	//Read cell data
	char *inp=data;
	int 
		icell,datasize=strlen(data),
		mvar=getVarBufSize(cells);
	DCell	*current_cell,*last_cell,**cell;
	if(nn<=0&&nc<=0) return 0;
	if(dnode_root==NULL)
	{	fprintf
		(	stderr,
			"Domain::parseStringCells:Can't parse cells since node array is empty\n"
		);
		return 0;
	}
	DNode	**nodelist=new DNode*[nn];
	nodelist[0]=dnode_root;
	for (int i=1; i<nn; i++)
	{	nodelist[i]=nodelist[i-1]->next;
		if(nodelist[i]==nodelist[0])
		{	fprintf
			(	stderr,
				"Domain::parseStringCells:Number of nodes %d greater than allocated: %d\n",
				nn,i
			);
			return 0;
		}
	}
	if (dcell_root!=NULL)
		deleteRing(dcell_root);
//	fread(&nc,sizeof(nc),1,inp);
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
///		fread(facetype,sizeof(facetype[0]),Nf,inp);
		//Read cell-node connectivity
		for (int i=0; i<Nv; i++)//TETRAHEDRAL ELEMENTS ONLY
		{	int ivert;
///			fread(&ivert,sizeof(int),1,inp);
			sscanf(inp,"%d",&ivert);
			if (--ivert>=0)
				current_cell->vert[i]=nodelist[ivert];
#ifdef DEBUG
			else
			{	fprintf(stderr,"loadGeom: cell-node connectivity error (ivert=%d)\n",ivert);
				exit(1);
			}
			if (current_cell->vert[i]==NULL)
			{	fprintf
				(	stderr,
					"loadGeom: icell=%d, node[%d(%d)]=%x\n",
					icell,ivert,i,nodelist[ivert]
				);
				exit(1);
			}
#endif
			while(!isspace(*inp)&&inp-data<datasize)inp++;
			if(inp-data==datasize)break;
			while( isspace(*inp)&&inp-data<datasize)inp++;
			if(inp-data==datasize)break;
		}
		current_cell=current_cell->next;
	}
	delete nodelist;
///	//Allocate temporal cells array
///	cell=new DCell*[nc];
///	last_cell=dcell_root;
///	for (int i=0; i<nc; i++)
///	{
///		cell[i]=last_cell;
///		last_cell=last_cell->next;
///	}
///	//Read cell-cell connectivity
///	last_cell=dcell_root;
///	for (int i=0; i<nc; i++)
///	{	for (int i=0; i<Nf; i++)
///		{	int	icell;
///			fread(&icell,sizeof(icell),1,inp);
///			if (icell>=0)
///				last_cell->neib[i]=cell[icell];
///			else
///				last_cell->neib[i]=NULL;
///		}
///		last_cell=last_cell->next;
///	}
///	delete cell;
/////	fclose(inp);
///	if (option.verbose)
///	{	printf("Cell-Node connectivity parsed\n");FLUSH;
///	}
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
		{	fprintf(stderr,"SaveGeom:: grid inconsistency: loose node %d\n",i++);
			exit(1);
		}
		node=node->next;
	}	while(node!=dnode_root);
	}
#endif
	return icell;
}
