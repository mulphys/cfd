#include <iostream>
#include <stdlib.h>
#include <sys/timeb.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <libxml/xmlmemory.h>
#include <libxml/parser.h>
#include "main.h"
#include "io.h"
#include "var.h"
#include "geom.h"
//#include "particles.h"
#include "tool.h"
#include "domain.h"
#include "input.h"
#include "output.h"
#include "gui.h"

#ifdef WITH_MPI
	#include "mpi.h"
	using namespace MPI;
#endif

using namespace std;

struct Option	option={0,0};

char
	programname[MAXLINLEN],
	configfile[MAXLINLEN];

#ifdef WITH_DOVE
/*
int
	iproc=0,//this processor number
	nproc=1;//total number of processors 
*/

int
	iproc,
	nproc;

#endif

uint mask[Nv];//used for cell-cell-node connectivity

struct RunTime	runtime;

using namespace	Gui;

int	ndomains=0;

Domain	*domain_root;

void	readcmdline(int argc, char *argv[])
{	int	i,iarg;
	char
		*gridprmflnam,
		*p,*q;
	for (p=argv[0];!isalpha(*p)&&*p!='\0';p++);
	strcpy(programname,p);
	option.verbose=0;
	option.debug=0;
	iarg=argc;
	for (i=1; i<argc; i++)
	{	if (*argv[i]=='-')
		switch (*(p=argv[i]+1))
		{	case 'f':
				if(i++>=argc)goto Usage;
				strcpy(::configfile,argv[i]);
				break;
			case 'v':
				option.verbose=1;
				break;
			case 'V':
				option.debug=1;
				break;
			case 'h':
				goto Usage;
				break;
		}
	}
	return;
	Usage:
	printf("Usage: %s -<option> [conf_file]\n",argv[0]);
	printf("where <option> can be\n");
	printf("\tv\t- verbose\n");
	printf("\tV\t- more verbose\n");
	printf("\tf file-name\t- specifies name for a configuration file\n\t\t(default: %s)\n",CONFIGFILE);
	exit(1);
}
#ifdef WITH_DOVE
void	connectall()
{//Finding domain overlaps
	for (int i=0; i<ndomains; i++)
	{	Domain	*dom=domain_root+i;
		//cout<<"domain_root["<<i<<"]->iproc = "<<dom->iproc<<endl;

		if(iproc==dom->iproc)//This processor
			dom->storeBoundaryCoordinates();
	}	
	for(int i=0;i<ndomains;i++)
	{	Domain	*dom=domain_root+i;
		if(iproc==dom->iproc)//This processor
			dom->ConnectThisProcDomains(iproc);
	}
#ifdef WITH_MPI
	//Other processors
	{
		for(int jproc=1;jproc<nproc;jproc++)
		{	int	
				to=(iproc+jproc)%nproc,
				from=(iproc+nproc-jproc)%nproc;
			// Send boundary
			for(int idom=0;idom<ndomains;idom++)
			{	Domain	*domi=domain_root+idom;
				if(iproc==domi->iproc)
				for(int jdom=1;jdom<ndomains;jdom++)
				{	Domain *domj=domain_root+(idom+jdom)%ndomains;
					if(to==domj->iproc)
						domi->sendBoundary(domi,domj);
				}
			}
			// Receive boundary
			for(int idom=0;idom<ndomains;idom++)
			{	Domain	*dom=domain_root+idom;
				if(iproc==dom->iproc)
				for(int jdom=1;jdom<ndomains;jdom++)
				{	Domain *domj=domain_root+(idom+jdom)%ndomains;
					if(from==domj->iproc)
						dom->receiveBoundary(dom, domj);
				}
			}
		}
	}
	for(int i=0;i<ndomains;i++)
	{	Domain	*dom=domain_root+i;
		if(iproc==dom->iproc)
			dom->sendOverlap();
	}
	for(int i=0;i<ndomains;i++)
	{	Domain	*dom=domain_root+i;
		if(iproc==dom->iproc)
			dom->receiveOverlap();
	}
#endif

	//Set boundary face-flags

	for (int i=0; i<ndomains; i++)
	{	Domain	*dom=domain_root+i;
		if(iproc==dom->iproc)
			dom->setBoundaryFaceTypes();
	}	

/*****/

}
#endif
void	setup()
{	using namespace Input;
	struct timeb	worldtime;
	ftime(&worldtime);
	//srandom(worldtime.time);//seed random numbers
	srandom(777);//DEBUG
	//set the mask parameter used to retrieve the neighbor cell-face index
	mask[0]=3;
	for (int i=1; i<Nv; i++)
		mask[i]=mask[i-1]<<2;
	//Times and time steps 
	getTime(::configfile);
	if(option.debug|option.verbose)printf("Initializing domains\n");
	ndomains=countObjects((char*)DOMAIN_KEYWORD,::configfile);
	if(ndomains<=0) ERROR1("No keyword found:",DOMAIN_KEYWORD);
	if(option.debug|option.verbose)
	{	printf("Allocating memory for %d domains\n",ndomains);fflush(stdout);
	}
	domain_root=new Domain[ndomains];
	if(option.debug|option.verbose){printf("Setting-up domains\n");fflush(stdout);}

	for (int i=0; i<ndomains; i++)
		domain_root[i].setDomain();
}
void	init()
{
	for (int i=0; i<ndomains; i++)
		domain_root[i].init(domain_root+i);
#ifdef WITH_DOVE
	connectall();
#endif
}
int	run(int niter)
{
	for (int i=0;i<ndomains;i++)
	{	Domain	*dom=domain_root+i;
#ifdef WITH_DOVE
		if(iproc==dom->iproc)
			if (dom->Steps(niter)==0)return 0;
#else
		if (dom->Steps(niter)==0)return 0;
#endif
	}
#ifdef WITH_DOVE
	//Updating variables on domain overlaps
	//This processor
	for(int idom=0;idom<ndomains;idom++)
	{	Domain	*dom=domain_root+idom;
		if(iproc==dom->iproc)
				dom->updateDove();
	}
	//Other processors
	//Send Dove
	for(int i=0;i<ndomains;i++)
	{	Domain	*dom=domain_root+i;
		int	maxvar=dom->getNoVariables();
		if(iproc==dom->iproc)
		for(int ivar=0;ivar<maxvar;ivar++)
			dom->sendDove(ivar);
	}
	//Receive Dove
	for(int i=0;i<ndomains;i++)
	{	Domain	*dom=domain_root+i;
		int	maxvar=dom->getNoVariables();
		if(iproc==dom->iproc)
		for(int ivar=0;ivar<maxvar;ivar++)
			dom->receiveDove(ivar);
	}
#endif
	return 1;//DBUG
}

/********************************
when there is no gui, get variable
 selection from variables.config*/
void	selectVariables()
{
	int	i,ivar,nvar, var,number[3],idom, icom;
	FILE *inp;

//	if (iproc==0)
//	{
		OPENREAD("variables.config",inp);
		SKIPLINE(inp);
//	}	
	for(int j=0; j<ndomains; j++)
	{
	//	if(iproc==0)
	//	{
			fscanf(inp, "%d\t\t%d\t\t%d\n", &number[0],&number[1], &number[2]);
			printf( "variables.config :%d\t\t%d\t\t%d\n",  number[0],number[1], number[2]);
	//	}

	//	COMM_WORLD.Bcast(number, 3, MPI_INT, 0);

		idom = number[0]-1;
		ivar = number[1]-1;
		icom = number[2]-1;
		domain_root[idom].setDisplayVariable(ivar);

		if (ivar<0||domain_root[idom].variable[ivar].dimension==0)
			domain_root[idom].setDisplayVariableComp(0);
		else
			domain_root[idom].setDisplayVariableComp(icom);
	}
	fclose(inp);
}

void	advance1step()
{
	run(1);
	printf("\rTime =%11.4f ",runtime.current);fflush(stdout);

	/* this block is using process 0 to save data of all domains
		commented out since all domains can do saving directly
		without send data back to process 0

	if(iproc==0)
	{
		printf("\rTime =%11.4f ",runtime.current);fflush(stdout);

		double x = fabs(runtime.current-(int)runtime.current);

		if(x<runtime.step)
		{	char outfilename[30];
			for(int i =0; i<ndomains; i++)
			{
				Domain *dom = domain_root + i;
				int ivar = dom->getDisplayVariable();
				sprintf(outfilename, "%s_%s", dom->name, dom->variable[ivar].name);
				printf("filename: %s\n", outfilename);
				dom->saveDispVarDefault(ivar,outfilename);
			}
		}
	}
	*/
	/*replace the above*/
#ifndef WITH_GUI
		//double x = fabs(runtime.current-(int)runtime.current);
		double x = fmod(runtime.current, 1.0);
		if(x<runtime.step)
		{	char outfilename[30];
			for(int i =0; i<ndomains; i++)
			{
				Domain *dom = domain_root + i;
#ifdef WITH_DOVE
				if(dom->iproc==iproc)
#endif
				{
					int ivar = dom->getDisplayVariable();
					sprintf(outfilename, "%s_%s", dom->name, dom->variable[ivar].name);
//					sprintf(outfilename, "%s_%s_%f_", dom->name, dom->variable[ivar].name, runtime.current);
					printf("Saving data in file: %s\n", outfilename);
					dom->saveDispVarDefault(ivar,outfilename);
				}
			}
		}
#endif
}

/*************************************/

int	main(int argc, char *argv[])
{
#ifdef WITH_MPI
	MPI::Init(argc, argv);	
	iproc=COMM_WORLD.Get_rank(),//this processor number
	cout<<"my rank:"<<iproc<<endl;
	nproc=COMM_WORLD.Get_size();//total number of processors 
#else
#ifdef WITH_DOVE
	iproc=0;
	nproc=1;
#endif
#endif

	strcpy(::configfile,CONFIGFILE);
	readcmdline(argc, argv);

	if(option.verbose)printf("%s\n",ABOUT);
	if(option.debug){printf("Debugging output is ON\n");fflush(stdout);}

printf("main:Gui:setup:variable=%d\n",Gui::variable);fflush(stdout);//-
	setup();
printf("main:Gui:init:variable=%d\n",Gui::variable);fflush(stdout);//-
	init();

#ifdef WITH_GUI
#ifdef WITH_DOVE
	if (iproc==0)
#endif
	{
printf("main:Gui:initgui:variable=%d\n",Gui::variable);fflush(stdout);//-
		initgui();
		guirun(argc,argv);
		cleanup();
	}
#ifdef WITH_DOVE
	else
#endif
	{
#ifdef WITH_MPI
		int *vars = new int[ndomains];
		int *coms = new int[ndomains];
		COMM_WORLD.Recv(vars, ndomains, MPI_INT, 0, 100);
		COMM_WORLD.Recv(coms, ndomains, MPI_INT, 0, 101);
		for(int j=0; j<ndomains; j++)
		{
//printf("bcast recvd vars[%d]=%d, coms[%d]=%d", j, vars[j], j, coms[j]);
			domain_root[j].setDisplayVariable(vars[j]);
			domain_root[j].setDisplayVariableComp(coms[j]);
		}
		COMM_WORLD.Barrier();
		while(1) advance1step();
#endif
	}
#else
	selectVariables();
	while(1) advance1step(); //should use a criteria here; not a infinite loop
#endif
	
#ifdef WITH_MPI
	MPI::Finalize();
#endif

	return 0;
}
void	cleanup()
{
//	delete particles;
//	delete field;
}
void	randvec
//Returns a random unit vector
(
	double	*e //Random unit vector
)
{	double	r,r2;
	do
	{	r2=0.0;
		for (int i=0; i<DIM; i++)
		{	double	rnd=2.*RND-1.;
			r2+=rnd*rnd;
			e[i]=rnd;
		}
	}	while (r2>1.0&&r2<SMALL);
	r=sqrt(r2);
	for (int i=0; i<DIM; i++) e[i]/=r;
}

