#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "io.h"
#include "main.h"
#include "geom.h"
#include "var.h"
//#include "particles.h"
#include "tool.h"
#include "domain.h"

void	Variable::setinit(int dim)
{
	dimension=dim;
	ALLOC(initval,double,dim);
	ALLOC(min,double,dim);
	ALLOC(max,double,dim);
	for (int i=0; i<dim; i++)
		initval[i]=min[i]=max[i]=0.0;
}
void	Variable::delinit()
{
	dimension=0;
	free(initval);
	free(min);
	free(max);
}
void	Variable::init(int size)
{
	if (dimension<=0||size<=0)
		ERROR("init: Rank or size of a variable is undefined");
	if (option.verbose|option.debug)
	{	INDENT(1);
		printf
		(	"Allocating %d*%d*%d=%d bytes\n",
			size,dimension,sizeof(double),size*dimension*sizeof(double)
		);
	}
	ALLOC(val,double,size*dimension);
	reset();
}
void	Variable::reset()
{
	if (dimension<=0||size<=0)
		ERROR("reset: Rank or size of a variable is undefined");
	for (int i=0; i<size; i++)
	for (int j=0; j<dimension; j++)
		val[i*dimension+j]=initval[j];
}
//	void	Variable::setvar
//	(	char *word, 
//		int *Level,
//		int	parent_procedure_type,
//		int	*Ne,//number of grid-span elements of different types
//		FILE *inp,
//		char	*inpfilename
//	)
//	{	//This function is invoked for grid-variables
//		int	level0;
//		char	*p;
//		level=*Level,level0=level;
//		if (parent_procedure_type==-1) ERROR1("Parent type undefined for",word);
//		//span=grid;
//		type=-1;
//		rank=nrank=-1;
//		size=-1;
//		initval=NULL;
//		//GETWORD(word,inp);
//		if ((p=strchr(word,'['))==NULL) rank=0;
//		else *p='\0';
//		if (rank==-1)
//		{	rank=1;
//			for (p++; *p!='\0'; p++)
//				if (*p==',') rank++;
//		}
//		SCOPY(name,word);
//		if(option.verbose|option.debug)
//		{	INDENT(0);
//			printf("Variable='%s'\n",name);fflush(stdout);
//			INDENT(1);
//			printf("rank=%d\n",rank);fflush(stdout);
//		}
//		GETWORD(word,inp);
//		while (*word=='!'&&!feof(inp)) {SKIPLINE(inp);GETWORD(word,inp);}
//		if (*word!='{') ERROR("Symbol '{' missing after "VARIABLE_KEYWORD);
//		level++;
//		//Scan variable description
//		while (!feof(inp)&&level>level0)
//		{	GETWORD(word,inp);
//			if (*word=='!') {SKIPLINE(inp);continue;}
//			if (*word=='{') {level++;continue;}
//			if (*word=='}') {level--;continue;}
//			if (level<=level0)break;
//			if (strcmp(word,VARIABLE_TYPE_KEYWORD)==0)
//			{	if(option.verbose|option.debug)
//				{	INDENT(0);
//					printf("%s=",VARIABLE_TYPE_KEYWORD);fflush(stdout);
//				}
//				GETWORD(word,inp);
//				if (isdigit(*word))
//					sscanf(word,"%d",&type);
//				else
//				{	for (int itype=0; itype<DIM+1; itype++)
//						if (strcmp(word,elementype[itype])==0) type=itype;
//				}
//				if(option.verbose|option.debug)
//					printf("%s\n",elementype[type]);fflush(stdout);
//				continue;
//			}
//			//initialize
//			if (strcmp(word,"input")==0)
//			{	GETWORD(word,inp);
//				if (word[0]=='"'&&word[1]=='"')
//					inpfile=NULL;
//				else
//				{	char	*p=word+1,*s=strchr(p,'"');
//					if (s==NULL) ERROR("Malformed input file name");
//					*s='\0';
//					inpfile=strdup(p);
//					if(option.verbose|option.debug)
//					{	INDENT(0);
//						printf("inpfile='%s'\n",inpfile);
//					}
//				}
//				continue;
//			}
//			if (strcmp(word,"output")==0)
//			{	GETWORD(word,inp);
//				if (word[0]=='"'&&word[1]=='"')
//					outfile=NULL;
//				else
//				{	char	*p=word+1,*s=strchr(p,'"');
//					if (s==NULL) ERROR("Malformed output file name");
//					*s='\0';
//					outfile=strdup(p);
//					if(option.verbose|option.debug)
//					{	INDENT(0);
//						printf("outfile='%s'\n",outfile);
//					}
//				}
//				continue;
//			}
//			if (strcmp(word,VARIABLE_INIT_KEYWORD)==0)
//			{	if(option.verbose|option.debug)
//				{	INDENT(0);
//					printf("initialized to ");fflush(stdout);
//				}
//				if (type<0||rank<0)
//				{	fprintf
//					(	stderr,
//						"Variable %s and %s should be defined before\n",
//						VARIABLE_TYPE_KEYWORD,
//						VARIABLE_RANK_KEYWORD,
//						VARIABLE_INIT_KEYWORD
//					);exit(1);
//				}
//				nrank=(int)pow(DIM,rank);
//				ALLOC(initval,double,nrank);
//				for (int i=0; i<nrank; i++)
//				{	GETWORD(word,inp);
//					if (!isdigit(*word))
//					{	fprintf(stderr,"\nNumber of initializatin values less than %d\n",nrank);
//						exit(1);
//					}
//					sscanf(word,"%lf",initval+i);
//					if(option.verbose|option.debug)printf(" %lg",initval[i]);fflush(stdout);
//				}
//				if(option.verbose|option.debug)printf("\n");
//				continue;
//			}
//		}
//		if (type==-1) type=parent_procedure_type;
//		size=Ne[type];
//		if (size<=0)
//			ERROR1("Can't determine the size for variable",name);
//		if (nrank<0) setinit(rank);
//		init(size*nrank);
//		*Level=level;
//	}
//	void	Variable::setvar
//	(	char *word, 
//		int *Level,
//		int	parent_procedure_type,
//		FILE *inp,
//		char	*inpfilename
//	)
//	{	int	level0;
//		char	*p;
//		level=*Level,level0=level;
//		if (parent_procedure_type==-1) ERROR1("Parent type undefined for",word);
//		type=-1;
//		rank=nrank=-1;
//		//span=-1;
//		size=-1;
//		//GETWORD(word,inp);
//		if ((p=strchr(word,'['))==NULL) rank=0;
//		else *p='\0';
//		if (rank==-1)
//		{	rank=1;
//			for (p++; *p!='\0'; p++)
//				if (*p==',') rank++;
//		}
//		SCOPY(name,word);
//		if(option.verbose|option.debug)
//		{	INDENT(0);
//			printf("Variable='%s'\n",name);fflush(stdout);
//			INDENT(1);
//			printf("rank=%d\n",rank);fflush(stdout);
//		}
//		GETWORD(word,inp);
//		while (*word=='!'&&!feof(inp)) {SKIPLINE(inp);GETWORD(word,inp);}
//		if (*word!='{') ERROR("Symbol '{' missing after "VARIABLE_KEYWORD);
//		level++;
//		//Scan variable description
//		while (!feof(inp)&&level>level0)
//		{	GETWORD(word,inp);
//			if (*word=='!') {SKIPLINE(inp);continue;}
//			if (*word=='{') {level++;continue;}
//			if (*word=='}') {level--;continue;}
//			if (level<=level0)break;
//			if (strcmp(word,VARIABLE_TYPE_KEYWORD)==0)
//			{	if(option.verbose|option.debug)
//				{	INDENT(0);
//					printf("%s=",VARIABLE_TYPE_KEYWORD);fflush(stdout);
//				}
//				GETWORD(word,inp);
//				if (isdigit(*word))
//					sscanf(word,"%d",&type);
//				else
//				{	for (int itype=0; itype<DIM+1; itype++)
//						if (strcmp(word,elementype[itype])==0) type=itype;
//				}
//				if(option.verbose|option.debug)
//					printf("%s\n",elementype[type]);fflush(stdout);
//				continue;
//			}
//			//determine variable scope
//			if (strcmp(word,VARIABLE_SPAN_KEYWORD)==0)
//			{	if(option.verbose|option.debug)
//				{	INDENT(0);
//					printf("%s=",VARIABLE_SPAN_KEYWORD);fflush(stdout);
//				}
//				GETWORD(word,inp);
//				if (isdigit(*word))
//					sscanf(word,"%d",&span);
//				else
//				{	for (int ispan=0; ispan<maxspantypes; ispan++)
//						if (strcmp(word,spantype[ispan])==0) span=ispan;
//				}
//				if(option.verbose|option.debug)printf("'%s'\n",spantype[span]);fflush(stdout);
//				continue;
//			}
//		}
//		if (type==-1) type=parent_procedure_type;
//		if (span==-1)
//		{	switch (level)
//			{	case 1:
//					ERROR1("Can't define a local variable on level 1:",name);
//					//span=grid;
//					break;
//				case 2:
//					if (parent_procedure_type==type) span=0;
//					else span=1;
//				break;
//			}
//		}
//		size=getvarsize(type,parent_procedure_type,span);
//		if (size<=0)
//			ERROR1("Can't determine the size for variable",name);
//		if (nrank<0) setinit(rank);
//		init(size*nrank);
//		*Level=level;
//	}
//	void	LocalVariable::setref
//	(	char *word, 
//		int	*Level,
//		int	parent_procedure_type,
//		int	nvar,
//		Variable	*variable,
//		FILE *inp,
//		char	*inpfilename
//	)
//	{	int	iref,level=*Level,level0=level;
//		char	*p;
//		Variable	*var;//Referenced variable to be located
//		if (parent_procedure_type==-1) ERROR1("Parent type undefined for",word);
//		type=-1;
//		rank=nrank=-1;
//		span=-1;
//		size=-1;
//		initval=NULL;
//		if ((p=strchr(word,'['))==NULL) rank=0;
//		else *p='\0';
//		if (rank==-1)
//		{	rank=1;
//			for (p++; *p!='\0'; p++)
//				if (*p==',') rank++;
//		}
//		//Locate the name in the list of parent variables
//		for (iref=0; iref<nvar; iref++)
//		if (strcmp(word,variable[iref].name)==0)
//		{	if (rank!=variable[iref].rank)
//				ERROR1("Declared ranks do not match for reference and reffered varaible",word);
//			var=variable+iref;//ifoundvar=1;
//			break;
//		}
//		if (iref==nvar)
//			ERROR1("No match found for reference:",word);
//		SCOPY(name,word);
//		type=var->type;
//		rank=var->rank;
//		initval=NULL;
//		if(option.verbose|option.debug)
//		{	INDENT(0);printf("Variable='%s'\n",name);fflush(stdout);
//			INDENT(1);printf("%s=%d\n",VARIABLE_TYPE_KEYWORD,type);fflush(stdout);
//			INDENT(1);printf("rank=%d\n",rank);fflush(stdout);
//		}
//		//Skip to the body-part (after '{')
//		GETWORD(word,inp);
//		while (*word=='!'&&!feof(inp)) {SKIPLINE(inp);GETWORD(word,inp);}
//		if (*word!='{') ERROR("Symbol '{' missing after "VARIABLE_KEYWORD);
//		level++;
//		//Scan variable description
//		while (!feof(inp)&&level>level0)
//		{	GETWORD(word,inp);
//			if (*word=='!') {SKIPLINE(inp);continue;}
//			if (*word=='{') {level++;continue;}
//			if (*word=='}') {level--;continue;}
//			if (level<=level0)break;
//			if (strcmp(word,VARIABLE_TYPE_KEYWORD)==0)
//				ERROR1("Can't define %s for referece variable",name);
//			//initialize
//			if (strcmp(word,"input")==0)ERROR1("Can't define input for referece variable",name);
//			if (strcmp(word,"output")==0)ERROR1("Can't define output for referece variable",name);
//			if (strcmp(word,VARIABLE_INIT_KEYWORD)==0)
//				ERROR1("Can't initialize referece variable",name);
//			if (strcmp(word,VARIABLE_SPAN_KEYWORD)==0)
//			{	if(option.verbose|option.debug)
//				{	INDENT(0);
//					printf("%s=",VARIABLE_SPAN_KEYWORD);fflush(stdout);
//				}
//				GETWORD(word,inp);
//				if (isdigit(*word))
//					sscanf(word,"%d",&span);
//				else
//				{	for (int ispan=0; ispan<maxspantypes; ispan++)
//						if (strcmp(word,spantype[ispan])==0) span=ispan;
//				}
//				if(option.verbose|option.debug)printf("'%s'\n",spantype[span]);fflush(stdout);
//				if (span==grid)ERROR("Span for the referece variables can not be 'grid'");
//				continue;
//			}
//		}
//		if (span==-1)
//		{	span=thiselement;
//			if(option.verbose|option.debug)
//				printf("Span not defined - assumed 'element' for variable",name);
//		}
//		switch(span)
//		{
//			case thiselement: size=1; break;
//			case neighbors: size=ConnectivitySize[parent_procedure_type][type];break;
//			default:
//				fprintf
//				(	stderr,
//					"Wrong span %d for a reference variable %s\n",
//					span,name
//				);exit(1);
//		}
//		if (size<=0)
//			ERROR1("Can't determine the size for variable",name);
//		ref=var->val;//points to the referenced variable data
//		ALLOC(val,double*,size);
//		//val[0:size][0:nrank]: accesses the values of the variable
//		*Level=level;
//	}
//	void	Procedure::setzero()
//	{
//		*name='\0';
//		nfunc=0;
//		nvar=0;
//		niter=1;
//	}
//	void	Function::setind(Variable *var[3])
//	{
//		int	i,ni,*ii[4],irank[4],
//			fclass;
//		ni=narg+1;
//		if (ni>4) 
//		{	fprintf
//			(	stderr,
//				"Number of arguments %d is too large for function %s\n",
//				ni,name
//			);exit(1);
//		}
//		//collect both function and argument indexes together
//		for (i=0; i<ni-1; i++)
//		{	ii[i]=vind[i];
//			irank[i]=var[i]->rank;
//		}
//		ii[ni-1]=ind;//function indexes
//		irank[ni-1]=rank;
//		//check on repeated indexes in the first argument
//		for (i=0; i<irank[0]; i++)
//		for (int j=0; j<i; j++)
//		if (ii[0][i]==ii[0][j])
//		{	fprintf
//			(
//				stderr,
//				"Dupplicate indexes are not allowed for the target argument %s of function %s\n",
//				var[0]->name,
//				name
//			);exit(1);
//		}
//		if 
//		(	  strcmp(name,ADD_KEYWORD)==0
//			||strcmp(name,SUB_KEYWORD)==0
//		)
//		{	//Additive class
//			if (rank!=0)
//			{	fprintf
//				(	stderr,
//					"Function %s may only have rank 0: specified rank is %d\n",
//					name,rank
//				);exit(1);
//			}
//			for (int i=1; i<narg-1; i++)
//			if (irank[i]!=irank[0])
//			{	fprintf
//				(	stderr,
//					"All arguments of function %s should have the same rank and type.\n\tVariable %s, rank=%d, type=%s\n\tVariable %s, rank=%d, type=%s\n",
//					name,
//					var[0]->name,irank[0],elementype[var[0]->type],
//					var[i]->name,irank[i],elementype[var[i]->type]
//				);exit(1);
//			}
//			for (int i=1; i<ni-1; i++)
//			for (int j=0; j<irank[i]; j++)
//			if (ii[0][j]!=ii[i][j])
//			{	fprintf(stderr,"ERROR in function %s: variables have different indexes:\n",name);
//				fprintf(stderr,"\t%s: ind=",var[0]->name);
//				for (int k=0; k<irank[0]; k++) printf(" %c",ii[0][k]);printf("\n");
//				fprintf(stderr,"\t%s: ind=",var[i]->name);
//				for (int k=0; k<irank[i]; k++) printf(" %c",ii[i][k]);printf("\n");
//				exit(1);
//			}
//			return;
//		}
//		else
//		if 
//		(	  strcmp(name,GET_KEYWORD)!=0
//			&&strcmp(name,MUL_KEYWORD)!=0
//			&&strcmp(name,DIV_KEYWORD)!=0
//		)
//			ERROR1("Can't recognize function",name);
//		//Multiplicative class
//		//check dummy indexes for add,sub functions
//		//replace ascii indexes: i,j,k,... with numbers: 0,1,2,...
//		//first (target) argument
//		for (i=0; i<irank[0]; i++)
//		{	int	m=ii[0][i],mfound=0;
//			for (int j=1; j<ni; j++)
//			for (int k=0; k<irank[j]; k++)
//			if (m==ii[j][k]) 
//			{	ii[0][i]=ii[j][k]=i;mfound++;
//			}
//			if (mfound==0)
//			{	fprintf
//				(	stderr,
//					"Free index '%c' for argument '%s' not matched in function %s\n",
//					m,var[0]->name,name
//				);exit(1);
//			}
//			if (mfound>1)
//			{	fprintf
//				(	stderr,
//					"Free index '%c' for argument '%s' has more than one matches in function %s\n",
//					m,var[0]->name,name
//				);exit(1);
//			}
//		}
//		for (int j=1; j<ni; j++)
//		for (int k=0; k<irank[j]; k++)
//		{	int	mfound=0,
//			m=ii[j][k];
//			for (int p=1; p<j; p++)
//			for (int q=0; q<irank[p]; q++)
//			if (m==ii[p][q])
//			{	ii[p][q]=i; mfound++; }
//			if (mfound>1)
//			{	fprintf
//				(	stderr,
//					"Dummy index '%c' for argument '%s' has more than one matches in function %s\n",
//					m,var[j]->name,name
//				);exit(1);
//			}
//			if (mfound==0)
//			{	//check if the index is present in the target arg
//				int	n;
//				for (int n=0; n<irank[0]; n++)
//				if (m==ii[0][n]) break;
//				if (m==irank[0]) //the index is not there
//				{	fprintf
//					(	stderr,
//						"Dummy index '%c' for argument '%s' has no matches in function %s\n",
//						m,var[j]->name,name
//					);exit(1);
//				}
//			}
//			if (mfound==1)ii[j][k]=i++;
//		}
//	printf("Function %s: rank=%d, ind=",name,rank);
//	for (int i=0; i<rank; i++)printf(" %d",ind[i]);printf("\n");
//	printf("\tTarget %s: rank=%d, ind=",var[0]->name,var[0]->rank);
//	for (int i=0; i<var[0]->rank; i++)printf(" %d",vind[0][i]);printf("\n");
//	for (int j=1; j<narg; j++)
//	{
//	printf("\tSource %s: rank=%d, ind=",var[j]->name,var[j]->rank);
//	for (int i=0; i<var[j]->rank; i++)printf(" %d",vind[j][i]);printf("\n");
//	}
//	}
//	void	Function::setfun(void)
//	{
//		func=NULL;
//	//	v0=&var[0]->val;
//	//	v1=&var[1]->val;
//	//	if (narg==3)v2=&var[2]->val;
//		if (strcmp(name,ADD_KEYWORD)==0)
//		{	if (var[0]->rank==0)
//			{	if (var[0]->span==grid)
//					func=&add_scalar_to_grid;
//				else
//					func=&add_scalar_to_local;
//			}
//			else
//			{	p0=var[0]->nrank;
//				if (var[0]->span==grid)
//					func=&add_to_grid;
//				else
//					func=&add_to_local;
//			}
//		}
//		if (strcmp(name,SUB_KEYWORD)==0)
//		{	if (var[0]->rank==0)
//			{	if (var[0]->span==grid)
//					func=&sub_scalar_to_grid;
//				else
//					func=&sub_scalar_to_local;
//			}
//			else
//			{	p0=var[0]->nrank;
//				if (var[0]->span==grid)
//					func=&sub_to_grid;
//				else
//					func=&sub_to_local;
//			}
//		}
//		if (func==NULL)
//			ERROR1("Function not implemented:",name);
//	}
//	void	add_scalar
//	(	//Add variable var to this variable
//		int	i, int j, int k,//not used
//		double	*a,
//		double	*b,
//		double	*c
//	)
//	{
//		*a=*b+*c;
//	}
//	void	add
//	(	//Add variable var to this variable
//		int	size, int dummy1, int dummy2,
//		double	*a,
//		double	*b,
//		double	*c
//	)
//	{	for (int i=0; i<size; i++)
//			a[i]=b[i]+c[i];	
//	}
//	void	sub
//	(	//Subtract
//		int	size, int dummy1, int dummy2,
//		double	*a,
//		double	*b,
//		double	*c
//	)
//	{
//		for (int i=0; i<size; i++)
//			a[i]=b[i]-c[i];	
//	}
//	void	sub_scalar
//	(	//Add variable var to this variable
//		int	i, int j, int k,//not used
//		double	*a,
//		double	*b,
//		double	*c
//	)
//	{
//		*a=*b-*c;
//	}
//	void	get
//	(
//		int	type,
//		int	rank,
//		double	*val[3]
//	)
//	{
//	
//	}
void	Variable::computeLimits()
{	
	for (int i=0; i<dimension; i++)
	{	min[i]=LARGE;
		max[i]=-LARGE;
	}
	for (int i=0; i<size; i++)
	{	double	*v=val+dimension*i;
		for (int j=0; j<dimension; j++)
		{	if (min[j]>v[j])min[j]=v[j];
			if (max[j]<v[j])max[j]=v[j];
		}
	}
}
void	Variable::getLimits(int icomp, double &vmin, double &vmax)
{
	if (icomp<0||icomp>=dimension) 
	{	fprintf
		(	stderr,
			"getLimits: Variable %s has %d component, but %d requested\n",
			name,dimension,icomp+1
		);exit(1);
	}
	vmin=min[icomp];
	vmax=max[icomp];
}
