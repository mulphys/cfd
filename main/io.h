#define	DOCTYPE	"mulphys"
#define	TAB	"   "
#define	CONFIGFILE	"mulphys.xml"
#define	GUI_KEYWORD	"GUI"
#define	DOMAIN_KEYWORD	"domain"
#define	DOMAIN_NAME_TAG	"id"
#define	MODEL_KEYWORD	"model"
#define	MODELS_KEYWORD	"models"
#define	ID_TAG	"id"
#define	MODEL_ID_TAG	"id"
#define	TYPE_TAG	"type"
#define	BOUNDS_TAG	"bounds"
#define	DOMAIN_TYPE_KEYWORD	"type"
#define	DOMAIN_PROC_TAG	"proc"
#define	POINTS_KEYWORD	"points"
#define	POINTS_TAG	"particles"
#define	DOMAIN_OPERATION_KEYWORD	"Do"
#define	TRANSFORM_TAG	"transform"
#define	ROTATE_TAG	"rotate"
#define	ROTATEX_TAG	"rotatex"
#define	ROTATEY_TAG	"rotatey"
#define	ROTATEZ_TAG	"rotatez"
#define	TRANSLATE_TAG	"translate"
#define	TIME_KEYWORD	"time"
#define	TIME_START_KEYWORD	"start"
#define	TIME_END_KEYWORD	"end"
#define	TIME_STEP_KEYWORD	"step"
#define	NUMBER_TAG	"number"
#define	GEOMETRY_TAG	"geometry"
#define	MESH_TAG	"mesh"
#define	GRID_TAG	"grid"
#define	SHAPE_TAG	"shape"
#define	NODES_TAG	"nodes"
#define	CELLS_TAG	"cells"
#define	COORDINATES_TAG	"coordinates"
#define	VARIABLE_KEYWORD	"Var"
#define	VARIABLE_TYPE_KEYWORD	"type"
#define	VARIABLE_RANK_KEYWORD	"rank"
#define	VARIABLE_SPAN_KEYWORD	"span"
#define	VARIABLE_INIT_KEYWORD	"init"
//#define	LOCAL_VARIABLE_KEYWORD	"var"
#define	PROCEDURE_KEYWORD	"Procedure"
#define	PROCEDURE_TYPE_KEYWORD	"type"
#define	PROCEDURE_NUMBER_KEYWORD	"num"
#define	DIMENSIONS_TAG	"dimensions"
#define	GRIDTYPE_KEYWORD	"grid.type"
#define	GRIDSIZE_KEYWORD	"grid.size"
#define	GRIDDIMENSIONS_TAG	"grid.dimensions"
#define	GRIDFILE_KEYWORD	"grid.file"
#define	FILENAME_TAG	"file"
#define	BOUNDARY_CONDITIONS_TAG	"boundary"
#define	TOOL_TAG	"tool"
///#define	MAXNODES_KEYWORD	"maxnodes"
///#define	MAXCELLS_KEYWORD	"maxcells"
#define	FRAME_KEYWORD	"Frame"
#define	CONNECTIVITY_KEYWORD	"connect"
#define	ADD_KEYWORD	"add"
#define	SUB_KEYWORD	"sub"
#define	GET_KEYWORD	"get"
#define	MUL_KEYWORD	"mul"
#define	DIV_KEYWORD	"div"
#define	OUTPUT_TYPE_KEYWORD	"output.type"
#define	OUTPUT_TAG	"output"
#ifndef double
#define	double	double
#endif
#ifndef	MAXLINLEN
#define	MAXLINLEN	510
#endif
/* Error checking malloc.  NOTE THE EXIT CALL ON FAILURE! */
#ifndef ALLOC
#define ALLOC(ptr,type,_N)	\
	if((ptr = (type *)malloc((unsigned)(_N)*sizeof(type))) == NULL) \
	{ fprintf(stderr,"malloc failed on request for %d bytes (at line %d in %s)\n",\
	((int)((_N)*sizeof(type))),__LINE__,__FILE__); exit(1);}
#define ALLOC0(ptr,type,_N)                                            \
	if((ptr = (type *)malloc((unsigned)(_N)*sizeof(type))) == NULL)     \
	{  fprintf                                                          \
		(	stderr,                                                       \
			"malloc failed on request for %d bytes (at line %d in %s)\n", \
			((int)((_N)*sizeof(type))),__LINE__,__FILE__                  \
		);                                                               \
		exit(1);                                                         \
   };                                                                  \
	for (int i=0; i<_N; i++)                                            \
		ptr[i]=(type) 0;

#endif
#define	FLUSH	fflush(stdout)
#define	FILESHORT \
	{fprintf(stderr,"FILE %s IS TOO SHORT\n",inpfilename); exit(1);}
#define	OPENREAD(flname,inp)	\
	if ((inp=fopen(flname,"r"))==NULL) \
	{	fprintf(stderr,"CAN'T OPEN %s\n",flname); exit(1); }
#define	OPENWRITE(flname,out)	\
	if ((out=fopen(flname,"w"))==NULL) \
	{	fprintf(stderr,"CAN'T OPEN %s\n",flname); exit(1); }

/*	while((c=(char)fgetc(fp))!=EOF&&(isspace(c)||strchr(",;:",c)==NULL)); */

#define	GETWORD(buf,fp)	\
	{	char c; int i=0; 	\
	while((c=(char)fgetc(fp))!=EOF&&(isspace(c)||c==':'||c=='='||c=='('||c==')')); \
	do	{	buf[i++]=c;	} \
 	while((c=(char)fgetc(fp))!=EOF&&!isspace(c)&&c!='='&&c!=':'&&c!='('&&c!=')'&&i<MAXLINLEN); \
	buf[i]='\0'; }
#define	GETFWORD(buf,fp)	\
	{	char c; int i=0; 	\
	while((c=(char)fgetc(fp))!=EOF&&(isspace(c)||c==':'||c=='='||c==')')); \
	do	{	buf[i++]=c;	} \
 	while((c=(char)fgetc(fp))!=EOF&&!isspace(c)&&c!=':'&&c!='='&&c!=')'&&i<MAXLINLEN); \
	buf[i]='\0'; }

#define	SKIPWORD(fp)	\
	{	char c; int i=0; 	\
	while((c=(char)fgetc(fp))!=EOF&&isspace(c)); \
	do	{	i++;	} \
 	while((c=(char)fgetc(fp))!=EOF&&!isspace(c)&&i<MAXLINLEN); \
	if ((int)c==EOF) FILESHORT;}
#define	GETINT(i,fp)	{char buf[MAXLINLEN];GETWORD(buf,fp);sscanf(buf,"%d",&i);}
#define	GETREAL(r,fp)	{char buf[MAXLINLEN];GETWORD(buf,fp);sscanf(buf,"%lg",&r);}
#define	SKIPLINE(fp) \
	{	char c; int i=0; 	\
		while((c=(char)fgetc(fp))!=EOF&&c!='\n'&&i<MAXLINLEN); \
	}

/*		if ((int)c==EOF) FILESHORT; */
#define	GETLINE(buf,fp) \
	{	char c; int i=0; 	\
		while((c=(char)fgetc(fp))!=EOF&&isspace(c)); \
		do { buf[i++]=c; }	\
		while((c=(char)fgetc(fp))!=EOF&&c!='\n'&&i<MAXLINLEN); \
		buf[i]='\0'; \
		if (i==MAXLINLEN) \
		{	fprintf (stderr,"\nLINE TOO LONG: %s <-- TRUNCATED\n",buf); \
		}	\
	}
#define PARSEWORD(line,word) \
		{word=line;while (!isspace(*word)&&*word!='['&&*word!='=')word++;*word++='\0';}
#define	SCOPY(a,b)	strncpy(a,b,MAXLINLEN)

#define ERROR(message) \
		{	fprintf(stderr,"\nERROR: %s\n",message);exit(1);}
#define ERROR1(message,line) \
		{	fprintf(stderr,"\nERROR: %s '%s'\n",message,line);exit(1);}
#define BUG(message) \
		{	fprintf(stderr,"\nERROR: %s\nPLEASE, SUBMIT BUG REPORT\n",message);exit(1);}
#define WARNING(message) \
		{	fprintf(stderr,"\nWARNING: %s\n",message);}
#define WARNING1(message,name) \
		{	fprintf(stderr,"\nWARNING: %s\n",message,name);}
#define	CHOPIND(word)	{char *p;for (p=word; *p!='\0'&&*p!='['; p++); *p='\0';}

//#define	INDENT	for (int i=0;i<level;i++)printf(TAB);
#define	INDENT(i)	for (int j=0;j<level+(i);j++)printf(TAB);
enum	OutputTypes
{
	default_output=0,
	user_output,
	tecplot_output,
	ensight_output,
	maxoutypes
};
extern int	getrank(char	*word);
extern void	setrank
(
	char	*word,
	int	*Rank,
	int	**ind
);
extern void	setind
(
	char	*word,
	int	rank,
	int	**ind
);
extern int	debug;///DDD


