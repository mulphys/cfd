#define	ABOUT	"\nProgram HEDRA: Contiuum and Discrete Media Solver (demo-version)\n \
	Author: Andrei SMIRNOV, andrei@smirnov.mae.wvu.edu\n"

#define DEBUG

#define VERSION	1

#define GUI
//#define GRIDGEN
#define	FEM

#ifdef PARTICLE_INTERACTION
#define MOMENTUM_EXCHANGE
/*#define COALESCENCE*/
#endif

#define	TETRA

#ifndef	MAXLINLEN
#define	MAXLINLEN	510
#endif
/*
#ifndef	INT
#define	INT	int
#endif
#ifndef	REAL
#define	REAL	double
#endif
*/
#ifndef	DIM
#define	DIM	3
#endif
#define	SMALL	10e-17
#define	LARGE	10e+30

#define	RND	random()/RAND_MAX

#ifndef ALLOC
#define ALLOC(ptr,type,_N)  \
	if((ptr = (type *)malloc((unsigned)(_N)*sizeof(type))) == NULL) {  \
	fprintf(stderr,"malloc failed on request for %d bytes (at line %d in %s)\n", \
	((int)((_N)*sizeof(type))),__LINE__,__FILE__); \
	exit(1);  \
	}

#define ALLOC0(ptr,type,_N)                                            \
	if((ptr = (type *)malloc((unsigned)(_N)*sizeof(type))) == NULL)     \
	{  fprintf                                                          \
		(	stderr,                                                       \
			"malloc failed on request for %d bytes (at line %d in %s)\n", \
			((int)((_N)*sizeof(type))),__LINE__,__FILE__                  \
		);                                                               \
		exit(1);                                                         \
   };                                                                  \
	{int i; for (i=0; i<_N; i++)                                        \
		ptr[i]=(type) 0; }

#endif

struct	Option
{
	unsigned int	verbose	:1;
	unsigned int	debug	:1;
};
struct RunTime
{	double	start,end,prev,current,step,step0;
};
extern int	ndomains;
#ifdef WITH_DOVE
extern int	iproc;//number of process of this domain
#endif
extern class Domain	*domain;
extern struct RunTime	runtime;
extern char	programname[];
extern char	configfile[];
extern Option option;
extern void	animate();
//?extern int	countObjects
//?	(	char	*keyword,
//?		char	*filename
//?	);
void	display();
void	reshape(int w, int h);
void	cleanup();
void	init();
void	randvec(double *e);
int	run(int niter);
//?void	getTime(char *filename);
