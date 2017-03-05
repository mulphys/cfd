
/*
 *		    Copyright (c) 1994, 1995, 1996
 *	    Computational Engineering International, Inc.
 *			 All Rights Reserved.
 */

#ifndef double
#define double	double
#endif

/* Output types */
#define ASCII	1
#define BINARY	2	/* default */
extern int EnSightFileType;

/* EnSight element types */
#define BAR2	0
#define BAR3	1
#define TRIA3	2
#define TRIA6	3
#define QUAD4	4
#define QUAD8	5
#define TETRA4	6
#define TETRA10	7
#define HEXA8	8
#define HEXA20	9
#define PENTA6	10
#define PYRMD5	11
#define PYRMD13	12

/* Node/element ID numbering schemes */
#define ID_UNKNOWN	0
#define ID_OFF		1
#define ID_GIVEN	2
#define ID_ASSIGN	3
#define ID_IGNORE	4

/* Error checking malloc.  NOTE THE EXIT CALL ON FAILURE! */
#define allocate(ptr,type,_N)  \
    if((ptr = (type *)malloc((unsigned)(_N)*sizeof(type))) == NULL) {  \
	fprintf(stderr,"malloc failed on request for %d bytes (at line %d in %s)\n",((int)((_N)*sizeof(type))),__LINE__,__FILE__);  \
        exit(1);  \
    }

typedef struct {
    char *name;
    int numnodes;
} ElemType;
extern ElemType elemTypes[];

/* Coodinates */
typedef struct {
    int num;			/* number of nodes */
    int *id;			/* node IDs (null if ID_OFF or ID_ASSIGN) */
    double *coords;		/* coordinates */
} Coords;

/* ElementSet: a set of elements of the same type */
typedef struct {
    int type;			/* element type */
    int numelems;		/* number */
    int *id;			/* element IDs (null if ID_OFF or ID_ASSIGN) */
    int *nindx;			/* array of node indices */
} ElementSet;

/* Part structure */
typedef struct {
    char *descrip;		/* part description */
    int numelemsets;		/* number of element sets present */
    ElementSet **elemSets;	/* element sets */
} Part;

/* Geometry structure: header, coords, parts */
typedef struct {
    char *descrip1,*descrip2;
    int nodeID;			/* node ID numbering scheme */
    int elemID;			/* element ID numbering scheme */
    Coords *coords;		/* coordinates */
    int numparts;		/* number of parts */
    Part **parts;		/* array of part pointers */
} Geometry;

/* Particle geometry structure: header, coords */
typedef struct {
    char *descrip1,*descrip2;
    Coords *coords;		/* coordinates (includes IDs) */
} ParticleGeometry;

typedef struct {
    char *descrip;
    char *varname;
    int num;
    double *scl;
} Scalar;

typedef struct {
    char *descrip;
    char *varname;
    int num;
    double *vec;
} Vector;

typedef struct {
    int nscl,nvec;		/* number of scalars and vectors */
    int geochange;		/* geometry changing flag */
    int nsteps;			/* number of time steps */
    double *time;		/* actual time values */
    int start,inc;		/* start and increment for file name numbers */
    char *geofname;		/* geometry (if changing) wildcard name */
    char **sclfname;		/* scalar file name(s) */
    char **sclvname;		/* corresponding scalar variable name */
    char **vecfname;		/* vector file name(s) */
    char **vecvname;		/* corresponding vector variable name */
} Result;


/* Callable routines in in.c */

Geometry *ReadGeoHeader(char *);
Geometry *ReadGeometry(char *);
ParticleGeometry *ReadParticleGeometry(char *);
Coords *ReadCoords(int, FILE *, int);
Part **ReadParts(FILE *, int *, int);


/* Callable routines in out.c */

void SetFileType(int);

int WriteGeometry(char *, Geometry *);
int WriteGeoHeader(char *, char *, char *, int, int);
void WriteGeoCoords(int, int, int *, double *);
void WriteGeoPart(char *);
void WriteGeoElem(int, int, int *, int *);

int WriteParticleGeometry(char *, ParticleGeometry *);
int WriteParticleGeoHeader(char *, char *);


/* Callable routines in results.c */

Scalar *ReadScalar(char *, int);
Vector *ReadVector(char *, int);

int WriteScalar(char *, Scalar *);
int WriteRawScalar(char *, char *, char *, int, double *);

int WriteVector(char *, Vector *);
int WriteRawVector(char *, char *, char *, int, double *);

Result *ReadResults(char *);
int WriteResults(char *, Result *);
int dumpResults(char *);

Scalar *newScalar(char *, char *, int);
void freeScalar(Scalar *);

Vector *newVector(char *, char *, int);
void freeVector(Vector *);


/* Useful stuff */
#ifndef True
#define True	1
#endif /* True */
#ifndef False
#define False	0
#endif /* False */


/* end eio.h */
