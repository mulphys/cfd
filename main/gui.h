/*
 * OpenGL(TM) is a trademark of Silicon Graphics, Inc.
 */

//#include <stdlib.h>
//#include <stdio.h>
//#include <ctype.h>
//#include <math.h>
//#include <string.h>
//#include <GL/glut.h>
//#include "vecalg.h"
//#include "io.h"
//#include "main.h"
//#include "geom.h"
//#include "var.h"
//#include "domain.h"
////#include "gui.h"
//#include "output.h"

#define NEEDLE

#define	WINDOW_SIZE	800
#define	MAX_WINDOW_WIDTH	800
#define	MAX_WINDOW_HEIGHT	600

#define	LEFT_BUTTON	1
#define	MIDDLE_BUTTON	2
#define	RIGHT_BUTTON	3

namespace	Gui
{
	enum	Color
	{
		red=0,
		green,
		blue,
		skyblue,
		brown,
		magenta,
		color6,
		color7,
		color8,
		color9,
		color10,
		color11,
		color12,
		color13,
		color14,
		color15,
		maxcolor
	};
	enum	PointParameters
	{
		variable=maxcolor, 
		vmin,vmax,size,
		sizevar,massvar,
		maxpointprm
	};
	enum	LineParameters
	{
		thickness=maxcolor, 
		length,
		maxlineprm
	};
	enum	ToolParameters
	{
		maxtoolprm=maxcolor
	};
	enum	AxesParameters
	{
		axeswidth=0,
		xaxislength,
		yaxislength,
		zaxislength,
		arrowheight,
		arrowidth,
		maxaxesprm
	};
	enum	SurfDispType
	{
		gridlines=0,
		solidsurface,
		maxsurfdisptypes
	};
	extern char	configfile[MAXLINLEN];
	extern int
		idomain;//current domain number
		//nvar,//maximum number of variables
	extern double
		xmin[DIM],xmax[DIM],
		lx,ly,lz,lmin,lmax;
	extern float
		step,
		vecval[maxlineprm],
		axes[maxaxesprm],
		rgbcolor[maxcolor][DIM];

	struct	ColorScale
	{
		int	relative	:1;
	};
	struct	Grid
	{
		float
			node[maxpointprm],
			line[maxlineprm];
	};
	struct	Particle
	{
		int
			m,     /* maximum number of particles           */
			n;     /* number of particles                   */
			//nv, iv;/* number of variables, current variable */
		float
			disp[maxpointprm];
		double
			*X;    /* coordinates */
//			*V,    /* vectors     */
//			*S;    /* scalars     */
		struct Point	*first,*last,*origin;/* pointers to the first and last particles */
	};

//	public:
	enum Movement
	{	stay=0,
		rotate,
		moveuv,
		movew
	};
	enum ShowPars
	{	showRun=0,
		showAxes,
		showTool,
		showParticles,
		showSpheres,
		showBonds,
		showGrid,
		showNodes,
		showVariables,
		showBoundaryVertexes,
		showBoundaryVectors,
		showToolVertexes,
		showBoundaryFaceCenters,
		showBoundaryFaces,
		showBoundaryGrid,
		showToolGrid,
		showFrame,
		showCellCenters,
		showFaceCenters,
		showIsoSurfaces,
		dumpWindow,
		maxshowpars
	};
	extern int
		showpar[maxshowpars],
		finished,
		animation;
	extern float
		xo[DIM], //origin of coordinate system
		dx,dy,dz,
		lastx, lasty,	lastz;
	struct	WindowGeom
	{
		int
			width,height;
	};

	void	inigui();
	void	Exit();
	void	Quit();
	void	finish();
	void	readconf();
	void	helpDisplay();
	int	readparam
	(
		char	*s,
		char	*param[],
		int	maxparam,
		char	*filename,
		float	val[]
	);
	void	initdisp();
	void	myinit(int argc, char *argv[]);
	void	InitMaterials(void);
	void	display();
	void	displayAxes();
	void	displayDomain(int idomain);
	void	writeData();
	void	helpCommand();
	void	printGridXLimits();
	void	printGridVecLimits();
	void	getXLimits(double	*xmin, double *xmax);
	void	getScaling(double	&lmin, double &lmax);
	void	printPartVelLimits();
	void	printDomains();
	void	selectCurrentDomain();
	void	printParticleVariables();
	void	showGridElements
	(	
		int ne, 
		double	*X,//coordinates
		float color[]
	);
	void	showVectorComponent
	(
		int	ivar,
		int	icomp,
		double	*X,//coordinates
		double	vmn,
		double	vmx
	);
	void	showVector
	(
		int	ivar,
		double	*X//coordinates
	);
	void	commandMode();
	void	printVariables
	(
		int	n
//		struct Variables	*V
	);

extern struct Particle	P;

void	initgui();
void	finish();
void	helpDisplay();
int	readparam
(
	char	*s,
	char	*param[],
	int	maxparam,
	char	*filename,
	float	val[]
);
void	readconf();
void	initdisp();
void	InitMaterials(void);
/*  Initialize material property and depth buffer.
 */
void	myinit(int argc, char *argv[]);
void	showGridElements
(	int	ne,//number of elements to show
	double	*X,//coordinates
	float	color[]
);
void	showVector
(
	int	ivar,
	double	*X//coordinates
);
void	showVectorComponent
(
	int	ivar,
	int	icomp,
	double	*X,//coordinates
	double	vmn,
	double	vmx
);
void	displayAxes();
void	displayDomain(int idomain);
/* ARGSUSED1 */

void	getScaling
(
	double	&l0,
	double	&l1
);
void	getXLimits
(	double	*x0,
	double	*x1
);
#ifndef GRIDGEN
void	printPartVelLimits();
#endif
void	printVariables
(
	int	n
//	struct Variables	*V
);
void	Exit();
void	Quit();
void	printDomains();
void	selectCurrentDomain();
void	selectVariable();
void	printParticleVariables();
void	setBackgroundRun();
void	setForegroundRun();
void	commandMode();

/* Main */

/* Notes:

Alternative grid color:

quadric.c:
   glDisable(GL_LIGHTING);
   glColor3f(0.0, 1.0, 1.0);

Isosurface is based on Examples
isosurf.c
smooth.c

aapoly.c
tri.c
wave.c

 */

/******** GLOBAL GUI-FUNCTIONS **********/

void	display();
void	animate();
void	writeData();
void	printGridXLimits();
void	printGridVecLimits();
void	helpCommand();
void	reshape(int w, int h);
void	mouse(int button, int state, int x, int y);
void	motion(int x, int y);
/* ARGSUSED3 */
void	keyboard(unsigned char key, int x, int y);
void	menu(int value);
void	guirun(int argc, char* argv[]);
void	drawSegment
(
	float	*x,
	float	*y 
);
#ifdef NEEDLE
struct Needle
{
	float	x[DIM],y[DIM];
};
#endif
void	displaymessage(char *);
} //END NAMESPACE
