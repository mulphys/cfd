/*
 * Graphic User Interface for MESS
 * Author: Andrei Smirnov http://smirnov.mae.wvu.edu
 */

#ifdef WITH_MPI
#include "mpi.h"
using namespace MPI;
#endif

#include<fstream>
#include<iostream>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
//#include <GL/glut.h>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <GL/glx.h>
#include <GL/glu.h>
#include <X11/keysym.h>
#include "vecalg.h"
#include "io.h"
#include "main.h"
#include "geom.h"
#include "var.h"
#include "tool.h"
#include "domain.h"
#include "gui.h"
#include "output.h"

namespace	Gui
{
	char	configfile[MAXLINLEN],
		windowname[MAXLINLEN];	
	int
		idomain,//current domain number
		//nvar,//maximum number of variables
		showpar[maxshowpars],
		finished,
		animation,
		nwdump,iwdump,
		attributeList[] = 
		{	GLX_RGBA, 
			GLX_RED_SIZE, 1, 
			GLX_GREEN_SIZE, 1,
			GLX_BLUE_SIZE, 1, 
			GLX_DOUBLEBUFFER, 
			GLX_DEPTH_SIZE, 1, 
			None 
		};
	float
		step,
		vecval[maxlineprm],
		axes[maxaxesprm],
		toolval[maxtoolprm],
		rgbcolor[maxcolor][DIM],
		wdtime,
		xo[DIM], //origin of coordinate system
		dx,dy,dz,
		lastx, lasty,	lastz;
	double
		xmin[DIM],xmax[DIM],
		lx,ly,lz,lmin,lmax;

	enum Movement	movement;

	struct	ColorScale	colorscale;
	struct	Grid	G;
	struct	Particle	P;
	struct	WindowGeom	window;

	struct	BFaceList	*bface_root;

	Display *dpy;
	Window win;

//	unsigned int key_recved=0, isRunRecved=0;
	int firstR = 1;
	int *vars, *coms;

///	GLUquadricObj *cyl;

///	static GLfloat segcol[] = { 1.0, 0.0, 1.0, 1.0 };
/********** added by zhang, get existing data used to show as vector ****************/
void	getDispVarDefault(Domain *dom, int ivar, char *filename)
{
	char tecfilename[MAXLINLEN];
	DNode	*node_root=dom->dnode_root,
		*node=node_root;
	char line[100];
	size_t len;
	int	n;
	if(ivar<0) return;
	/*
	do//Count the nodes
	{	n++;
		node=node->next;
	}	while(node!=node_root);
	*/

	sprintf(tecfilename,"%s%dtec.dat",filename,ivar);
	std::ifstream f_in(tecfilename);
	if(!f_in) return ;
	if(node_root==NULL) return;
	//Read header
	f_in.getline(line,100);
	printf("\r%s\n", line);
	f_in.getline(line, 100);
	//Read data
	do//Output nodal variables
	{	double
			*x=node->x,
			*var=node->var;

		f_in>>n;
		f_in>>x[0];
		f_in>>x[1];
		f_in>>x[2];
		//for(int i=0; i<nvar; i++)
			int	i=ivar,
				varloc=dom->variable[i].loc,
				lenvar=dom->variable[i].dimension;
			double	*v=var+varloc;
		
			for(int j=0; j<lenvar; j++)
			{
				f_in>>v[j];
			}
		node=node->next;
	}	while(node!=node_root);
	f_in.close();
	return;
}
/*********/

int	query_extension(char* extName) 
{
	char *p = (char *) glGetString(GL_EXTENSIONS);
	char *end = p + strlen(p);
	while (p < end) 
	{
		int n = strcspn(p, " ");
		if ((strlen(extName) == n) && (strncmp(extName, p, n) == 0))
			return GL_TRUE;
		p += (n + 1);
	}
	return GL_FALSE;
}
void	initgui()
{
///	cyl = gluNewQuadric();
///	gluQuadricDrawStyle( cyl, GLU_FILL );
///	gluQuadricNormals( cyl, GLU_SMOOTH );

	showpar[dumpWindow] = 0;
	showpar[showRun] = 1;
	showpar[showAxes] = 1;
	showpar[showParticles] = 1;
	showpar[showSpheres] = 0;
	showpar[showNodes] = 0;
	showpar[showVariables] = 0;
//	showpar[showVectors] = 0;
	showpar[showBoundaryVertexes] = 0;
	showpar[showBoundaryVectors] = 0;
	showpar[showToolVertexes] = 0;
	showpar[showFrame] = 0;
	showpar[showBoundaryFaceCenters] = 0;
	showpar[showBoundaryFaces] = 0;
	showpar[showBoundaryGrid] = 1;
	showpar[showToolGrid] = 0;
	showpar[showCellCenters] = 0;
	showpar[showFaceCenters] = 0;
	showpar[showBonds] = 0;
	showpar[showGrid] = 0;
	showpar[showTool] = 1;
	finished=0;
	animation=0;
	lastz=1.0;
	SCOPY(configfile,"gui.cfg");//DEBUG: should be one cfg file for all
	wdtime=-1.0;
}
void	finish()
{	//Dump all unsaved data and close files
	
}
void	helpDisplay()
{
	printf("MOUSE CONTROL:\n");
	printf("left button & drag - rotate\n");
	printf("right button & drag - translate\n");
	printf("KEYBOARD CONTROL:\n");
	printf("?     - show help\n");
	printf("a     - show/hide axes\n");
	printf("c     - show/hide cell centers\n");
	printf("f     - show/hide face-centers\n");
	printf("G(g)  - show/hide (boundary) grid\n");
	printf("i     - initialize domains\n");
	printf("l     - load configuration from %s\n",configfile);
	printf("n     - show/hide nodes\n");
	printf("o     - toggle output\n");
	printf("q     - quit without saving\n");
	printf("r     - run\n");
	printf("s     - advance one iteration step\n");
	printf("t(T)  - toggle tool nodes (surface)\n");
	printf("+     - move forward along the arrow-axis\n");
	printf("-     - move backward along the arrow-axis\n");
	printf(":     - command mode\n esc: exit\n");
}
int	readparam
(
	char	*s,
	char	*param[],
	int	maxparam,
	char	*filename,
	int	val[]
)
{	int	i;
	for (i=0; i<maxparam; i++)
	{	int lparam;
		if(param[i]==NULL) continue;
		lparam=strlen(param[i]);
		if 
		(	lparam>0 &&
			(*(s+lparam)==' ' || *(s+lparam)=='=')
			&& memcmp(s,param[i],lparam)==0
		)	break;
	}
	if (i<maxparam)
	{	if ((s=strchr(s,'='))==NULL)
		{	fprintf
			(	stderr,
				"WARNING: wrong format at '%s' in file %s\n\t-> IGNORED\n",
				s,filename
			);
			return -1;
		}
		while(isspace(*++s));
		sscanf(s,"%d",&val[i]);
	}
	return i;
}
int	readparam
(
	char	*s,
	char	*param[],
	int	maxparam,
	char	*filename,
	float	val[]
)
{	int	i;
	for (i=0; i<maxparam; i++)
	{	int lparam=strlen(param[i]);
		if 
		(	(*(s+lparam)==' ' || *(s+lparam)=='=')
			&& memcmp(s,param[i],lparam)==0
		)	break;
	}
	if (i<maxparam)
	{	if ((s=strchr(s,'='))==NULL)
		{	fprintf
			(	stderr,
				"WARNING: wrong format at '%s' in file %s\n\t-> IGNORED\n",
				s,filename
			);
			return -1;
		}
		while(isspace(*++s));
		sscanf(s,"%g",&val[i]);
	}
	return i;
}
void	readconf()
{	int	i;
	char
		*colorname[maxcolor],
		*showparnam[maxshowpars],
		*pointparam[maxpointprm],
		*lineparam[maxlineprm],
		*axesprmnam[maxaxesprm],
		*toolprmnam[maxtoolprm],
		*filename=configfile,
		buf[MAXLINLEN],*s;
	float	a;
	FILE	*file;
	colorname [red      ]=strdup("color.red");
	colorname [green    ]=strdup("color.green");
	colorname [blue     ]=strdup("color.blue");
	colorname [skyblue  ]=strdup("color.skyblue");
	colorname [brown    ]=strdup("color.brown");
	colorname [magenta  ]=strdup("color.magenta");
	for (int i=magenta+1;i<maxcolor;i++)
	{	char	name[MAXLINLEN];
		sprintf(name,"color%d",i);
		colorname[i]=strdup(name);
	}
	for (int i=0;i<maxshowpars;i++)
		showparnam[i]=NULL;
	showparnam[showSpheres]=strdup("spheres"    );
	pointparam[variable ]=strdup("variable"     );
	pointparam[size     ]=strdup("size"         );
	pointparam[sizevar  ]=strdup("sizevar"      );
	pointparam[massvar  ]=strdup("massvar"      );
	pointparam[vmin     ]=strdup("value.min"    );
	pointparam[vmax     ]=strdup("value.max"    );
	lineparam [thickness]=strdup("thickness"    );
	lineparam [length   ]=strdup("length"       );
	axesprmnam[axeswidth]=strdup("width"        );
	axesprmnam[xaxislength]=strdup("length.x"   );
	axesprmnam[yaxislength]=strdup("length.y"   );
	axesprmnam[zaxislength]=strdup("length.z"   );
	axesprmnam[arrowidth  ]=strdup("arrow.width");
	axesprmnam[arrowheight]=strdup("arrow.length");
//#ifndef GRIDGEN
//	isoparnam[isosetvarf]=strdup("variable");
//	isoparnam[isodispvarf]=strdup("display.variable");
//	isoparnam[isolevel]=strdup("level");
//	isoparnam[isovarmin]=strdup("display.value.min");
//	isoparnam[isovarmax]=strdup("display.value.max");
//	isoparnam[isodisptype]=strdup("display.type");
//#endif
	for (i=0; i<maxcolor; i++)
	{	pointparam[i]=strdup(colorname[i]);
		lineparam [i]=strdup(colorname[i]);
		toolprmnam [i]=strdup(colorname[i]);
	}
	if ((file=fopen(filename,"r"))==NULL)
	{	fprintf
		(	stderr,
			"CAN'T OPEN CONFIGURATION FILE %s: using default settings\n",
			filename
		);return;
	}
	if(option.verbose)printf("Reading configuration from %s ...",filename);fflush(stdout);
	while (fgets(buf,MAXLINLEN,file)!=NULL)
	{
		for (s=buf;isspace(*s);s++);
		if (memcmp(s,"Translation.step",16)==0)
		{	if ((s=strchr(s,'='))==NULL)
			{
				fprintf(stderr,"MISSING '=' in Step difinition in %s\n",filename);
				continue;
			}
			while(isspace(*++s));
			sscanf(s,"%g",&step);
			dx=step*lx;
			dy=step*ly;
			dz=step*lz;
		}
		else
		if (memcmp(s,"Display.",8)==0)
		{
			if (readparam(s+8,showparnam,maxshowpars,filename,showpar)==-1)
				continue;
		}
		else
		if (memcmp(s,"Vector.",7)==0)
		{
			if (readparam(s+7,lineparam,maxlineprm,filename,vecval)==-1)
				continue;
		}
		else
		if (memcmp(s,"Grid.",5)==0)
		{	s+=5;
			if (memcmp(s,"node.",5)==0)
			{	if 
				(	(i=readparam(s+5,pointparam,maxpointprm,filename,G.node))==-1
				)	continue;
			}
			else
			if (memcmp(s,"line.",5)==0)
			{	if 
				(	(i=readparam(s+5,lineparam,maxlineprm,filename,G.line))==-1
				)	continue;
			}
		}
		else
		if (memcmp(s,"Particles.",10)==0)
		{	if 
				(	(i=readparam(s+10,pointparam,maxpointprm,filename,P.disp))==-1
				)	continue;
		}
		else
		if (memcmp(s,"Axes.",5)==0)
		{	if 
			(	(i=readparam(s+5,axesprmnam,maxaxesprm,filename,axes))==-1
			)	continue;
		}
		else
		if (memcmp(s,"Tool.",5)==0)
		{	if 
			(	(i=readparam(s+5,toolprmnam,maxtoolprm,filename,toolval))==-1
			)	continue;
		}
	}
	fclose(file);
	if(option.verbose)printf(" done\n");fflush(stdout);
}
void	refresh(int idomain)
{	Domain	*dom=domain_root+idomain;
	if
	(	showpar[showBoundaryFaces]||showpar[showBoundaryGrid]
	)
	{
		if(dom->getNoVariables(boundary_faces)==0)
		{	dom->deleteRing(bface_root);
			dom->createBoundaryFaceList(bface_root);
		}
		else bface_root=dom->bface_root;
	}
}
void	initdisp()
{
	if(option.debug){printf("Initializing display\n");fflush(stdout);}
	idomain=0;
	rgbcolor[red    ][red      ]=1.0;
	rgbcolor[red    ][green    ]=0.0;
	rgbcolor[red    ][blue     ]=0.0;
	rgbcolor[green  ][red      ]=0.0;
	rgbcolor[green  ][green    ]=1.0;
	rgbcolor[green  ][blue     ]=0.0;
	rgbcolor[blue   ][red      ]=0.0;
	rgbcolor[blue   ][green    ]=0.0;
	rgbcolor[blue   ][blue     ]=1.0;
	rgbcolor[skyblue][red      ]=0.0;
	rgbcolor[skyblue][green    ]=1.0;
	rgbcolor[skyblue][blue     ]=1.0;
	rgbcolor[brown  ][red      ]=1.0;
	rgbcolor[brown  ][green    ]=1.0;
	rgbcolor[brown  ][blue     ]=0.0;
	rgbcolor[magenta][red      ]=1.0;
	rgbcolor[magenta][green    ]=0.0;
	rgbcolor[magenta][blue     ]=1.0;
	colorscale.relative=0;
	if (maxcolor>6)
	{	int
			ncolor=(int)pow((double)maxcolor,1./3.),
			icolor=5;
		for (int ired  =1; ired  <ncolor; ired++   )
		for (int igreen=1; igreen<ncolor; igreen++ )
		for (int iblue =1; iblue <ncolor; iblue++  )
		{	if (++icolor>=maxcolor) break;
			rgbcolor[icolor][ired  ]=(float)ired  /(float)ncolor;
			rgbcolor[icolor][igreen]=(float)igreen/(float)ncolor;
			rgbcolor[icolor][iblue ]=(float)iblue /(float)ncolor;			
		}
	}
	vecval[length   ]=.2;
	vecval[thickness]=1.5;
	vecval[red      ]=0.5;
	vecval[green    ]=1.0;
	vecval[blue     ]=1.0;
	axes[axeswidth  ]=5.0;
	axes[xaxislength]=0.1;
	axes[yaxislength]=0.1;
	axes[zaxislength]=0.1;
	axes[arrowheight]=0.1;
	axes[arrowidth  ]=0.03;
//	G.ivar=-1;
	G.node[size ]=1.0;
	G.node[red  ]=1.;
	G.node[green]=1.;
	G.node[blue ]=0.;
	G.line[thickness]=1.0;
	G.line[red      ]=1.0;
	G.line[green    ]=0.0;
	G.line[blue     ]=0.0;
//	G.isoline[thickness]=1.0;
//	G.isoline[red      ]=0.0;
//	G.isoline[green    ]=1.0;
//	G.isoline[blue     ]=1.0;
#ifndef GRIDGEN
//	G.isovar=-1;
//	G.isodispvar=-1;
//	P.iv=-1;
	P.disp[size ]=3.0;
	P.disp[red  ]=0;
	P.disp[green]=1.0;
	P.disp[blue ]=0.0;
	toolval[red  ]=0.5;
	toolval[green]=0.5;
	toolval[blue ]=0.5;
#endif
	readconf();
	for (int i=0; i<DIM; i++)
	{	xmin[i]=LARGE; xmax[i]=-LARGE;}
	for 
	(	Domain	*d=domain_root; d-domain_root<ndomains; d++
	)
	{	int	idomain=d-domain_root;
		double	x0[DIM],x1[DIM];
		//Get global grid limits
		d->getGridLimits(x0,x1);
		for (int i=0; i<DIM; i++)
		{	if (xmin[i]>x0[i])xmin[i]=x0[i];
			if (xmax[i]<x1[i])xmax[i]=x1[i];
		}
		if(option.debug)
		{	printf
			(	"Grid limits: xmin={%g,%g,%g}, xmax={%g,%g,%g}\n",
				xmin[0],xmin[1],xmin[2],xmax[0],xmax[1],xmax[2]
			); fflush(stdout);
		}
		d->setDisplay();
		for (int eltype=0; eltype<maxelements; eltype++)
			d->setElementColor(eltype,rgbcolor[(idomain+eltype)%maxcolor]);
		d->setElementColor(nodes,rgbcolor[green]);
		d->setElementColor(edges,rgbcolor[red]);
		d->setElementColor(faces,rgbcolor[brown]);
		d->setElementColor(cells,rgbcolor[skyblue]);
	}
	lmin= LARGE;
	lmax=-LARGE;
	lx=xmax[0]-xmin[0];
	ly=xmax[1]-xmin[1];
	lz=xmax[2]-xmin[2];
	if (lx>lmax) lmax=lx;
	if (ly>lmax) lmax=ly;
	if (lz>lmax) lmax=lz;
	if (lx<lmin) lmin=lx;
	if (ly<lmin) lmin=ly;
	if (lz<lmin) lmin=lz;
	xo[0]=xmin[0]+.5*lx;
	xo[1]=xmin[1]+.5*ly;
	xo[2]=xmin[2]+.5*lz;
	dx=step*lx;
	dy=step*ly;
	dz=step*lz;
	if (option.debug|option.verbose)
		printf("Origin of coordinate axes: %g\t%g\t%g\n",xo[0],xo[1],xo[2]);
	{	float	window_size;
		window_size=1.5*WINDOW_SIZE/sqrt(lx*lx+ly*ly);
		window.width=(int)(window_size*lx); 
		if (window.width>MAX_WINDOW_WIDTH)window.width=MAX_WINDOW_WIDTH; 
		window.height=(int)(window_size*ly);
		if (window.height>MAX_WINDOW_HEIGHT)window.height=MAX_WINDOW_HEIGHT; 
	}
}
void	InitMaterials(void)
{
//		static float ambient[] =
//		{0.1, 0.1, 0.1, 1.0};
//		static float diffuse[] =
//		{0.5, 1.0, 1.0, 1.0};
//		static float position0[] =
//		{0.0, 0.0, 20.0, 0.0};
//		static float position1[] =
//		{0.0, 0.0, -20.0, 0.0};
//		static float front_mat_shininess[] =
//		{60.0};
//		static float front_mat_specular[] =
//		{0.2, 0.2, 0.2, 1.0};
//		static float front_mat_diffuse[] =
//		{0.5, 0.28, 0.38, 1.0};
//		static float lmodel_ambient[] =
//		{1.0, 1.0, 1.0, 1.0};
//		static float lmodel_twoside[] =
//		{GL_FALSE};
//	
//		glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
//		glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
//		glLightfv(GL_LIGHT0, GL_POSITION, position0);
//		glEnable(GL_LIGHT0);
//	
//		glLightfv(GL_LIGHT1, GL_AMBIENT, ambient);
//		glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse);
//		glLightfv(GL_LIGHT1, GL_POSITION, position1);
//		glEnable(GL_LIGHT1);
//	
//		glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
//		glLightModelfv(GL_LIGHT_MODEL_TWO_SIDE, lmodel_twoside);
//		glEnable(GL_LIGHTING);
//	
//		glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, front_mat_shininess);
//		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, front_mat_specular);
//		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, front_mat_diffuse);


//		GLfloat diffuseMaterial[4] = { 0.5, 0.5, 0.5, 1.0 };
//	
//	
//		 GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
//		 GLfloat light_position[] = { 1.0, 1.0, 1.0, 0.0 };
//	
//		 glClearColor (0.0, 0.0, 0.0, 0.0);
//		 glShadeModel (GL_SMOOTH);
//		 glEnable(GL_DEPTH_TEST);
//		 glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuseMaterial);
//		 glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
//		 glMaterialf(GL_FRONT, GL_SHININESS, 25.0);
//		 glLightfv(GL_LIGHT0, GL_POSITION, light_position);
//		 glEnable(GL_LIGHTING);
//		 glEnable(GL_LIGHT0);
//	
//		 glColorMaterial(GL_FRONT, GL_DIFFUSE);
//		 glEnable(GL_COLOR_MATERIAL);


    GLfloat ambient[] = { 0.0, 0.0, 0.0, 1.0 };
    GLfloat diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
    GLfloat position[] = { 0.0, 3.0, 3.0, 0.0 };

    GLfloat lmodel_ambient[] = { 0.2, 0.2, 0.2, 1.0 };
    GLfloat local_view[] = { 0.0 };

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

    glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
    glLightfv(GL_LIGHT0, GL_POSITION, position);

    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
    glLightModelfv(GL_LIGHT_MODEL_LOCAL_VIEWER, local_view);

    glFrontFace (GL_CW);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_AUTO_NORMAL);
    glEnable(GL_NORMALIZE);

    glMatrixMode (GL_MODELVIEW);
    glLoadIdentity ();

    glClearColor(0.0, 0.0, 0.0, 0.0);
    glClearAccum(0.0, 0.0, 0.0, 0.0);

}
/*  Initialize material property and depth buffer.
 */
void	myinit(int argc, char *argv[])
{
	GLfloat mat[4] = { 0.8, 0.7, 0.5, 1.0 },
		shine=0.6;
	if(option.debug)
	{	printf("myinit: argc=%d",argc);FLUSH;
		for (int i=0; i<argc; i++)
		printf(" '%s'",argv[i]);
		printf("\n");FLUSH;
	}
//	GLfloat fogColor[4] = {0.0, 0.0, 0.0, 1.0};
//	  GLfloat mat_diffuse[] = { 0.7, 0.7, 0.7, 1.0 };
//	  GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
//	  GLfloat mat_shininess[] = { 100.0 };
//	GLfloat light_ambient[] = { 2.0, 1.0, 0.0, 1.0 };
//	GLfloat light_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
//	GLfloat light_specular[] = { 1.0, 2.0, 0.0, 1.0 };
//	GLfloat light_position[] = { 0.0, 100.0, 1000.0, 1.0 };

//	glClearColor (0.0, 0.0, 0.0, 1.0);

//	glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
//	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
//	glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);

//	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
//	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
//	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
//	glLightfv(GL_LIGHT0, GL_POSITION, light_position);

//	glEnable(GL_LIGHTING);
//	glEnable(GL_LIGHT0);

//	glDepthFunc(GL_LESS);
//	glEnable(GL_DEPTH_TEST);
//	glEnable(GL_AUTO_NORMAL);
//	glEnable(GL_NORMALIZE);
//	glShadeModel (GL_SMOOTH);
//	glShadeModel (GL_FLAT);

// Lines
//	glEnable (GL_LINE_SMOOTH);
//	glEnable (GL_BLEND);
//	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//	glHint (GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
//	glLineWidth (G.line[thickness]);

//	glClearColor(0.0, 0.0, 0.0, 0.0);
//	glClearAccum(0.0, 0.0, 0.0, 0.0);

//	  theNurb = gluNewNurbsRenderer();
//	  gluNurbsProperty(theNurb, GLU_SAMPLING_TOLERANCE, 25.0);
//	  gluNurbsProperty(theNurb, GLU_DISPLAY_MODE, GLU_FILL);
//	  glMatrixMode(GL_MODELVIEW);
//	  glLoadIdentity();

//	glEnable(GL_FOG);
//	glFogi (GL_FOG_MODE, GL_LINEAR);
//	glHint (GL_FOG_HINT, GL_NICEST);  /*  per pixel   */
//	glFogf (GL_FOG_START, 4*lz  );
//	glFogf (GL_FOG_END, 10.*lmax);
//	glFogfv (GL_FOG_COLOR, fogColor);

	InitMaterials();
		//glEnable(GL_LIGHTING);
		//glDisable(GL_LIGHTING);
    glMaterialfv (GL_FRONT, GL_AMBIENT, mat);
    glMaterialfv (GL_FRONT, GL_DIFFUSE, mat);
    glMaterialfv (GL_FRONT, GL_SPECULAR, mat);
    glMaterialf (GL_FRONT, GL_SHININESS, shine*128.0);

	glTranslatef(-xo[0], -xo[1], -xo[2]-3*lz);
}
void	showGridElements
(	int	ne,//number of elements to show
	double	*X,//coordinates
	float	color[]
)
{	int	i;
	glPointSize(G.node[size]);
	glColor3f(color[red], color[green], color[blue]);
	glDisable(GL_LIGHTING);
	glBegin(GL_POINTS);
	for (i=0; i<ne; i++)
	{//	float	r,g,b;
		//double	*x=element[i].x;
		double	*x=X+DIM*i;
		glVertex3f((GLfloat)x[0],(GLfloat)x[1],(GLfloat)x[2]);
	}
	glEnd();
	glEnable(GL_LIGHTING);
}
void	showVector
(
	double	*x,//coordinates
	double	*v //vector
)
{
	if (vecval[length]*(fabs(v[0])+fabs(v[1])+fabs(v[2]))>1.e6)
	{	fprintf
		(	stderr,
			"VALUE AT (%g,%g,%g) TOO LARGE (not displayed)\n",
			x[0],x[1],x[2]
		);
		return;
	}
	glBegin(GL_LINES);
	glVertex3f
	(	(GLfloat)x[0],
		(GLfloat)x[1],
		(GLfloat)x[2]
	);
	glVertex3f
	(	(GLfloat)(x[0]+vecval[length]*v[0]),
		(GLfloat)(x[1]+vecval[length]*v[1]),
		(GLfloat)(x[2]+vecval[length]*v[2])
	);
	glEnd();
}
template <class Element>
void	showDVectors
(	int loc,
	Element *root
)
{
	if (root!=NULL)
	{	Element	*e=root;
		do
		{	double
				*x=e->x,
				*v=e->var+loc;
			showVector(x,v);
//if((x[0]==0)&&(x[1]==0)) printf("x: %g, %g, %g; v: % g, %g, %g \n", x[0],x[1],x[2], v[0],v[1],v[2]);
			e=e->next;
		}	while(e!=root);
	}
}
void	showVectors
(
	int	ivar,
	Domain	*dom
)
{
	Variable	*var=dom->variable+ivar;
	int
		type=var->type,
		size=var->size;
	double
		*X=dom->coordinates[type].val,//coordinates
		*x,*v,*u;
	glLineWidth (vecval[thickness]);
	glDisable(GL_LIGHTING);
	glColor3f(vecval[red], vecval[green], vecval[blue]);
	u=var->val; //G.V+G.var[G.ivec].i*G.n*DIM;
	if (type==points)
	{	int	n=dom->getNoPoints();
		struct Point	*p,*origin;
		dom->getPoints(P.m,P.n,P.origin,P.first,P.last,P.X);
		origin=P.origin;
		p=P.first;
		for (int i=0; i<n; i++)
		{	int	ip=p-origin,jp=DIM*ip;
			x=X+jp;v=u+jp;
			showVector(x,v);
			p=p->next;
		}
	}
	else
	if (dom->type==dynamic)
	{	int	loc=var->loc;
		switch(type)
		{	case nodes:
			showDVectors(loc,dom->dnode_root);
			break;
			case cells:
	///		showDVectors(loc,dom->dcell_root);
	///		break;
			default:
			fprintf(stderr,"Can't display vector variable type %d",type);
		}
	}
	else
		for (x=X,v=u;x-X<size*DIM;x+=DIM,v+=DIM)
			showVector(x,v);
	glEnable(GL_LIGHTING);
}
void	showPoint
(
	double	vmn,
	double	vmx,
	double	val,
	double	*x
)
{
		float	r,g,b,
			vav=.5*(vmn+vmx),
			alow =PI/(2.*(vav-vmn)),
			ahigh=PI/(2.*(vmx-vav)),
			a;
		if (val<vmn)
		{	r=g=0.0; b=1.0; }
		else
		if (val>vmx)
		{	r=1.0; g=b=0.0; }
		else
		if (val<vav)
		{	a=(float)(val-vmn)*alow;
			b=cos(a); g=sin(a);
			r=0.0;
		}
		else
		{	a=(float)(val-vav)*ahigh;
			g=cos(a); r=sin(a);
			b=0.0;
		}
		glColor3f(r,g,b);
		glVertex3f((GLfloat)x[0],(GLfloat)x[1],(GLfloat)x[2]);
}
void	showSphere
(
	double	vmn,
	double	vmx,
	double	val,
	double	rad,
	double	*x
)
{
		GLfloat mat[4],/// = { 0.6, 0.8, 0.3, 1.0 },
				shine=0.6;
		float	r,g,b,
			vav=.5*(vmn+vmx),
			alow =PI/(2.*(vav-vmn)),
			ahigh=PI/(2.*(vmx-vav)),
			a;
		if (val<vmn)
		{	r=g=0.0; b=1.0; }
		else
		if (val>vmx)
		{	r=1.0; g=b=0.0; }
		else
		if (val<vav)
		{	a=(float)(val-vmn)*alow;
			b=cos(a); g=sin(a);
			if(g<0.0)g=0.0;
			if(b<0.0)b=0.0;
			r=0.0;
		}
		else
		{	a=(float)(val-vav)*ahigh;
			g=cos(a); r=sin(a);
			if(g<0.0)g=0.0;
			if(r<0.0)r=0.0;
			b=0.0;
		}
		mat[0]=r;mat[1]=g;mat[2]=b;mat[3]=1.0;
		glEnable(GL_LIGHTING);
		glMaterialfv (GL_FRONT, GL_AMBIENT, mat);
		glMaterialfv (GL_FRONT, GL_DIFFUSE, mat);
		glMaterialfv (GL_FRONT, GL_SPECULAR, mat);
		glMaterialf (GL_FRONT, GL_SHININESS, shine*128.0);
		glTranslatef (x[0], x[1], x[2]); 
		glutSolidSphere (rad, 15, 15);
		glTranslatef (-x[0],-x[1],-x[2]); 
}
template <class Element>
void	showDScalar
(	int loc,
	Element *root,
	double	vmn,
	double	vmx
)
{	if (root==NULL)return;
	{	Element	*e=root;
		glBegin(GL_POINTS);
		do
		{	double value=e->var[loc];
			showPoint(vmn,vmx,value,e->x);
			e=e->next;
		}	while(e!=root);
		glEnd();
	}
}
void	showScalar
(
	int	ivar,
	int	icomp,
	double	vmn,
	double	vmx,
	Domain	*dom
)
{
	Variable	*var=dom->variable+ivar;
	char	*name=var->name;
	int
		type=var->type,
		dim=var->dimension,
		size=var->size,
		jcomp=icomp;
	double
		*X=dom->coordinates[type].val,
		*value=var->val;
	if (icomp>DIM)
	{	fprintf
		(	stderr,"Can't display variable component %d for variable %s\n",
			icomp,name
		);
		return;
	}
	if (dim==DIM&&icomp<0)
	{	fprintf
		(	stderr,"A vector component should first be selected for variable %s\n",
			name
		);
		return;
	}
	if (icomp>=dim)
	{	fprintf
		(	stderr,
			"Can't display component %d of variable %s of dim=%d\n",
			icomp,name,dim
		);
		return;
	}
	if (icomp<0) jcomp=0;
	if(type==points)
	{	int	n=dom->getNoPoints();
		struct Point	*p,*origin;
		dom->getPoints(P.m,P.n,P.origin,P.first,P.last,P.X);
		origin=P.origin;
		p=P.first;
		if (showpar[showSpheres])
		{	GLfloat mat[4] = { 0.6, 0.8, 0.3, 1.0 },
				shine=0.6;
			glEnable(GL_LIGHTING);
			glMaterialfv (GL_FRONT, GL_AMBIENT, mat);
			glMaterialfv (GL_FRONT, GL_DIFFUSE, mat);
			glMaterialfv (GL_FRONT, GL_SPECULAR, mat);
			glMaterialf (GL_FRONT, GL_SHININESS, shine*128.0);
			for (int i=0; i<n; i++)
			{	int	ip=p-origin;
				double val=value[ip*dim+jcomp],
					rad=(fabs(val)-vmn)/(vmx-vmn),///=0.2
					*x=X+DIM*ip;
				glTranslatef (x[0], x[1], x[2]); 
				if (val>=0.0)
					glutSolidSphere (rad, 15, 15);
				else
				{	//color different if val<0
					mat[0]=.8;mat[1]=0.7;mat[2]=0.5;
					glMaterialfv (GL_FRONT, GL_AMBIENT, mat);
					glMaterialfv (GL_FRONT, GL_DIFFUSE, mat);
					glMaterialfv (GL_FRONT, GL_SPECULAR, mat);
					glutSolidSphere (rad, 15, 15);
					mat[0] = 0.6;
					mat[1] = 0.8;
					mat[2] = 0.3;
					glMaterialfv (GL_FRONT, GL_AMBIENT, mat);
					glMaterialfv (GL_FRONT, GL_DIFFUSE, mat);
					glMaterialfv (GL_FRONT, GL_SPECULAR, mat);
				}
				glTranslatef (-x[0],-x[1],-x[2]); 
				//showPoint(vmn,vmx,rad,x);
				p=p->next;
			}
		}
		else
		{	glBegin(GL_POINTS);
			for (int i=0; i<n; i++)
			{	int	ip=p-origin;
				showPoint(vmn,vmx,value[ip*dim+jcomp],X+DIM*ip);
				p=p->next;
			}
			glEnd();
		}
	}
	else
	if (dom->type==dynamic)
	{	int	loc=var->loc;
		switch(type)
		{	case nodes:
				showDScalar(loc,dom->dnode_root,vmn,vmx);
				break;
			case cells:
			///	showDScalar(loc,dom->dcell_root,vmn,vmx);
			///	break;
			default:
			fprintf(stderr,"Can't display scalar variable type %d",type);
		}
	}
	else
	{
		glBegin(GL_POINTS);
		for (int i=0; i<size; i++)
			showPoint(vmn,vmx,value[i*dim+jcomp],X+DIM*i);
		glEnd();
	 }
}
void	displayVariables(Domain *dom)
{	int
		ivar=dom->getDisplayVariable(),
		icomp=dom->getDisplayVariableComp();
	class Variable	*var;
	if (ivar>=0)
	{	var=dom->variable+ivar;
		if (var->type==points)
			glPointSize(P.disp[size]);
		else
			glPointSize(G.node[size]);
		glDisable(GL_LIGHTING);
		if(icomp<0) showVectors(ivar,dom);
		else
		{	double	vmn,vmx;
			if (colorscale.relative)
			{	var->computeLimits();
				var->getLimits(icomp,vmn,vmx);
				if(var->type==points)
				{	P.disp[vmin]=vmn;
					P.disp[vmax]=vmx;
				}
				else
				{	G.node[vmin]=vmn;
					G.node[vmax]=vmx;
				}
			}
			else
			{	if(var->type==points)
				{	vmn=P.disp[vmin];
					vmx=P.disp[vmax];
				}
				else
				{	vmn=G.node[vmin];
					vmx=G.node[vmax];
				}
			}
			if(vmx>vmn)
				showScalar
				(	ivar,icomp,
					vmn,vmx,
					dom
				);
			else
			{	if (var->dimension>1)
				fprintf
				(	stderr,
					"Component %d of ",icomp
				);
				if (vmx==vmn)
				fprintf
				(	stderr,
					"%s has a constant value of %g in the whole domain %s\n",
					var->name,vmn,dom->name
				);
				else
				{	fprintf
					(	stderr,
						"%s limits are in a wrong range: min=%g, max=%g\n",
						var->name,G.node[vmin],G.node[vmax]
					);BUG("Internal inconsistency");
				}
			}
		}
		glEnable(GL_LIGHTING);
	}
}
void	displayAxes()
{
	GLfloat mat[4] = { 0.6, 0.8, 0.3, 1.0 },
		shine=0.6;
	mat[0]=.8;mat[1]=0.7;mat[2]=0.5;
	glMaterialfv (GL_FRONT, GL_AMBIENT, mat);
	glMaterialfv (GL_FRONT, GL_DIFFUSE, mat);
	glMaterialfv (GL_FRONT, GL_SPECULAR, mat);
	glMaterialf (GL_FRONT, GL_SHININESS, shine*128.0);
	if (!showpar[showAxes]) return;
	glEnable(GL_LIGHTING);
	glLineWidth (axes[axeswidth]);
	glBegin (GL_LINES);
	glVertex3f(xo[0],xo[1],xo[2]);
	glVertex3f(xo[0]+axes[xaxislength]*lx,xo[1],xo[2]);
	glVertex3f(xo[0],xo[1],xo[2]);
	glVertex3f(xo[0],xo[1]+axes[yaxislength]*ly,xo[2]);
	glVertex3f(xo[0],xo[1],xo[2]);
	glVertex3f(xo[0],xo[1],xo[2]+axes[zaxislength]*lz);
	glEnd ();
	glTranslatef (xo[0], xo[1], xo[2]+axes[zaxislength]*lz);
	glutSolidCone (axes[arrowidth]*lz, axes[arrowheight]*lz, 4, 4);
	glTranslatef (-xo[0], -xo[1], -(xo[2]+axes[zaxislength]*lz));
}
void	displayParticles(Domain *dom)
{	float	color[3];
	dom->getElementColor(points,color);
	dom->getPoints(P.m,P.n,P.origin,P.first,P.last,P.X);
	if(P.n>0) 
	{	int icolorvar=(int)P.disp[variable];
		struct Point	*p,*origin;
		glColor3f(P.disp[red], P.disp[green], P.disp[blue]);
		origin=P.origin;
		p=P.first;
		if (showpar[showSpheres])
		{	GLfloat mat[4] = { 0.6, 0.8, 0.3, 1.0 },
				shine=0.6;
			double 
				psize=P.disp[size],
				vmn=P.disp[vmin],
				vmx=P.disp[vmax],
				dv=vmx-vmn,
				///rad=(fabs(psize)-vmn)/(vmx-vmn);
				rad;//particle radius
			mat[0]=P.disp[red];mat[1]=P.disp[green];mat[2]=P.disp[blue];
			glEnable(GL_LIGHTING);
			glMaterialfv (GL_FRONT, GL_AMBIENT, mat);
			glMaterialfv (GL_FRONT, GL_DIFFUSE, mat);
			glMaterialfv (GL_FRONT, GL_SPECULAR, mat);
			glMaterialf (GL_FRONT, GL_SHININESS, shine*128.0);
			if(psize>0)//scale radius by psize
			{	rad=psize*(rad-vmn)/dv;//calibrate the radius
				do 
				{	int ip=(int)(p-origin);
					double	*x=P.X+ip*DIM;
					if(icolorvar>=0)
					{	double colorvar=dom->variable[icolorvar].val[ip];
						showSphere(vmn,vmx,colorvar,rad,x);
					}
					else
					{	glTranslatef (x[0], x[1], x[2]); 
						glutSolidSphere (rad, 15, 15);
						glTranslatef (-x[0],-x[1],-x[2]); 
					}
					p=p->next;
				}	while(p!=P.last->next);
			}
			else //get size from a variable
			{	if(P.disp[sizevar]>=0)
				{	int ivar=(int)P.disp[sizevar];
					do 
					{	int ip=(int)(p-origin);
						double	val=dom->variable[ivar].val[ip],
							*x=P.X+ip*DIM;
						rad=(-psize)*(fabs(val)-vmn)/dv;//calibrate the radius
						if(icolorvar>=0)
						{	double colorvar=dom->variable[icolorvar].val[ip];
							showSphere(vmn,vmx,colorvar,rad,x);
						}
						else
						{	glTranslatef (x[0], x[1], x[2]); 
							glutSolidSphere (rad, 15, 15);
							glTranslatef (-x[0],-x[1],-x[2]); 
						}
						p=p->next;
					}	while(p!=P.last->next);
				}
				else
				if(P.disp[massvar]>=0)
				{	int ivar=(int)P.disp[massvar];
					const double pi=4.0*atan(1.);
					double	density=1.0,///DDD
						factor=0.75/(density*pi);
					do 
					{	int ip=(int)(p-origin);
						double	val=dom->variable[ivar].val[ip],
							*x=P.X+ip*DIM;
						rad=pow(factor*fabs(val),1.0/3.0);
						rad=(-psize)*(rad-vmn)/dv;//calibrate the radius
						if(icolorvar>=0)
						{	double colorvar=dom->variable[icolorvar].val[ip];
							showSphere(vmn,vmx,colorvar,rad,x);
						}
						else
						{	glTranslatef (x[0], x[1], x[2]); 
							glutSolidSphere (rad, 15, 15);
							glTranslatef (-x[0],-x[1],-x[2]); 
						}
						p=p->next;
					}	while(p!=P.last->next);
				}
				else
				{	fprintf
					(	stderr,
						"WARNING: Either Particles.size or Particles.sizevar or Particles.massvar should be specified in gui.cfg to display particle sizes.\n"
					);
				}
			}
		}
		else
		{
			glDisable(GL_LIGHTING);
			glPointSize(P.disp[size]);
			glBegin(GL_POINTS);
			do 
			{	double	*x=P.X+(int)(p-origin)*DIM;
				glVertex3f((GLfloat)x[0],(GLfloat)x[1],(GLfloat)x[2]);
				p=p->next;
			}	while(p!=P.last->next);
			glEnd();
			glEnable(GL_LIGHTING);
		}
	}
}
void	displayDynamicDomain(int idomain)
{
	Domain	*dom=domain_root+idomain;
	Tool	*tool=dom->tool;
	refresh(idomain);
	if (showpar[showToolVertexes]&&tool!=NULL)
	{
		GLfloat mat[4] = { 0.6, 0.8, 0.3, 1.0 },
			shine=0.6;
		DNodeList
			*first_tool_boundary_dnode=dom->tool->first_tool_boundary_dnode,
			*first_tool_dnode=dom->tool->first_tool_dnode;
		//SHOW BOUNDARY NODES
		glEnable(GL_LIGHTING);
		glMaterialfv (GL_FRONT, GL_AMBIENT, mat);
		glMaterialfv (GL_FRONT, GL_DIFFUSE, mat);
		glMaterialfv (GL_FRONT, GL_SPECULAR, mat);
		glMaterialf (GL_FRONT, GL_SHININESS, shine*128.0);
		if (first_tool_boundary_dnode!=NULL)
		for 
		(	DNodeList *n=first_tool_boundary_dnode->next;
			;	n=n->next
		) 
		{	double	*x=n->node->x;
			glTranslatef (x[0], x[1], x[2]); 
			glutSolidSphere (0.2, 15, 15);
			glTranslatef (-x[0],-x[1],-x[2]); 
			if (n==first_tool_boundary_dnode) break;
		}
		//SHOW NEAR-BOUNDARY NODES
		mat[0]=.8;mat[1]=0.7;mat[2]=0.5;
		glMaterialfv (GL_FRONT, GL_AMBIENT, mat);
		glMaterialfv (GL_FRONT, GL_DIFFUSE, mat);
		glMaterialfv (GL_FRONT, GL_SPECULAR, mat);
		glMaterialf (GL_FRONT, GL_SHININESS, shine*128.0);
		if (first_tool_dnode!=NULL)
		for 
		(	DNodeList *n=first_tool_dnode->next;
			;
			n=n->next
		) 
		{	double	*x=n->node->x;
			glTranslatef (x[0], x[1], x[2]); 
			glutSolidSphere (0.2, 7, 7);
			glTranslatef (-x[0],-x[1],-x[2]); 
			if (n==first_tool_dnode) break;
		}
		//RETURN TO THE ORIGINAL COLORS
		mat[0]=.8;mat[1]=0.7;mat[2]=0.5;
		glMaterialfv (GL_FRONT, GL_AMBIENT, mat);
		glMaterialfv (GL_FRONT, GL_DIFFUSE, mat);
		glMaterialfv (GL_FRONT, GL_SPECULAR, mat);
		glMaterialf (GL_FRONT, GL_SHININESS, shine*128.0);
		glEnable(GL_LIGHTING);
	}
	if (showpar[showToolGrid]&&tool!=NULL)
	{
		DCellList
			*first_tool_boundary_dcell=dom->tool->first_tool_boundary_dcell,
			*first_relax_dcell=dom->tool->first_relax_dcell;
		glDisable(GL_LIGHTING);
		glColor3f(vecval[red], vecval[green], vecval[blue]);
		glLineWidth (vecval[thickness]);
		if (first_tool_boundary_dcell!=NULL)
		for 
		(	DCellList *c=first_tool_boundary_dcell->next;
			; c=c->next
		) 
		{	DNode	**vert=c->cell->vert;
			double	*x;
			glBegin(GL_LINES);
			for (int i=0; i<Nv; i++)
			{	x=vert[i]->x;
				
				glVertex3f(x[0],x[1],x[2]);
				x=vert[(i+1)%Nv]->x;
				glVertex3f(x[0],x[1],x[2]);
			}
			x=vert[0]->x;
			glVertex3f(x[0],x[1],x[2]);
			x=vert[2]->x;
			glVertex3f(x[0],x[1],x[2]);
			x=vert[1]->x;
			glVertex3f(x[0],x[1],x[2]);
			x=vert[3]->x;
			glVertex3f(x[0],x[1],x[2]);
			glEnd();
			if (c==first_tool_boundary_dcell) break;
		}
		if (first_relax_dcell!=NULL)
		for 
		(	DCellList *c=first_relax_dcell->next;
			; c=c->next
		) 
		{	DNode	**vert=c->cell->vert;
			double	*x;
			glBegin(GL_LINES);
			for (int i=0; i<Nv; i++)
			{	x=vert[i]->x;
				glVertex3f(x[0],x[1],x[2]);
				x=vert[(i+1)%Nv]->x;
				glVertex3f(x[0],x[1],x[2]);
			}
			x=vert[0]->x;
			glVertex3f(x[0],x[1],x[2]);
			x=vert[2]->x;
			glVertex3f(x[0],x[1],x[2]);
			x=vert[1]->x;
			glVertex3f(x[0],x[1],x[2]);
			x=vert[3]->x;
			glVertex3f(x[0],x[1],x[2]);
			glEnd();
			if (c==first_relax_dcell) break;
		}
		glEnable(GL_LIGHTING);
	}
	if (showpar[showBoundaryGrid]&&bface_root!=NULL)
	{
		float	color[3];
		///BFaceList	*bface_root=dom->bface_root;	
		glDisable(GL_LIGHTING);
		glLineWidth (G.line[thickness]);
		dom->getElementColor(edges,color);
		glColor3f(color[red], color[green], color[blue]);
		for
		(	BFaceList	*c=bface_root->next;
			; c=c->next
		)
		{	DCell	*cell=c->cell;
			DNode	**vert=cell->vert;
			int	iv=c->iface;
			//if (cell->facetype[iv]==boundary)
			//{	
				glBegin(GL_LINE_LOOP);
				for (int j=0; j<Nfv; j++)
				{	double *x=vert[(iv+j+1)%Nv]->x;
					glVertex3f((GLfloat)x[0],(GLfloat)x[1],(GLfloat)x[2]);
				}
				glEnd();
			//}
			if (c==bface_root) break;
		}
	}
	if (showpar[showGrid])
	{
		float	color[3];
		//int showboundarygrid=showpar[showBoundaryGrid];
		DCell *first=dom->dcell_root;
		glDisable(GL_LIGHTING);
		glLineWidth (G.line[thickness]);
		dom->getElementColor(edges,color);
		glColor3f(color[red], color[green], color[blue]);
		if (first!=NULL)
		for 
		(	DCell *c=first->next;
			;c=c->next
		) 
		{	DNode	**vert=c->vert;
			DCell	**neib=(DCell**)c->neib;
			double	*x;
			//ElementState	*facetype=c->facetype;
			//if (showboundarygrid)
			//{	for (int i=0; i<Nv; i++)
			//	if (neib[i]==NULL)
			//	{	glBegin(GL_LINE_LOOP);
			//		for (int j=1; j<=Nfv1; j++)
			//		{	x=vert[(i+j)%Nv]->x;
			//			glVertex3f(x[0],x[1],x[2]);
			//			x=vert[(i+j+1)%Nv]->x;
			//			glVertex3f(x[0],x[1],x[2]);
			//		}
			//		glEnd();
			//	}
			//}
			//else
			//{
			glBegin(GL_LINES);
			for (int i=0; i<Nv; i++)
			{	x=vert[i]->x;
				glVertex3f(x[0],x[1],x[2]);
				x=vert[(i+1)%Nv]->x;
				glVertex3f(x[0],x[1],x[2]);
			}
			x=vert[0]->x;
			glVertex3f(x[0],x[1],x[2]);
			x=vert[2]->x;
			glVertex3f(x[0],x[1],x[2]);
			x=vert[1]->x;
			glVertex3f(x[0],x[1],x[2]);
			x=vert[3]->x;
			glVertex3f(x[0],x[1],x[2]);
			glEnd();
			//}
			if (c==first) break;
		}
		glEnable(GL_LIGHTING);
	}
	if(showpar[showNodes]) 
	{
		float	ncolor[3],bcolor[3];
		int
			showboundarynodes=showpar[showBoundaryVertexes];
		DNode	*dnode_root=dom->dnode_root;
		dom->getElementColor(nodes,ncolor);
		dom->getElementColor(faces,bcolor);
		if(showpar[showSpheres])
		{	GLfloat	radius=G.node[size];
			GLfloat mat[4],//- = { 0.6, 0.8, 0.3, 1.0 },
				shine=0.6;
			DNode	*node=dnode_root;
			mat[3]=1.0;
			glEnable(GL_LIGHTING);
			glMaterialf (GL_FRONT, GL_SHININESS, shine*128.0);
			if (dnode_root!=NULL)
			do 
			{	double	*x=node->x;
				glTranslatef ((GLfloat)x[0], (GLfloat)x[1], (GLfloat)x[2]); 
				if (node->state.boundary)
				{
					mat[0]=bcolor[red];mat[1]=bcolor[green];mat[2]=bcolor[blue];
				}
				else
				{
					mat[0]=ncolor[red];mat[1]=ncolor[green];mat[2]=ncolor[blue];
				}
				glMaterialfv (GL_FRONT, GL_AMBIENT, mat);
				glMaterialfv (GL_FRONT, GL_DIFFUSE, mat);
				glMaterialfv (GL_FRONT, GL_SPECULAR, mat);
				glutSolidSphere (radius, 12, 12);
				glTranslatef (-(GLfloat)x[0],-(GLfloat)x[1],-(GLfloat)x[2]); 
				node=node->next;
			}	while(node!=dnode_root);
		}
		else
		{	double
				vmn=(double)maxElementStatus,
				vmx=-1.0,
				vav=.5*(vmn+vmx),
				alow =PI/(2.*(vav-vmn)),
				ahigh=PI/(2.*(vmx-vav));
			DNode	*node=dnode_root;
			glPointSize(G.node[size]);
			glDisable(GL_LIGHTING);
			glBegin(GL_POINTS);
			if (dnode_root!=NULL)
			do 
			{	double	*x=node->x;
				if (node->state.boundary)
				{///	glColor3f(bcolor[red], bcolor[green], bcolor[blue]);
					double	r,g,b;
					double	a,
						var=(double)node->type;
					if (var<vmn)
					{	r=g=0.0; b=1.0; }
					else
					if (var>vmx)
					{	r=1.0; g=b=0.0; }
					else
					if (var<vav)
					{	a=(double)(var-vmn)*alow;
						b=cos(a); g=sin(a);
						r=0.0;
					}
					else
					{	a=(double)(var-vav)*ahigh;
						g=cos(a); r=sin(a);
						b=0.0;
					}
					glColor3f(r,g,b);
					glVertex3f((GLfloat)x[0],(GLfloat)x[1],(GLfloat)x[2]);
				}
				else
				{	if (!showboundarynodes)
					{	glColor3f(ncolor[red], ncolor[green], ncolor[blue]);
						glVertex3f((GLfloat)x[0],(GLfloat)x[1],(GLfloat)x[2]);
					}
				}
				node=node->next;
			}	while(node!=dnode_root);
			glEnd();
			glEnable(GL_LIGHTING);
		}
	}
	if(showpar[showBonds])
	{
		float	edgescolor[3];
		DNode	*p,*dnode_root=dom->dnode_root;
		dom->getElementColor(edges,edgescolor);
		glDisable(GL_LIGHTING);
		glLineWidth (G.line[thickness]);
		glBegin(GL_LINES);
		p=dnode_root;
		if(dnode_root!=NULL)
		do 
		{	double	*x=p->x;
			Bond	*bond=p->bond;
			int	ibond=0;
			if(bond!=NULL)
			do 
			{	double	*y,z[DIM];
				DNode	*q=bond->node;
				if (q==NULL)
				{
					glColor3f(edgescolor[red], edgescolor[green], edgescolor[blue]);
					y=bond->x;
					for(int i=0;i<DIM;i++)z[i]=x[i]+y[i];
					y=z;
				}
				else
				{
					glColor3f(vecval[red], vecval[green], vecval[blue]);
					y=q->x;
				}
				glVertex3f
				(	(GLfloat)x[0],
					(GLfloat)x[1],
					(GLfloat)x[2]
				);
				glVertex3f
				(	(GLfloat)y[0],
					(GLfloat)y[1],
					(GLfloat)y[2]
				);
				bond=bond->next;ibond++;
			}	while(bond!=p->bond);
			p=p->next;
		}	while(p!=dnode_root);
		glEnd();
	}
	if(showpar[showFrame]&&dom->tool!=NULL) 
	{
		float	color[3];
		Frame	
			*first_framenode=dom->tool->first_framenode,
			*last_framenode=dom->tool->last_framenode;
		dom->getElementColor(cells,color);
		glColor3f(color[red], color[green], color[blue]);
		glPointSize(G.node[size]);
		glDisable(GL_LIGHTING);
		glBegin(GL_POINTS);
		if (first_framenode!=NULL)
		for 
		(	Frame	*node=first_framenode;
			; node=node->next
		)
		{	double	*x=node->x;
			glVertex3f((GLfloat)x[0],(GLfloat)x[1],(GLfloat)x[2]);
			if (node==last_framenode) break;
		}
		glEnd();
		glEnable(GL_LIGHTING);
	}
	if (showpar[showTool]&&tool!=NULL)
	{
		double	r,x[DIM];
		GLfloat mat[4] = { 0.8, 0.7, 0.5, 1.0 };
		mat[0]=(GLfloat)toolval[red];
		mat[1]=(GLfloat)toolval[green];
		mat[2]=(GLfloat)toolval[blue];
		glEnable(GL_LIGHTING);
		glMaterialfv (GL_FRONT, GL_AMBIENT, mat);
		glMaterialfv (GL_FRONT, GL_DIFFUSE, mat);
		glMaterialfv (GL_FRONT, GL_SPECULAR, mat);
		dom->tool->getpos(x);
		r=dom->tool->getradius();
		glTranslatef (x[0], x[1], x[2]); 
		glutSolidSphere (r, 15, 15);
		glTranslatef (-x[0],-x[1],-x[2]); 
	}
	if (showpar[showVariables])
		displayVariables(dom);
	if(showpar[showBoundaryFaces])
	{
		if(bface_root!=NULL)
		{	GLfloat mat[4] = { 0.8, 0.7, 0.5, 1.0 },
				shine=0.6;
			glEnable(GL_LIGHTING);
			glMaterialfv (GL_FRONT, GL_AMBIENT, mat);
			glMaterialfv (GL_FRONT, GL_DIFFUSE, mat);
			glMaterialfv (GL_FRONT, GL_SPECULAR, mat);
			glMaterialf (GL_FRONT, GL_SHININESS, shine*128.0);
			///BFaceList	*bface_root=dom->bface_root;	
			//glShadeModel (GL_FLAT);
			glShadeModel (GL_SMOOTH);
			glBegin(GL_TRIANGLES);
			for
			(	BFaceList	*c=bface_root->next;
				; c=c->next
			)
			{	DCell	*cell=c->cell;
				DNode	**vert=cell->vert;
				int	iv=c->iface;
				GLfloat	norm[DIM];
				if (cell->facetype[iv]>=boundary)
				{	for (int i=0; i<DIM; i++)norm[i]=(GLfloat)c->norm[i];
					glNormal3fv(norm);
					for (int j=0; j<Nfv; j++)
					{	double *x=vert[(iv+j+1)%Nv]->x;
						glVertex3f((GLfloat)x[0],(GLfloat)x[1],(GLfloat)x[2]);
					}
				}
				if (c==bface_root) break;
			}
			glEnd();
			if (showpar[showBoundaryVectors])
			{	//Show normal vectors
				glDisable(GL_LIGHTING);
				glColor3f(vecval[red], vecval[green], vecval[blue]);
				glBegin(GL_LINES);
				for
				(	BFaceList	*c=bface_root->next;
					; c=c->next
				)
				{	DCell	*cell=c->cell;
					int	iv=c->iface;
					double
						veclen=sqrt(c->area),
						*x=c->x;
					GLfloat	y[DIM];
						for (int i=0; i<DIM; i++)y[i]=x[i]+vecval[length]*veclen*(GLfloat)c->norm[i];
					glVertex3f
					(	(GLfloat)x[0],
						(GLfloat)x[1],
						(GLfloat)x[2]
					);
					glVertex3f
					(	(GLfloat)(y[0]),
						(GLfloat)(y[1]),
						(GLfloat)(y[2])
					);
					if (c==bface_root) break;
				}
				glEnd();
			}
		}
	}
	if(showpar[showParticles]) 
		displayParticles(dom);
}
//?void	displayStaticDomain(int idomain)
//?{
//?	Domain	*dom=domain_root+idomain;
//?	struct Node	*node=dom->node;
//?	struct Face	*face=dom->face;
//?	struct Cell	*cell=dom->cell;
//?
//?//#ifndef GRIDGEN
//?//	if (showIsoSurfaces)
//?//	if (S==NULL)
//?//	{	printf("No isosurfaces to show\n");fflush(stdout);
//?//		showIsoSurfaces=0;
//?//	}
//?//	else
//?//	{	int i,j;
//?/////		glEnable(GL_LIGHTING);
//?/////		glShadeModel (GL_FLAT);
//?////		glDisable(GL_LIGHTING);
//?////		glLineWidth (G.line[thickness]);
//?////		glColor3f(G.line[red], G.line[green], G.line[blue]);
//?//		switch (G.surfdisptype)
//?//		{
//?//			case solidsurface:
//?//			for (i=0; i<S->ne; i++)
//?//			{	struct IsosurfElement	*e=S->e+i;
//?//				glBegin(GL_TRIANGLES);
//?//				for (j=0; j<Nfv; j++)
//?//				{	double *x=G.X+DIM*e->vert[j];
//?//					GLfloat
//?//						*n=e->norm+DIM*j;
//?//					glVertex3f((GLfloat)x[0],(GLfloat)x[1],(GLfloat)x[2]);
//?//					glNormal3fv(n);
//?//				}
//?//				glEnd();
//?//			}
//?//			break;
//?//			case gridlines:
//?//			default:
//?//			glDisable(GL_LIGHTING);
//?//			glLineWidth (G.line[thickness]);
//?//			glColor3f(G.isoline[red], G.isoline[green], G.isoline[blue]);
//?//			for (i=0; i<S->ne; i++)
//?//			{	struct IsosurfElement	*e=S->e+i;
//?//				glBegin (GL_LINE_LOOP);
//?//				for (j=0; j<Nfv; j++)
//?//				{	double *x=G.X+DIM*e->vert[j];
//?//					glVertex3f((GLfloat)x[0],(GLfloat)x[1],(GLfloat)x[2]);
//?//				}
//?//				glEnd();
//?//			}
//?//			glEnable(GL_LIGHTING);
//?//		}
//?//	  //	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//?//		//	glBegin(GL_TRIANGLE_STRIP);
//?//		//	for (i=0; i<S->ne; i++)
//?//		//	{	struct IsosurfElement	*e=S->e+i;
//?//		//		for (j=0; j<Nfv; j++)
//?//		//		{	int	k;
//?//		//			double
//?//		//				*x=G.X+DIM*e->vert[j];
//?//		//			GLfloat v[3],
//?//		//				*n=e->norm+DIM*j;
//?//		//			for (k=0; k<DIM; k++)
//?//		//				v[k]=(GLfloat)x[k];
//?//		//			glNormal3fv(n);
//?//		//			glVertex3fv(v);
//?//		//		}
//?//		//	}
//?//		//	glEnd();
//?////		glEnable(GL_LIGHTING);
//?//		glShadeModel (GL_SMOOTH);
//?//	}
//?//#endif
//?	if(showpar[showVariables])
//?		displayVariables(dom);
//?	if(showpar[showNodes]) 
//?	{	float	color[3];
//?		dom->getElementColor(nodes,color);
//?		showGridElements(dom->Ne[nodes],dom->coordinates[nodes].val,color);
//?	}
//?	if(showpar[showFaceCenters]) 
//?	{	float	color[3];
//?		dom->getElementColor(faces,color);
//?		showGridElements(dom->Ne[faces],dom->coordinates[faces].val,color);
//?	}
//?	if(showpar[showCellCenters]) 
//?	{	float	color[3];
//?		dom->getElementColor(cells,color);
//?		showGridElements(dom->Ne[cells],dom->coordinates[cells].val,color);
//?	}
//?	if(showpar[showBoundaryVertexes]) 
//?	{	int	nbv=dom->getNoBoundaryVertexes();
//?		double
//?			vmn=(float)dead,
//?			vmx=-1.0,
//?			vav=.5*(vmn+vmx),
//?			alow =PI/(2.*(vav-vmn)),
//?			ahigh=PI/(2.*(vmx-vav)),
//?			*X=dom->coordinates[nodes].val,
//?			*x;
//?			float	a;
//?		struct BoundaryVertex *i;
//?		glPointSize(G.node[size]);
//?		glDisable(GL_LIGHTING);
//?		glBegin(GL_POINTS);
//?		for (i=dom->Bv; i-dom->Bv<nbv; i++)
//?		{	float	r,g,b;
//?			double	var=(double)i->b;
//?			x=X+DIM*i->i;
//?			if (var<vmn)
//?			{	r=g=0.0; b=1.0; }
//?			else
//?			if (var>vmx)
//?			{	r=1.0; g=b=0.0; }
//?			else
//?			if (var<vav)
//?			{	a=(float)(var-vmn)*alow;
//?				b=cos(a); g=sin(a);
//?				r=0.0;
//?			}
//?			else
//?			{	a=(float)(var-vav)*ahigh;
//?				g=cos(a); r=sin(a);
//?				b=0.0;
//?			}
//?			glColor3f(r,g,b);
//?			glVertex3f((GLfloat)x[0],(GLfloat)x[1],(GLfloat)x[2]);
//?		}
//?		glEnd();
//?		glEnable(GL_LIGHTING);
//?	}
//?	if(showpar[showBoundaryFaceCenters])
//?	{	int	nbf=dom->getNoBoundaryFaces();
//?		double
//?			vmn=(float)boundary,
//?			vmx=-1.0,
//?			vav=.5*(vmn+vmx),
//?			alow =PI/(2.*(vav-vmn)),
//?			ahigh=PI/(2.*(vmx-vav)),
//?			*Y=dom->coordinates[faces].val;
//?			float	a;
//?		glPointSize(G.node[size]);
//?		glDisable(GL_LIGHTING);
//?		glBegin(GL_POINTS);
//?		for (int	i=0; i<nbf; i++)
//?		{	float	r,g,b;
//?			struct Face	*f=dom->face+i;
//?			double
//?				var=(double)f->cell[1],
//?				*y=Y+DIM*i;
//?			if (var<vmn)
//?			{	r=g=0.0; b=1.0; }
//?			else
//?			if (var>vmx)
//?			{	r=1.0; g=b=0.0; }
//?			else
//?			if (var<vav)
//?			{	a=(float)(var-vmn)*alow;
//?				b=cos(a); g=sin(a);
//?				r=0.0;
//?			}
//?			else
//?			{	a=(float)(var-vav)*ahigh;
//?				g=cos(a); r=sin(a);
//?				b=0.0;
//?			}
//?			glColor3f(r,g,b);
//?			glVertex3f((GLfloat)y[0],(GLfloat)y[1],(GLfloat)y[2]);
//?		}
//?		glEnd();
//?		glEnable(GL_LIGHTING);
//?	}
//?	if(showpar[showBoundaryFaces])
//?	{
//?		int
//?			lighting=dom->getLighting(faces),
//?			nbf=dom->getNoBoundaryFaces();
//?		float	color[3];
//?		double
//?			*X=dom->coordinates[nodes].val,	
//?			*x;
//?		if (lighting)
//?		{
//?			GLfloat mat[4] = { 0.8, 0.7, 0.5, 1.0 },
//?				shine=0.6;
//?			glMaterialfv (GL_FRONT, GL_AMBIENT, mat);
//?			glMaterialfv (GL_FRONT, GL_DIFFUSE, mat);
//?			glMaterialfv (GL_FRONT, GL_SPECULAR, mat);
//?			glMaterialf (GL_FRONT, GL_SHININESS, shine*128.0);
//?			glEnable(GL_LIGHTING);
//?			glShadeModel (GL_SMOOTH);
//?		}
//?		else
//?		{	dom->getElementColor(faces,color);
//?			glDisable(GL_LIGHTING);
//?			glLineWidth (G.line[thickness]);
//?			glColor3f(color[red], color[green], color[blue]);
//?		}
//?		glBegin(GL_TRIANGLES);
//?		for (int i=0; i<nbf; i++)
//?		{	struct Face	*f=face+i;
//?			double	*y=dom->coordinates[faces].val+DIM*i;
//?			GLfloat norm[DIM];
//?			for (int j=0; j<DIM; j++)
//?				norm[j]=(GLfloat)f->norm[j];
//?			glNormal3fv(norm);
//?			for (int j=0; j<Nfv; j++)
//?			{
//?				x=X+DIM*f->vert[j];
//?				glVertex3f((GLfloat)x[0],(GLfloat)x[1],(GLfloat)x[2]);
//?			}
//?		}
//?		glEnd();
//?		if (!lighting)
//?		glEnable(GL_LIGHTING);
//?	}
//?	if(showpar[showBoundaryGrid]) 
//?	{	int
//?			nbf=dom->getNoBoundaryFaces();
//?		float	color[3];
//?		double
//?			*X=dom->coordinates[nodes].val,
//?			*x;
//?		glDisable(GL_LIGHTING);
//?		glLineWidth (G.line[thickness]);
//?		dom->getElementColor(edges,color);
//?		glColor3f(color[red], color[green], color[blue]);
//?		for (struct Face *f=face; f-face<nbf; f++)
//?		{	glBegin (GL_LINE_LOOP);
//?			for (int j=0; j<Nfv; j++)
//?			{
//?				x=X+DIM*f->vert[j];
//?				glVertex3f((GLfloat)x[0],(GLfloat)x[1],(GLfloat)x[2]);
//?			}
//?			glEnd();
//?		}
//?		glEnable(GL_LIGHTING);
//?	}
//?#ifdef OLDVERSION
//?	if(showpar[showVectors]) 
//?	{	int	nvar=dom->getNoVariables();
//?		double	*x,*v,*u;
//?		if (G.ivec>=0&&G.ivec<nvar)
//?		{
//?			glLineWidth (vecval[thickness]);
//?			glDisable(GL_LIGHTING);
//?			glColor3f(vecval[red], vecval[green], vecval[blue]);
//?			switch (G.var[G.ivec].type)
//?			{	case node:
//?				u=G.V+G.var[G.ivec].i*G.n*DIM;
//?				for (x=G.X,v=u; v-u<G.n*DIM; x+=DIM,v+=DIM)
//?				{	glBegin(GL_LINES);
//?					glVertex3f
//?					(	(GLfloat)x[0],
//?						(GLfloat)x[1],
//?						(GLfloat)x[2]
//?					);
//?					glVertex3f
//?					(	(GLfloat)(x[0]+vecval[length]*v[0]),
//?						(GLfloat)(x[1]+vecval[length]*v[1]),
//?						(GLfloat)(x[2]+vecval[length]*v[2])
//?					);
//?					glEnd();
//?				}
//?				break;
//?				case face:
//?				break;
//?				case center:
//?				u=G.U+G.var[G.ivec].i*G.nc*DIM;
//?				for (x=G.Y,v=u; v-u<G.nc*DIM; x+=DIM,v+=DIM)
//?				{
//?					glBegin(GL_LINES);
//?					glVertex3f
//?					(	(GLfloat)x[0],
//?						(GLfloat)x[1],
//?						(GLfloat)x[2]
//?					);
//?					glVertex3f
//?					(	(GLfloat)(x[0]+vecval[length]*v[0]),
//?						(GLfloat)(x[1]+vecval[length]*v[1]),
//?						(GLfloat)(x[2]+vecval[length]*v[2])
//?					);
//?					glEnd();
//?				}
//?				break;
//?			}
//?			glEnable(GL_LIGHTING);
//?		}
//?		else
//?		{
//?			printf("Select a vector component to display (use \"Show vectors\" to see the options)\n");
//?		}
//?	}
//?#endif //OLDVERSION
//?	if(showpar[showParticles]) 
//?		displayParticles(dom);
//?}
void	displayDomain(int idomain)
{
	Domain	*dom=domain_root+idomain;
	switch(dom->type)
	{
		case dynamic:
			displayDynamicDomain(idomain);
			break;
		default:
		{	char msg[MAXLINLEN];
			sprintf
			(	msg,
				"Can not display domain %s: Unrecognized domain type: %d\n",
				dom->name,dom->type
			);
			displaymessage(msg);
		}
	}
}
/* ARGSUSED1 */
void	getScaling
(
	double	&l0,
	double	&l1
)
{
	l0=lmin,l1=lmax;
}
void	getXLimits
(	double	*x0,
	double	*x1
)
{
	for (int i=0; i<DIM; i++)
	{	x0[i]=xmin[i];
		x1[i]=xmax[i];
	}
}
#ifndef GRIDGEN
void	printPartVelLimits()
{	extern void	getPartVelLimits
	(
		double	*xmin,
		double	*xmax
	);
	double
		x0[DIM],x1[DIM];
	printf
	(	"Particle velocity limits:\n\tX: %g, %g\n\tY: %g, %g\n\tZ: %g, %g\n",
		x0[0],x1[0],x0[1],x1[1],x0[2],x1[2]
	);
}
#endif
void	printVariables
(
	int	n
)
{
	for (int idom=0;idom<ndomains; idom++)
	{	Domain	*dom=domain_root+idom;
		int	nvar=dom->getNoVariables();
		printf("Domain %d: '%s'\n",idom,dom->name);
		for (int ivar=0; ivar<nvar; ivar++)
		{	struct Variable *var=dom->variable+ivar;
			int	
				type=var->type,
				dim=var->dimension;
			double	vmin,vmax;
			printf
			(	"\t%d: %s, %s, type='%s', dim=%d:",
				ivar,var->name[0],var->name[1],elementype[type],dim
			);fflush(stdout);
			var->computeLimits();
			for (int icomp=0; icomp<dim; icomp++)
			{	var->getLimits(icomp,vmin,vmax);
				printf(" (%g,%g)",vmin,vmax);
			}
			printf("\n");fflush(stdout);
		}
	}
}
void	Exit()
{
	char	s[MAXLINLEN];
	if (animation) 
	{
		animation=0;
	}
	printf("Exit (y/n)? ");fflush(stdout);
	fgets(s,MAXLINLEN,stdin);
	if (*s!='y') return;
	if (!finished)finish();
	cleanup();
	exit(0);
}
void	Quit()
{
	char	s[MAXLINLEN];
	if (animation) 
	{
		animation=0;
	}
	printf("Quit (y/n)? ");fflush(stdout);
	fgets(s,MAXLINLEN,stdin);
	if (*s!='y') return;
	::cleanup();
	exit(0);
}
void	printDomains()
{
	printf("Domain list\n");
	for (int i=0; i<ndomains; i++)
	{	int	nvar=domain_root[i].getNoVariables();
		printf
		(	"%d. %s: variables: %d\n",
			i+1,domain_root[i].name,nvar
		);
	}
	if (idomain>=0)
		printf("Current domain: %d. %s\n",idomain+1,domain_root[idomain].name);
}
void	selectCurrentDomain()
{	int	i=-1;
	printDomains();
	if(ndomains==1)
	{	idomain=0; return;}
	while(i<0||i>ndomains)
	{	printf("Set current domain (1..%d domain number): ",ndomains);
		scanf("%d",&i);
/*
#ifdef WITH_MPI
		COMM_WORLD.Bcast(&i, 1, MPI_INT, 0);
#endif
*/
		if (i<1||i>ndomains) printf("Invalid domain number: %d\n",i);
	}
	idomain=i-1;
		printf("Current domain set to ");
//	if (idomain>0)
		printf("%d: %s\n",i,domain_root[idomain].name);
//	else
//		printf("none\n");
}
void	selectVariable()
{	int	i,ivar,nvar;

	if (idomain<0||idomain>ndomains)
	{	printf("Select current domain before selecting a variable\n");
		return;
	}
	nvar=domain_root[idomain].getNoVariables();
	domain_root[idomain].printVariables();
	ivar=domain_root[idomain].getDisplayVariable();
	if (ivar>=0)	
		printf
		(	"Current display variable: %d: %s\n",
			ivar+1,domain_root[idomain].variable[ivar].name
		);
	else
		printf("No variable is selected to display\n");
	while (1)
	{	char	buf[MAXLINLEN];
		printf("Select a new variable (0-none, 1..%d - new one): ",nvar);
		fgets(buf,MAXLINLEN,stdin);
/*
#ifdef WITH_MPI
		COMM_WORLD.Bcast(buf, MAXLINLEN, MPI_CHAR, 0);
#endif
*/
		sscanf(buf,"%d",&i);
		if (i<0||i>nvar)printf("Wrong variable number: %d\n",i);
		else break;
	}
	ivar=i-1;
	domain_root[idomain].setDisplayVariable(ivar);
	vars[idomain]=ivar;
	if (ivar>=0)
		printf
		(	"Variable %s of domain %s is selected for display\n",
			domain_root[idomain].variable[ivar].name,domain_root[idomain].name
		);
	else
		printf
		(	"No variables will be shown for domain %s\n",
			domain_root[idomain].name
		);
	if (ivar<0||domain_root[idomain].variable[ivar].dimension==0)
	{
		coms[idomain]=0;
		domain_root[idomain].setDisplayVariableComp(0);
	}
	else
	{
		while (1)
		{	int	dim=domain_root[idomain].variable[ivar].dimension;
			char	buf[MAXLINLEN];
			printf
			(	"Select variable component to display (0-vector plot, 1..%d-color plot): ",
				dim
			);
		fgets(buf,MAXLINLEN,stdin);
/*
#ifdef WITH_MPI
		COMM_WORLD.Bcast(buf, MAXLINLEN, MPI_CHAR, 0);
#endif
*/
			sscanf(buf,"%d",&i);
			if (i<0||i>DIM) printf("\nWrong variable component: %d\n",i);
			else 
			if (i==0&&dim!=DIM) 
				printf("\nCan't use vector plot for variables of dimension %d\n",dim);
			else break;
		}
		coms[idomain]=i-1;
		domain_root[idomain].setDisplayVariableComp(i-1);
	}
	if(showpar[showVariables]==0)showpar[showVariables]=1;
}
void	menu(int value)
{
	printf("Menu %d: ",value);
	switch (value) 
	{
		case 0:
			showpar[showBoundaryVertexes] = showpar[showBoundaryVertexes]==0?1:0;
			printf("Show Boundary Vertexes = %d\n",showpar[showBoundaryVertexes]);
			break;
		case 1:
			showpar[showBoundaryGrid] = showpar[showBoundaryGrid]==0?1:0;
			printf("Show Boundary Grid = %d\n",showpar[showBoundaryGrid]);
			break;
		case 2:
			showpar[showBoundaryFaceCenters] = showpar[showBoundaryFaceCenters]==0?1:0;
			printf("Show Boundary Faces = %d\n",showpar[showBoundaryFaceCenters]);
			break;
///#ifndef GRIDGEN
///		case 5:
//			showIsoSurfaces = showIsoSurfaces==0?1:0;
//			printf("Show IsoSurfaces = %d\n",showIsoSurfaces);
//			printf("Not impemented\n");fflush(stdout);
//			if (showIsoSurfaces==0)
//			{	//Delete isosurfaces
//				deleteIsoSurfaces(&S);
//			}
//			else
//			{	//Create isosurfaces
//				createIsoSurfaces(&S);
//				if (S==NULL) showIsoSurfaces = 0;
//			}
//			break;
//#else
		case 3:
			showpar[showNodes] = showpar[showNodes]==0?1:0;
			printf("Show Grid Nodes = %d\n",showpar[showNodes]);
			break;
		case 4:
			showpar[showFaceCenters] = showpar[showFaceCenters]==0?1:0;
			printf("Show Face Centers = %d\n",showpar[showFaceCenters]);
			break;
		case 5:
			showpar[showCellCenters] = showpar[showCellCenters]==0?1:0;
			printf("Show Cell Centers = %d\n",showpar[showCellCenters]);
			break;
		case 6:
			showpar[showGrid] = showpar[showGrid]==0?1:0;
			printf("Show Grid = %d\n",showpar[showGrid]);
			break;
///#endif
		case 7:
			showpar[showParticles] = showpar[showParticles]==0?1:0;
			printf("Show particles = %d\n",showpar[showParticles]);
			if (showpar[showParticles]==1)
			{
				//printf("No particles yet\n");P.n=0;
				domain_root[idomain].getPoints(P.m,P.n,P.origin,P.first,P.last,P.X);
				printf("Number of particles = %d\n",P.n);
			}
			break;
		case 8:
			showpar[showBonds] = showpar[showBonds]==0?1:0;
			printf("Show Bonds = %d\n",showpar[showBonds]);
			break;
		case 9:
			if (idomain<0||idomain>ndomains)
			{	printf("No domain is selected\n");
				break;
			}
			{	int
					  ivar=domain_root[idomain].getDisplayVariable(),
					  icomp=domain_root[idomain].getDisplayVariableComp();
				double min,max;
				if (ivar<0&&showpar[showVariables]==0)
				{	printf("No variable is selected to display\n");
					break;
				}
				if (idomain>=0&&ivar>=0&&icomp>=0)
				{
					domain_root[idomain].variable[ivar].getLimits
					(	icomp,
						min,max
					);
					if (colorscale.relative)
					{	G.node[vmin]=(float)min;
						G.node[vmax]=(float)max;
					}
				}
				if (option.verbose)
				{	int ivar;
					if (idomain>=0)
					{	ivar=domain_root[idomain].getDisplayVariable();
						printf("Current Domain = %d: %s\n",idomain+1,domain_root[idomain].name);
						printf("Current Variable = ");
						if (ivar>=0)
							printf
							(	"%d: %s",
								ivar+1,domain_root[idomain].variable[ivar].name
							);
						else
							printf("undefined");
						printf("\n");
					}
					else
					{	printf("Warning: no current domain was selected\n");
					}
				}
			}	
			showpar[showVariables] = showpar[showVariables]==0?1:0;
			printf("Show Variables = %d\n",showpar[showVariables]);
			break;
		case 10:
			showpar[showAxes] = showpar[showAxes]==0?1:0;
			printf("Show Axes = %d\n",showpar[showAxes]);
			break;
		case 11:
			selectCurrentDomain();
			break;
		case 12:
			selectVariable();
			break;
		case 13:
			init();
			break;
		case 14:
			if (finished) break;
			animation=animation==0?1:0;
			if (animation) 
				glutIdleFunc(animate);
			else
			{	glutIdleFunc(NULL);
				printf("\n");fflush(stdout);
			}
		break;
		case 15:
			helpDisplay();
		break;
		case 16:
		printf(ABOUT);
		break;
		case 17:
		printf("\nSend your comments to andrei@smirnov.mae.wvu.edu\n");
		cleanup();
		exit(0);
		break;
	}
	fflush(stdout);
	//glutPostRedisplay();
}
void displayMenu()
{
	printf("----------------- MENU -----------------\n");
	printf(" ACTION                   (key-shortcut)\n");
	printf(" 0. Quit menu             (<Enter>)\n");
	printf("Toggle the display of        \n");
	printf(" 1. Boundary vertexes        \n");
	printf(" 2. Boundary grid         (g)\n");
	printf(" 3. Boundary faces           \n");
	printf(" 4. Nodes                 (n)\n");
	printf(" 5. Faces                 (f)\n");
	printf(" 6. Cells                 (c)\n");
	printf(" 7. Edges                 (e)\n");
	printf(" 8. Particles             (p)\n");
	printf(" 9. Particle bonds        (b)\n");
	printf("10. Variables             (v)\n");
	printf("11. Axes                  (a)\n");
	printf("Actions:\n");
	printf("12. Select current domain (d)\n");
	printf("13. Select variable       (v)\n");
	printf("14. Initialize            (r)\n");
	printf("15. Run/Stop              (r)\n");
	printf("16. Help                  (?)\n");
	printf("17. About                    \n");
	printf("18. Quit program          (q)\n");
}
void consoleMenu()
{	int	selection=-1;
	do
	{	char	buf[MAXLINLEN];
		displayMenu();
		printf("Select: ");
		selection=0;
		fgets(buf,MAXLINLEN,stdin);
		sscanf(buf,"%d",&selection);
		if (selection==0) return; 
		if (selection>0&&selection<=15)
			menu(selection-1);
		else printf("Wrong selection: %d\n",selection);
		//else	break;
	}	while (selection!=0);
}
void	printParticleVariables()
{/*	extern void	getPartVarLimits
	(
		struct Variables **v
	);
	double
		x0[DIM],x1[DIM];
 */
}
void	setBackgroundRun()
{
	showpar[showRun]=0;
}
void	setForegroundRun()
{
	showpar[showRun]=1;
}
void	toggleBoundaryVectors()
{
	showpar[showBoundaryVectors]=showpar[showBoundaryVectors]?0:1;
}
void	toggleSpheres()
{
	showpar[showSpheres]=showpar[showSpheres]?0:1;
	printf("Show Spheres = %d\n",showpar[showSpheres]);
}
void	runmany()
{	int n=0;
	printf("Enter number or iterations: ");FLUSH;
	scanf("%d",&n);
	printf("Executing %d iterations ...",n);FLUSH;
	run(n);
	printf(" done\n");
}
void	setWriteFormat()
{	using namespace Output;
	int	selection=-1;
	do
	{	char	buf[MAXLINLEN];
		printf("Output formats:\n");
		for(int i=0;i<maxoutypes;i++)
			printf("%d - %s\n",i,Output::outypename[i]);
		printf("Select: ");
		selection=0;
		fgets(buf,MAXLINLEN,stdin);
		sscanf(buf,"%d",&selection);
		if (selection>=0&&selection<maxoutypes)break;
		printf("Wrong selection: %d\n",selection);
	}	while (selection!=0);
	outype=(OutputTypes)selection;
	printf("Selected output type: %s\n",Output::outypename[outype]);FLUSH;
}
void	writeData()
{
	char	*s,buf[MAXLINLEN+1],
		outfilename[MAXLINLEN];
	Domain	*dom=domain_root+idomain;
	using namespace	Output;
	sprintf(outfilename,"%s",dom->name);
#ifdef GRIDGEN
	printf("Write grid (g), boundary (b)? ");fflush(stdout);
#else
	printf("Write:\n\t grid (g)\n\t boundary (b)\n\t current-variable (v)\n\t field-data (f)\n\t particle-data (p)?\nEnter your choice: ");fflush(stdout);
#endif
	fgets(buf,MAXLINLEN,stdin);
	for (s=buf; isspace(*s); s++);
	switch (*s)
	{
		case 'g':
		dom->saveGeom(outfilename);
		break;
		case 'b':
		outputBoundary();
		break;
#ifndef GRIDGEN
		case 'v':
		dom->saveDispVarDefault(dom->getDisplayVariable(),outfilename);
		break;
		case 'f':
		dom->outputData(outfilename);
		break;
		case 'p':
		printf("Ouput particles not implemented yet\n");
//		outputParticles();
		break;
#endif
		default:
		printf("Operation skipped");
	}
	printf("\n");
}
void	toggleWindowDump()
{
	printf("Enter window dump time interval: ");FLUSH;
	scanf("%g",&wdtime);
	if(wdtime>0.0)
	{	showpar[dumpWindow]=1;
		iwdump=0;nwdump=(int)(wdtime/runtime.step+0.5);
		printf("Window dump = %d, dump time interval = %g\n",showpar[dumpWindow],wdtime);
	}
	else
	{	showpar[dumpWindow]=0;
		printf("Window dump disabled\n",showpar[dumpWindow],wdtime);
	}
}
void	dumpwindow()
{	using namespace Output;
	xwd(windowname);
}
void	commandMode()
{
	static struct Command
	{
		char	*name;
		void	(*f)();
	}	command[]	= 
	{
///#ifndef GRIDGEN
///		{"pv",showParticleVariables},
///		{"pvel",showPartVelLimits},
		{"bv",toggleBoundaryVectors},
		{"dom",selectCurrentDomain},
		{"v",selectVariable},
///#endif
		{"m",consoleMenu},
		{"w",writeData},
		{"wf",setWriteFormat},
		{"lim",printGridXLimits},
		{"gvec",printGridVecLimits},
		{"bg",setBackgroundRun},
		{"fg",setForegroundRun},
		{"?",helpCommand},
		{"h",helpCommand},
		{"help",helpCommand},
		{"r",runmany},
		{"sp",toggleSpheres},
		{"wd",toggleWindowDump},
		{"e",Exit},
		{"q",Quit},
		{".",NULL},
		{"\0",NULL}
	},	*cmd;
	char	buf[MAXLINLEN+1],*p,*s;
	printf("Command mode\n");
	while (1)
	{
		printf(":");fflush(stdout);
#ifdef WITH_DOVE
	if(iproc==0)
#endif
	{
		fgets(buf,MAXLINLEN,stdin);
	}
/*
#ifdef WITH_MPI
		COMM_WORLD.Bcast(buf, MAXLINLEN, MPI_CHAR, 0);
#endif
*/
		for (s=buf; isspace(*s); s++);
		if ((p=strchr(s,'\n'))!=NULL) *p='\0';
		for (cmd=command; *cmd->name!='\0'; cmd++)
			if (strcmp(s,cmd->name)==0) break;

	
		switch (*cmd->name)
		{
			case '\0':
			if (*s=='\0')
			{	printf("Display-mode\n");fflush(stdout);
				return;
			}
			else
				printf("Unknown command\n");
			break;
			case '.':
				printf("\nDisplay-mode\n");fflush(stdout);
				return;
			default:
			(*cmd->f)();
		}
	}
}

/* Main */


/******** GLOBAL GUI-FUNCTIONS **********/

void	display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glPushMatrix();
	displayAxes();
	for (int	idomain=0; idomain<ndomains; idomain++)
		displayDomain(idomain);
	/*
	{	Domain	*dom=domain_root+idomain;
		if(iproc==dom->iproc)	displayDomain(idomain);
	}
	*/
	glPopMatrix();
//	glFlush();
//glutSwapBuffers();
	glXSwapBuffers(dpy, win);
}
void	animate()
{
	if (run(1)==0)
	{	printf("\n");
		finish();
		finished=1;
		animation=0;
//-		glutIdleFunc(NULL);
		display();
		return;
	}
	printf("\rTime =%11.4f ",runtime.current);fflush(stdout);

	// Custom output:
	//double x = fmod (runtime.current, 1.0);
	//if(x<runtime.step)
	//{	char outfilename[30];
	//	for(int i =0; i<ndomains; i++)
	//	{
	//		Domain *dom = domain_root + i;
	//		int ivar = dom->getDisplayVariable();
	//		sprintf(outfilename, "%s_%s", dom->name, dom->variable[ivar].name);
	//		dom->saveDispVarDefault(ivar,outfilename);
	//	}
	//}

	if (showpar[showRun])
	{///	refresh();
		display();
	}
	if (showpar[dumpWindow]) 
	{
		if (iwdump++%nwdump==0)
			dumpwindow();
	}
	if (showpar[showParticles]) printf("np=%d",P.n);fflush(stdout);
}
void	printGridXLimits()
{
	double
		xmin[DIM],xmax[DIM];
	for (int i=0; i<DIM; i++)
	{	xmin[i]=LARGE;
		xmax[i]=-LARGE;
	}
	for 
	(	Domain	*d=domain_root; d-domain_root<ndomains; d++
	)
	{	int	idomain=d-domain_root;
		double x0[DIM],x1[DIM];
		d->getGridLimits(x0,x1);
		for (int i=0; i<DIM; i++)
		{	if (xmin[i]>x0[i])xmin[i]=x0[i];
			if (xmax[i]<x1[i])xmax[i]=x1[i];
		}
		printf
		(	"Domain %d: Xmin=%g, Xmax=%g\nYmin=%g, Ymax=%g\nZmin=%g, Zmax=%g\n",
			idomain,x0[0],x1[0],x0[1],x1[1],x0[2],x1[2]
		);
	}
	printf
	(	"Global limits: xmin=(%g,%g,%g), xmax=(%g,%g,%g)\n",
		xmin[0],xmin[1],xmin[2],
		xmax[0],xmax[1],xmax[2]
	);
}
void	printGridVecLimits()
{
	double
		x0[DIM],x1[DIM];
//#ifndef GRIDGEN
//	extern void	getGridVecLimits
//	(
//		double	*xmin,
//		double	*xmax
//	);
//	getGridVecLimits(x0,x1);
//#endif
//	printf("Vector field limits:\n\tX: %g, %g\n\tY: %g, %g\n\tZ: %g, %g\n",x0[0],x1[0],x0[1],x1[1],x0[2],x1[2]);
}
void	helpCommand()
{
	printf("dom   - select a domain\n");
	printf("e     - save and exit\n");
	printf("h     - show help\n");
	printf("bg    - set background run (no automatic redisplay)\n");
	printf("fg    - set foreground run (automatic redisplay)\n");
	printf("lim   - display domain limits\n");
	printf("m     - select from a menu\n");
	printf("gvec  - display grid vector limits\n");
	printf("gv    - display grid variables limits\n");
	printf("pv    - display particle variables limits\n");
	printf("q     - quit without saving\n");
	//printf("r     - run\n");
	printf("var   - select a variable\n");
	printf("w     - write data at this time step\n");
	printf(".     - quit command mode\n");
}
void	reshape(int w, int h)
{
	float	
		r0=5.,
		r1=10.0;
	double
		lmin,lmax,
		xmin[DIM],xmax[DIM];
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	if (w <= h) 
//		glOrtho 
//		(	-r0, r0, 
//			-r0*(GLfloat)h/(GLfloat)w, 
//			 r0*(GLfloat)h/(GLfloat)w, 
//			-r1, r1
//		);
//		else 
//		glOrtho 
//		(	-r0*(GLfloat)w/(GLfloat)h, 
//			 r0*(GLfloat)w/(GLfloat)h, 
//			-r0, r0, -r1, r1
//		);
//	gluPerspective (40., (GLdouble)w/(GLdouble)h, -3., 8.);
	getXLimits(xmin,xmax);
	glOrtho 
	(	xmin[0], 
		xmax[0], 
		xmin[1], 
		xmax[1], 
		xmin[2], 
		xmax[2]
	);
	glLoadIdentity ();
	getScaling(lmin,lmax);
//	gluPerspective (45., (GLdouble)w/(GLdouble)h, .1*lmin, 10*lmax);
  	gluPerspective (50., (GLdouble)w/(GLdouble)h, 0.1*lmin, 10*lmax);
	glMatrixMode(GL_MODELVIEW);
}
void	mouse(int button, int state, int x, int y)
{
//	if (button == LEFT_BUTTON) 
//	{
//-		if (state == GLUT_DOWN) 
		if (
			button == GLFW_MOUSE_BUTTON_RIGHT ||
			button == GLFW_MOUSE_BUTTON_LEFT
		)
		{	lastx=x;
			lasty=y;
//			lastz=y;
//			down = 1;
//		} else 
//		{
//			down = false;
			movement=stay;
		}
//	}
	switch (button)
	{
		case GLFW_MOUSE_BUTTON_LEFT: //-LEFT_BUTTON:
		movement=rotate;
		break;
		case GLFW_MOUSE_BUTTON_MIDDLE: //-MIDDLE_BUTTON:
		movement=moveuv;
		lastx=x; lasty=y;
		break;
		case GLFW_MOUSE_BUTTON_RIGHT: //-RIGHT_BUTTON:
		movement=movew;
		lastz=x;
		break;
		default:
		movement=stay;
	}
}
void	motion(int x, int y)
{
	switch(movement)
	{
		case	rotate:
		glTranslatef(xo[0],xo[1],xo[2]);
		glRotatef(x-lastx, 0, 1, 1);
		glRotatef(y-lasty, 1, 0, 1);
		glTranslatef(-xo[0],-xo[1],-xo[2]);
		lastx=x; lasty=y;
		break;
		case	moveuv:
		glTranslatef(dx*(float)(x-lastx), dy*(float)(lasty-y), 0.0);
		lastx=x; lasty=y;
		break;
		case	movew:
		glTranslatef(0.0, 0.0, lastz-x>=0?-dz:dz);
		lastz=x;
		break;
	}
}
/* ARGSUSED3 */
void	keyboard(unsigned int key)
{
	using namespace Output;
	if (option.debug) printf("key=%c (%d)\n",(char)key,key);
	if (key<32||key>127)return;
/*
#ifdef WITH_MPI
	if(iproc==0)
	{
		for(int i=1; i<COMM_WORLD.Get_size(); i++)
		{
			COMM_WORLD.Send(&key, 1, MPI_UNSIGNED, i, 1);
		}
	}
#endif
*/
	switch (key) 
	{
		case '+':
			glTranslatef(0.0, 0.0, dz);
			break;
		case '-':
			glTranslatef(0.0, 0.0, -dz);
			break;
		case 'A':
		case 'a':
			showpar[showAxes] = showpar[showAxes]==0?1:0;
			printf("Show Axes = %d\n",showpar[showAxes]);
			break;
		case 'B':
			showpar[showBoundaryFaces] = showpar[showBoundaryFaces]==0?1:0;
			printf("Show Boundary Faces = %d\n",showpar[showBoundaryFaces]);
			if(domain_root[idomain].getNoVariables(boundary_faces)==0)
			{	domain_root[idomain].deleteRing(bface_root);
				domain_root[idomain].createBoundaryFaceList(bface_root);
			}
			break;
		case 'b':
			showpar[showToolGrid] = showpar[showToolGrid]==0?1:0;
			printf("Show Tool Grid = %d\n",showpar[showToolGrid]);
			break;
		case 'D':
		case 'd':
			selectCurrentDomain();
			break;
		case 'E':
		case 'e':
			showpar[showBonds] = showpar[showBonds]==0?1:0;
			printf("Show Bonds = %d\n",showpar[showBonds]);
			break;
		case 'C':
		case 'c':
			showpar[showCellCenters] = showpar[showCellCenters]==0?1:0;
			printf("Show Cell Centers = %d\n",showpar[showCellCenters]);
			break;
		case 'G':
			showpar[showGrid] = showpar[showGrid]==0?1:0;
			printf("Show Grid = %d\n",showpar[showGrid]);
			break;
		case 'g':
			showpar[showBoundaryGrid] = showpar[showBoundaryGrid]==0?1:0;
			printf("Show Boundary Grid = %d\n",showpar[showBoundaryGrid]);
			break;
		case 'F':
			showpar[showFaceCenters] = showpar[showFaceCenters]==0?1:0;
			printf("Show Face Centers = %d\n",showpar[showFaceCenters]);
			break;
		case 'f':
			showpar[showFrame] = showpar[showFrame]==0?1:0;
			printf("Show Frame = %d\n",showpar[showFrame]);
			break;
		case 'I':
#ifdef WTC
			{using namespace Solver;
				veldump=veldump<SMALL?4.0:0.0;
				printf("veldump=%g\n",veldump);
			}
			break;
#endif
		case 'i':
			init();
			break;
		case 'M':
		case 'm':
			consoleMenu();
			break;
		case 'N':
			showpar[showNodes] = showpar[showNodes]==0?1:0;
			printf("Show Grid Nodes = %d\n",showpar[showNodes]);
			//glutPostRedisplay();
			break;
		case 'n':
			showpar[showBoundaryVertexes] = showpar[showBoundaryVertexes]==0?1:0;
			printf("Show Boundary Vertexes = %d\n",showpar[showBoundaryVertexes]);
			break;
		case 'l':
		case 'L':
			readconf();
///	#ifndef GRIDGEN
///				if (P.iv>=0&&P.iv<P.nv)
///					printf("Particle.variable[%d]=%s\n",P.iv,P.var[P.iv].name[0]);
///				if (G.ivar>=0&&G.ivar<G.nvar)
///					printf("Grid.variable[%d]=%s\n",G.ivar,G.var[G.ivar].name[0]);
///	#endif
			break;
		case 'P':
		case 'p':
			showpar[showParticles] = showpar[showParticles]==0?1:0;
			printf("Show particles = %d\n",showpar[showParticles]);
			if (showpar[showParticles]==1)
			{
				//printf("No particles yet\n");P.n=0;
				domain_root[idomain].getPoints
				(	P.m,P.n,P.origin,P.first,P.last,P.X);
				printf("Number of particles = %d\n",P.n);
			}
			break;
		case 'R':
		case 'r':
			if (finished) break;
#ifdef WITH_MPI
			if(firstR)
			{
				firstR=0;
				for(int j=1; j<COMM_WORLD.Get_size(); j++)
				{
					COMM_WORLD.Send(vars, ndomains, MPI_INT, j, 100);
					COMM_WORLD.Send(coms, ndomains, MPI_INT, j, 101);
				}
				COMM_WORLD.Barrier();
			}
#endif
				
			animation=animation==0?1:0;
			if (option.verbose|option.debug)
				printf("animation=%d\n",animation);
			//	if (animation) 
			//		glutIdleFunc(animate);
			//	else
			//	{	glutIdleFunc(NULL);
			//		printf("\n");fflush(stdout);
			//	}

			break;
		case 'S':
		case 's':
			if (finished) break;
			printf("Time = ");fflush(stdout);
			if (run(1)==0) 
			{	finish();
				finished=1;
			}
			///refresh();
			printf("%g\n",runtime.current);fflush(stdout);
			break;
		case 'T':
			showpar[showTool] = showpar[showTool]==0?1:0;
			//if (showpar[showVariables]==1&&ivar<0)selectVariable();
			printf("Show Tool = %d\n",showpar[showTool]);
			break;
		case 't':
			showpar[showToolVertexes] = showpar[showToolVertexes]==0?1:0;
			printf("Show Tool Vertexes = %d\n",showpar[showToolVertexes]);
			break;
		case 'O':
		case 'o':
			toggle();
			break;
		case 'V':
		case 'v':
			showpar[showVariables] = showpar[showVariables]==0?1:0;
			//if (showpar[showVariables]==1&&ivar<0)selectVariable();
			printf("Show Variables = %d\n",showpar[showVariables]);
			break;
		case 'W':
			domain_root[idomain].saveGeom(programname);
			break;
		case 'w':
			dumpwindow();
		break;
		case ':':
		case '.':
			commandMode();
			break;
		case 27:
			Exit();
			break;
		case '?':
			helpDisplay();
			break;
		case 'Z': //added by zhang, used to show existing data as vectors
			char infilename[30];
			for(int i =0; i<ndomains; i++)
			{
				Domain *dom = domain_root + i;
				int ivar = dom->getDisplayVariable();
				sprintf(infilename, "%s_%s", dom->name, dom->variable[ivar].name);
				getDispVarDefault(dom, ivar, infilename);
			}
			break;
	}
}
static int
events_loop(Display *dpy, Window win) 
{	XEvent event;
	do 
	{	char buf[31];
		KeySym keysym;
		XNextEvent(dpy, &event);
		switch(event.type) 
		{	case Expose:
				break;
			case ConfigureNotify: 
				reshape(event.xconfigure.width, event.xconfigure.height);
			//	{	/* this approach preserves a 1:1 viewport aspect ratio */
			//		int vX, vY, vW, vH;
			//		int eW = event.xconfigure.width, eH = event.xconfigure.height;
			//		if (eW >= eH) 
			//		{
			//			vX = 0;
			//			vY = (eH - eW) >> 1;
			//			vW = vH = eW;
			//		} 
			//		else 
			//		{
			//			vX = (eW - eH) >> 1;
			//			vY = 0;
			//			vW = vH = eH;
			//		}
			//		glViewport(vX, vY, vW, vH);
			//	}
				break;
			case KeyPress:
				(void) XLookupString(&event.xkey, buf, sizeof(buf), &keysym, NULL);
				keyboard((unsigned char)keysym);
				break;
			case ButtonPress:
				mouse(event.xbutton.button,GLUT_DOWN,event.xbutton.x,event.xbutton.y);
				break;
			case MotionNotify:
				motion(event.xbutton.x, event.xbutton.y);
			break;
			default:
			break;
		}
/*
#ifdef WITH_MPI
	//post a recv buffer to receive command from process 0
		if(iproc!=0)
		{
			if (isRunRecved==0)
			{
				COMM_WORLD.Recv(&key_recved, 1 , MPI_UNSIGNED, 0, 1);
				keyboard(key_recved);
				if((key_recved ==82)||(key_recved ==114)) //'R'=82;'r'=114
					isRunRecved = 1;
			}
		}
#endif
*/
		if (animation)
		do
		{	animate();
			//display();
		}	while(!XPending(dpy));
	} while (XPending(dpy));
	display();
	return 0;
}
void	guirun(int argc, char* argv[])
{
	vars = new int[ndomains];
	coms = new int[ndomains];
	XVisualInfo *vi;
	XSetWindowAttributes swa;
	GLXContext cx;
	glutInit(&argc,argv);
	sprintf(windowname,"%s",argv[0]);
	initdisp();
	dpy = XOpenDisplay(0);
	if (!dpy) ERROR("CAN'T OPEN DISPLAY");
	vi = glXChooseVisual(dpy, DefaultScreen(dpy), attributeList);
	if (!vi) ERROR("NO SUITABLE VISUAL");
	cx = glXCreateContext(dpy, vi, 0, GL_TRUE);
	swa.colormap = XCreateColormap
	(	dpy, 
		RootWindow(dpy, vi->screen),
		vi->visual, AllocNone
	);
	swa.border_pixel = 0;
	swa.event_mask = ExposureMask | StructureNotifyMask | KeyPressMask |
	ButtonPressMask | ButtonMotionMask;
	if(option.debug)
		printf
		(	"Open window: width=%d, height=%d\n",
			window.width,window.height
		);
	if(option.debug)
		printf("Window: width=%d, height=%d\n",window.width, window.height);
	win=XCreateWindow
	(	dpy, 
		RootWindow(dpy, vi->screen), 
		0, 0, window.width, window.height,
		0, vi->depth, InputOutput, vi->visual,
		CWBorderPixel|CWColormap|CWEventMask, &swa
	);
	XStoreName(dpy, win, windowname);
	XMapWindow(dpy, win);
	glXMakeCurrent(dpy, win, cx);
	/* check for the polygon offset extension */
#ifndef GL_VERSION_1_1
	if (!query_extension("GL_EXT_polygon_offset"))
		error(argv[0], "polygon_offset extension is not available");
#else
	(void) query_extension;
#endif
	myinit(argc, argv);
	helpDisplay();
	while (events_loop(dpy, win)==0);
}
void displaymessage(char *msg)
{
	puts(msg);
}
} //END NAMESPACE
