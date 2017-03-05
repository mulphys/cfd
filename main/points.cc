/*
#include <stdlib.h>
#include <stdio.h>
#include <iostream.h>
#include <string.h>
#include <math.h>
*/
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <math.h>
#include "main.h"
#include "io.h"
#include "vecalg.h"
#include "geom.h"
#include "var.h"
//#include "particles.h"
#include "tool.h"
#include "domain.h"
//#include "solver.h"
#include "force.h"

using namespace std;

//using namespace Solver;
//using namespace Force;

void	Domain::injectParticlesAtFrame
(	int ninj //number of particles to inject
)
{	//Release the boundary particles
}
void	Domain::injectParticles
(	int ninj //number of particles to inject
)
{	//Release the boundary particles
}
void	Domain::initp()
{
	if(option.verbose)
	{	printf("Max number of particles: mp=%d\n",mp);fflush(stdout);
	}
	for (int i=0; i<mp; i++)
	{	Point	*p=origin+i;
		origin[i].next=origin+(i+1)%mp;
		origin[i].prev=origin+(i+mp-1)%mp;
	}
	np=0;
	first=origin;
	last=origin+mp-1;
}
void	Domain::addp()
{
	//	if (np==mp)
	//	{	if(option.verbose)
	//			cerr << "Can't add more particles\n"<<flush;
	//		return;
	//	}
	last=last->next;
	if (np++>0&&last==first)
	{
		first=first->next; np--;
	}
	if (np>mp||np<0)
	{
		cerr << "ERROR:Domain::addp: mp=" << mp << ", np="<<np<<"\n";
		exit(1);
	}
}
int	Domain::delp
(
	Point	*&p
)
{	if (np<=0) return 0;
	np--;
	if (p==first) 
		first=first->next;
	else
	if (p==last) 
	{	last=last->prev;
		p=last;
	}
	else
	{	Point *q;
		q=p;
		p=p->prev;
		p->next=q->next;
		q->next->prev=p;
		q->next=first;
		q->prev=first->prev;
		first->prev->next=q;
		first->prev=q;
	}
	return 1;
}
int	Domain::putp
(//Put point at a certain posisition 
	double	*x //coordinates
)
{	int	ip;
	double	*y;
	if (np==mp)
	{
		if(option.verbose) 
			cerr<<"WARNING: Can't add more particles: np="<<np<<", mp="<<mp<<endl;
		return -1;
	}
	ip=last->next-origin;
	y=coordinates[points].val+ip*DIM;
	for (int i=0; i<DIM; i++) y[i]=x[i];
//		variable[iscl].val[ip]=scl;
	addp();
	return ip;
}
void	Domain::setsclp
(//Initialize a single scalar variable
	int	ip,  //index to the particle
	int	iscl,//index to the scalar variable
	double	scl//initialization value for the scalar
)
{
	if (ip<0||ip>=mp)
	{
		if(option.verbose) 
			cerr<<"WARNING: Can't set scalar for particle: ip="<<ip<<", mp="<<mp<<endl;
		return;
	}
	variable[iscl].val[ip]=scl;
}
void	Domain::setvecp
(//Initialize a single vector variable
	int	ip,  //index to the particle
	int	ivec,//index to the vector variable
	double	*vec//initialization value for the vector
)
{	double	*v;
	if (ip<0||ip>=mp)
	{
		if(option.verbose) 
			cerr<<"WARNING: Can't set scalar for particle: ip="<<ip<<", mp="<<mp<<endl;
		return;
	}
	v=variable[ivec].val+DIM*ip;
	for (int i=0; i<DIM; i++)
		v[i]=vec[i];
}
int	Domain::injp
(//Inject a point at a certain posisition 
 // with a given velocity and initialize
 // a single scalar variable
	int	nvar,//number of initialized variables
	int	*ivar,//indexes to the variables to be initialized
	double	*var,//variables initialization array
	double	*x //coordinates
)//RETURNS: the number of injected particle
{	int	ip;
	double	*y,*v;
	if (np==mp)
	{
		if(option.verbose) 
			cerr<<"WARNING: Can't add more particles: np="<<np<<", mp="<<mp<<endl;
		return -1;
	}
	ip=last->next-origin;
	//Initializing coordinates
	y=coordinates[points].val+ip*DIM;
	for (int i=0; i<DIM; i++) y[i]=x[i];
	//Initializing variables
	v=var;
	for (int i=0; i<nvar; i++)
	{	int
			iv=ivar[i],
			dim=variable[iv].dimension;
		for (int k=0; k<dim; k++)
			variable[iv].val[ip+k]=v[k];
		v+=dim;
	}
	addp();
	return ip;
}
void	Domain::getPoints
(
	int &m,
	int	&n,
	struct Point *&origin0,
	struct Point *&first0,
	struct Point *&last0,
	double	*&X
)
{
	m=mp;
	n=np;
	X=coordinates[points].val;
	first0=first;
	last0=last;
	origin0=origin;
}
int	Domain::getNoPoints()
{
	return np;
}
int	Domain::getMaxNoPoints()
{
	return mp;
}
