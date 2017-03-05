/* 
 *  MODEL: Generic Model Template
 *
 */
#include <stdlib.h>
#include <stdio.h>
//#include <iostream.h>
#include <string.h>
#include <math.h>
#include "main.h"
#include "io.h"
#include "vecalg.h"
#include "geom.h"
#include "func.h"
#include "var.h"
//#include "particles.h"
#include "tool.h"
#include "domain.h"


namespace Model
{
	enum Variables
	{
		density,
		pressure,
		velocity,
		convection,
		diffusion,
		maxVar
	};
};

using namespace Model;

int	getNvarMod0()
{
	return maxVar;
}
void	defVarMod0
(
	int	ivar,
	char	*name,
	int	&type,
	int	&rank,
	Domain	*dom
)
{
	int idomain=dom-domain_root;
	*name='\0';
	switch(ivar)
	{
		case pressure:
			strcpy(name,"pressure");
			type=nodes;
			rank=0;
			break;
		case density:
			strcpy(name,"density");
			type=nodes;
			rank=0;
			break;
		case velocity:
			strcpy(name,"velocity");
			type=cells;
			rank=1;
			break;
		case diffusion:
			strcpy(name,"diffusion");
			type=cells;
			rank=1;
			break;
		case convection:
			strcpy(name,"convection");
			type=cells;
			rank=1;
			break;
		default:
			fprintf
			(	stderr,
				"Specification of variable %d domain %d (%s) is incomplete\n",
				ivar+1,idomain+1,name
			);
			return;
	}
}
void	initVarMod0
(	int	ivar,
	char	*name,
	int	type,
	int	size,
	int	rank,
	double	*val,
	Domain	*dom
)
{	switch(ivar)
	{
		case pressure:
		break;
		default:
		;
	}
	if (option.debug)
	{	printf
		(	"Variable initialized\n"
		);FLUSH;
	}
}
void	initMod0(Domain*)
{

}
void	stepMod0(double dt,Domain *dom)
{
	time.current+=time.step;
}
