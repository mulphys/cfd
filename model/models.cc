/* 
 *  MODEL: Viscous Flow
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
#include "var.h"
#include "tool.h"
#include "domain.h"
#include "models.h"

void	Domain::setModel()
{
	if(option.debug){printf("Setting model %d\n",model);FLUSH;}
	switch(model)
	{
		case 0:
			nvar=getNvarModCoral();
			defvar=defVarModCoral;
			initvar=initVarModCoral;
			init=initModCoral;
			step=stepModCoral;
			break;
		case 1:
			nvar=getNvarModGasFlow();
			defvar=defVarModGasFlow;
			initvar=initVarModGasFlow;
			init=initModGasFlow;
			step=stepModGasFlow;
			break;
		case 2:
			nvar=getNvarModMembrane();
			defvar=defVarModMembrane;
			initvar=initVarModMembrane;
			init=initMicrobe;
			step=stepMicrobe;
			break;
		case 3:
			nvar=getNvarModLattice();
			defvar=defVarModLattice;
			initvar=initVarModLattice;
			init=initModLattice;
			step=stepModLattice;
			break;
		case 4:
			ERROR("MODEL NOT DONE YET");
			nvar=getNvarModBonding();
			defvar=defVarModBonding;
			initvar=initVarModBonding;
			init=initModBonding;
			step=stepModBonding;
			break;
		case 5:
			ERROR("MODEL NOT DONE YET");
			nvar=getNvarModCollapse();
			defvar=defVarModCollapse;
			initvar=initVarModCollapse;
			init=initModCollapse;
			step=stepModCollapse;
			break;
		case 6:
			nvar=getNvarModPoisson();
			defvar=defVarModPoisson;
			initvar=initVarModPoisson;
			init=initModPoisson;
			step=stepModPoisson;
			break;
		case 7:
			nvar=getNvarModIncomFlow();
			defvar=defVarModIncomFlow;
			initvar=initVarModIncomFlow;
			init=initModIncomFlow;
			step=stepModIncomFlow;
			break;
		case 8:
			nvar=getNvarModParticles();
			defvar=defVarModParticles;
			initvar=initVarModParticles;
			init=initModParticles;
			step=stepModParticles;
			break;
		case 9:
			nvar=getNvarModCrystal();
			defvar=defVarModCrystal;
			initvar=initVarModCrystal;
			init=initModCrystal;
			step=stepModCrystal;
			break;
		case 10:
			nvar=getNvarModRFG();
			defvar=defVarModRFG;
			initvar=initVarModRFG;
			init=initModRFG;
			step=stepModRFG;
			break;
		case 11:
			nvar=getNvarModParticleFlow();
			defvar=defVarModParticleFlow;
			initvar=initVarModParticleFlow;
			init=initModParticleFlow;
			step=stepModParticleFlow;
			break;
		case 12:
			nvar=getNvarModPRFG();
			defvar=defVarModPRFG;
			initvar=initVarModPRFG;
			init=initModPRFG;
			step=stepModPRFG;
			break;
		case 13:
			nvar=getNvarModMembrane();
			defvar=defVarModMembrane;
			initvar=initVarModMembrane;
			init=initBreath;
			step=stepBreath;
			break;
		case 14:
//-			nvar=getNvarModPStatic();
//-			defvar=defVarModPStatic;
//-			initvar=initVarModPStatic;
//-			init=initModPStatic;
//-			step=stepModPStatic;
//-			break;
//-		case 15:
//?			nvar=getNvarModMembrane();
//?			defvar=defVarModMembrane;
//?			initvar=initVarModMembrane;
//?			init=initMembrane;
//?			step=stepMembrane;
			nvar=getNvarModElastic();
			defvar=defVarModElastic;
			initvar=initVarModElastic;
			init=initElastic;
			step=stepElastic;
			break;
		default:
			fprintf(stderr,"Model %d not implemented\n",model);
			exit(1);
	}
}
