#include "main.h"

namespace Mod1
{
	enum VariablesMod1
	{
		pressure,
		pressure_sourceA,
		pressure_sourceB,
		velocity,
		boundary_velocity,
		boundary_scalar,
		boundary_source,
		convection,
		diffusion,
		maxVarMod1
	};
};

using namespace Mod1;

int	getNvarMod1()
{
	return maxVarMod1;
}
void	defVarMod1
(
	int	ivar,
	char	*name,
	int	&type,
	int	&rank,
	Domain	*D
)
{

}
void	initVarMod1
(	int	ivar,
	char	*name,
	int	type,
	int	size,
	int	rank,
	double	*val,
	Domain	*dom
)
{

}
void	initMod1(Domain*)
{

}
void	stepMod1(double dt,Domain *d)
{

}
