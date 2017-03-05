#define REAL	double
//#define LOCALDIFFUSION

namespace ModGasFlow
{
extern double
		Rgas,//Gas-constant of air [work/(mass*temp)
		visc,//Viscosity, [length^2/time]
		Pamb,//Ambient pressure [mass/(length*time^2)]
		Tamb,//Ambient temperature
		Pin,//Inlet pressure
		Tin,//Inlet temperature
		Cwall;//Wall-friction coefficient [mass/(length^2*time)]
void	AdvanceFlow(double denin, double denout, double dt,Domain *dom);
double	InletDensity(double time);
double	OutletDensity(double time);
double	Pressure(double	den);
double	Pressure(double	den, double temp);
}
