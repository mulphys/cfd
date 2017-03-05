extern	int	getNvarModCoral();
extern	void	defVarModCoral
(
	int	ivar,
	char	*name,
	int	&type,
	int	&rank,
	Domain*	D
);
extern	void	initVarModCoral
(	int	ivar,
	char	*name,
	int	type,
	int	size,
	int	rank,
	double	*val,
	Domain	*dom
);
extern	void	initModCoral(Domain*);
extern	void	stepModCoral(double dt, Domain *d);
//----------------------------------------------------
extern	int	getNvarModGasFlow();
extern	void	defVarModGasFlow
(
	int	ivar,
	char	*name,
	int	&type,
	int	&rank,
	Domain	*D
);
extern	void	initVarModGasFlow
(	int	ivar,
	char	*name,
	int	type,
	int	size,
	int	rank,
	double	*val,
	Domain	*dom
);
extern	void	initModGasFlow(Domain*);
extern	void	stepModGasFlow(double dt,Domain *d);
//----------------------------------------------------
extern	int	getNvarModMembrane();
extern	void	defVarModMembrane
(
	int	ivar,
	char	*name,
	int	&type,
	int	&rank,
	Domain	*D
);
extern	void	initVarModMembrane
(	int	ivar,
	char	*name,
	int	type,
	int	size,
	int	rank,
	double	*val,
	Domain	*dom
);
extern	void	initMembrane(Domain*);
extern	void	stepMembrane(double dt,Domain *d);
extern	void	initMicrobe(Domain*);
extern	void	stepMicrobe(double dt,Domain *d);
extern	void	initBreath(Domain*);
extern	void	stepBreath(double dt,Domain *d);
//----------------------------------------------------
extern	int	getNvarModElastic();
extern	void	defVarModElastic
(
	int	ivar,
	char	*name,
	int	&type,
	int	&rank,
	Domain	*D
);
extern	void	initVarModElastic
(	int	ivar,
	char	*name,
	int	type,
	int	size,
	int	rank,
	double	*val,
	Domain	*dom
);
extern	void	initElastic(Domain*);
extern	void	stepElastic(double dt,Domain *d);
//----------------------------------------------------
extern	int	getNvarModLattice();
extern	void	defVarModLattice
(
	int	ivar,
	char	*name,
	int	&type,
	int	&rank,
	Domain	*D
);
extern	void	initVarModLattice
(	int	ivar,
	char	*name,
	int	type,
	int	size,
	int	rank,
	double	*val,
	Domain	*dom
);
extern	void	initModLattice(Domain*);
extern	void	stepModLattice(double dt,Domain *d);
//----------------------------------------------------
extern	int	getNvarModBonding();
extern	void	defVarModBonding
(
	int	ivar,
	char	*name,
	int	&type,
	int	&rank,
	Domain	*D
);
extern	void	initVarModBonding
(	int	ivar,
	char	*name,
	int	type,
	int	size,
	int	rank,
	double	*val,
	Domain	*dom
);
extern	void	initModBonding(Domain*);
extern	void	stepModBonding(double dt,Domain *d);
//----------------------------------------------------
extern	int	getNvarModCollapse();
extern	void	defVarModCollapse
(
	int	ivar,
	char	*name,
	int	&type,
	int	&rank,
	Domain	*D
);
extern	void	initVarModCollapse
(	int	ivar,
	char	*name,
	int	type,
	int	size,
	int	rank,
	double	*val,
	Domain	*dom
);
extern	void	initModCollapse(Domain*);
extern	void	stepModCollapse(double dt,Domain *d);
//----------------------------------------------------
extern	int	getNvarModPoisson();
extern	void	defVarModPoisson
(
	int	ivar,
	char	*name,
	int	&type,
	int	&rank,
	Domain	*D
);
extern	void	initVarModPoisson
(	int	ivar,
	char	*name,
	int	type,
	int	size,
	int	rank,
	double	*val,
	Domain	*dom
);
extern	void	initModPoisson(Domain*);
extern	void	stepModPoisson(double dt,Domain *d);
//----------------------------------------------------
extern	int	getNvarModIncomFlow();
extern	void	defVarModIncomFlow
(
	int	ivar,
	char	*name,
	int	&type,
	int	&rank,
	Domain	*D
);
extern	void	initVarModIncomFlow
(	int	ivar,
	char	*name,
	int	type,
	int	size,
	int	rank,
	double	*val,
	Domain	*dom
);
extern	void	initModIncomFlow(Domain*);
extern	void	stepModIncomFlow(double dt,Domain *d);
//-------------------------------------------------
extern	int	getNvarModParticles();
extern	void	defVarModParticles
(
	int	ivar,
	char	*name,
	int	&type,
	int	&rank,
	Domain	*D
);
extern	void	initVarModParticles
(	int	ivar,
	char	*name,
	int	type,
	int	size,
	int	rank,
	double	*val,
	Domain	*dom
);
extern	void	initModParticles(Domain*);
extern	void	stepModParticles(double dt,Domain *d);
//-------------------------------------------------
extern	int	getNvarModCrystal();
extern	void	defVarModCrystal
(
	int	ivar,
	char	*name,
	int	&type,
	int	&rank,
	Domain	*D
);
extern	void	initVarModCrystal
(	int	ivar,
	char	*name,
	int	type,
	int	size,
	int	rank,
	double	*val,
	Domain	*dom
);
extern	void	initModCrystal(Domain*);
extern	void	stepModCrystal(double dt,Domain *d);
//-------------------------------------------------
extern	int	getNvarModRFG();
extern	void	defVarModRFG
(
	int	ivar,
	char	*name,
	int	&type,
	int	&rank,
	Domain	*D
);
extern	void	initVarModRFG
(	int	ivar,
	char	*name,
	int	type,
	int	size,
	int	rank,
	double	*val,
	Domain	*dom
);
extern	void	initModRFG(Domain*);
extern	void	stepModRFG(double dt,Domain *d);
//-------------------------------------------------
extern	int	getNvarModParticleFlow();
extern	void	defVarModParticleFlow
(
	int	ivar,
	char	*name,
	int	&type,
	int	&rank,
	Domain	*D
);
extern	void	initVarModParticleFlow
(	int	ivar,
	char	*name,
	int	type,
	int	size,
	int	rank,
	double	*val,
	Domain	*dom
);
extern	void	initModParticleFlow(Domain*);
extern	void	stepModParticleFlow(double dt,Domain *d);
//-------------------------------------------------
extern	int	getNvarModPRFG();
extern	void	defVarModPRFG
(
	int	ivar,
	char	*name,
	int	&type,
	int	&rank,
	Domain	*D
);
extern	void	initVarModPRFG
(	int	ivar,
	char	*name,
	int	type,
	int	size,
	int	rank,
	double	*val,
	Domain	*dom
);
extern	void	initModPRFG(Domain*);
extern	void	stepModPRFG(double dt,Domain *d);
//-------------------------------------------------
//- extern	int	getNvarModPStatic();
//- extern	void	defVarModPStatic
//- (
//- 	int	ivar,
//- 	char	*name,
//- 	int	&type,
//- 	int	&rank,
//- 	Domain	*D
//- );
//- extern	void	initVarModPStatic
//- (	int	ivar,
//- 	char	*name,
//- 	int	type,
//- 	int	size,
//- 	int	rank,
//- 	double	*val,
//- 	Domain	*dom
//- );
//- extern	void	initModPStatic(Domain*);
//- extern	void	stepModPStatic(double dt,Domain *d);
