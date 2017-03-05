namespace ModPoisson
{
	void	setDiagCoeff
	(	//Sets the diagonal coefficients c1ii
		//for the stiffness matrix \ref{}
		int	ivolc,
		int	idiag,
		DNode	*node_root,
		DCell	*cell_root
	);
/* 
 * Assemble source terms for
 * FEM discretized Poissson equation
 *
 */
	void	assemblePoissonSource
	(
		int	iscr,// index to the source variable to be assembled
		int	ifor,// index to the force variable 
		int	ivolc,//index to the cell-volume variables
		DNode	*node_root,
		DCell	*cell_root
	);
	void	assemblePoissonSourceGrad
	(
		int	iscr,// index to the source variable to be assembled
		int	ifor,// index to the force variable 
		int	ivolc,//index to the cell-volume variables
		DNode	*node_root,
		DCell	*cell_root
	);
	void	assemblePoissonSourceNeuman
	// Is called after assemblePoisssonSource
	// if the non-zero gradinet boundary conditions
	// are present
	(
		int	isrc,// index to the boundary source variable to be assembled
		int	ifor,// index to the force variable 
		int	ibnd,// index to the boundary vector variable 
//-		int	idiag,// boundary diagonal coeffs
		int	ivolc,//index to the cell-volume variables
		BFaceList	*bface_root,
		DNode	*node_root,
		DCell	*cell_root
	);
/*
 * One step of Gauss-Seidel iteration of a Poisson FEM solver
 *
 */
	void	stepPoisson
	(
		int	iscl,//index to the variable to be updated
		int	isrc,//source term variable-index
		int	idif,//diffusion variable-index
		int	idiag,//diagonal matrix elements
		int	ivoln,
		int	ivolc,
		DNode	*node_root,
		DCell	*cell_root
	);
};
