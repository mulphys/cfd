void	assemblePoissonNode
(
	int	nc,//number of cells
	int	nn,//number of nodes
	int	nbf,//number of boundary faces 
	int	nbv,//number of boundary vertexes
	struct Cell	*cell,
	struct Face	*face,
	struct Node	*node,
	struct BoundaryVertex	*Bv,
	double	*Vol,//cell-volumes array
	double	*S,//source terms
//	double	*B,//boundary source terms (nodal)
	double	*A //coefficients to be assembled
);
void	stepPoissonNodeGS
//One Gauss-Seidel iteration of a Poisson equation in FEM discretization
(
	int	nc,//number of cells
	int	nn,//number of nodes
	struct Cell	*cell,
	struct Face	*face,
	struct Node	*node,
	double	*Volc,//cell-volumes array
	double	*A,//Source terms (assembpled coeffs ../doc/fem.tex:\ref{poisiter})
	double	*B,//Boundary source terms: vector whose face-normal component
	         // determines a face-normal derivative of P
	//RETURNS:
	double	*P //Node-scalar variable to be solved
);
void	stepPoissonNodeGSNeuman
//One Gauss-Seidel iteration of a Poisson equation in FEM discretization
(
	int	nc,//number of cells
	int	nn,//number of nodes
	int	nbc,//number of boundary cells
	struct Cell	*cell,
	struct Face	*face,
	struct Node	*node,
	double	*Volc,//cell-volumes array
	double	*Voln,//node-volumes array
	double	*Y,//face-centers
	double	*Z,//cell-centers
	double	*A,//Source terms (assembpled coeffs ../doc/fem.tex:\ref{poisiter})
	double	*B,//Boundary source terms: vector whose face-normal component
	         // determines a face-normal derivative of P
	double	*C,//Neighbor's contributions
	//RETURNS:
	double	*P //Node-scalar variable to be solved
);
void	stepPoissonFaceGS
//One iteration of cell-centered poisson solver
(
	int	nf,
	int	nbf,
	struct Cell	*cell,
	struct Face	*face,
	double	*Volc,//cell-center volumes
	double	*Volf,//cell-center volumes
	double	*X,
	double	*Y, //face-coordinates
	double	*Z, //cell-coordinates
	double	*S, //source terms (face scalar)
	double	*P  //pressure (face scalar)
);
void	stepPoissonCellGS
(
	int	nc,//number of cells
	struct Cell	*cell,
	struct Face	*face,
	double	*Vol,//cell-volume
	double	*Y,//face-center nodes
	double	*Src, //source-term
	double	*Val//variable to solve 
);
void	gradNodeScl2CellVec
//Computes gradient of a nodal scalar at the cell centers
(
	int	nc,//number of cells
	struct Cell	*cell,
	struct Face	*face,
	double	*X,//nodal coordinates
	double	*S,//nodal scalar
	double	*G //cell-center gradient
);
void gradFaceScl2CellVec
//Computes gradient of a face-scalaer at the cell-centers
// delatation=(\nabla\,,scalar)
(
	int	nc,//number of cells
	struct Cell	*cell,
	struct Face	*face,
	double	*Vol,//cell volumes
	double	*P,//face-scalar variable which gradient is to be computed
	double	*G //cell-center vector representing the gradient of P
);
void gradCellScl2CellVec
//Computes gradient of a face-scalaer at the cell-centers
// delatation=(\nabla\,,scalar)
(
	int	nc,//number of cells
	struct Cell	*cell,
	struct Face	*face,
	double	*Vol,//cell volumes
	double	*P,//cell-scalar variable which gradient is to be computed
	double	*G //cell-center vector representing the gradient of P
);
void	fluxCellVec2NodeScl
//Computes delatation of a face-vector at the nodes
// delatation=(\nabla\,,variable)
(	
	int	nc,//number of cells
	int nn,//number of nodes
	struct Cell	*cell,
	struct Face	*face,
	struct Node	*node,
	struct BoundaryVertex	*Bv,
	double	*U,//cell-velocity to compute boundary convection terms
	double	*V,//face-vector variable which gradient is to be computed
	double	*D //scalar node-variable representing the delatation of V
);
void	fluxFaceVec2NodeScl
//Computes delatation of a face-vector at the nodes
// delatation=(\nabla\,,variable)
(	
	int	nc,//number of cells
	int nn,//number of nodes
	struct Cell	*cell,
	struct Face	*face,
	double	*V,//face-vector variable which gradient is to be computed
	double	*D //scalar node-variable representing the delatation of V
);
void gradFaceVec2NodeScl
//Computes delatation of a face-vector at the nodes
// delatation=(\nabla\,,variable)
(	
	int	nc,//number of cells
	int	nn,//number of nodes
	struct Cell	*cell,
	struct Face	*face,
	double	*Vol,//node-volumes
	double	*V,//face-vector variable which gradient is to be computed
	double	*D //scalar node-variable representing the delatation of V
);
void fluxFaceVec2CellScl
//Computes delatation of a face-vector at the cell-centers
// delatation=(\nabla\,,variable)
(
	int	nc,//number of cells
	struct Cell	*cell,
	struct Face	*face,
	double	*V,//face-vector variable which gradient is to be computed
	double	*D //scalar cell-variable representing the delatation of V
);
void	delatFaceVec2NodeScl
//Computes delatation of a face-vector at the nodes
// delatation=(\nabla\,,variable)
(	
	int	nc,//number of cells
	int	nn,//number of nodes
	struct Cell	*cell,
	struct Face	*face,
	double	*Vol,//node-based volume
	double	*V,//face-vector variable which gradient is to be computed
	double	*D //scalar node-variable representing the delatation of V
);
void delatFaceVec2CellScl
//Computes delatation of a face-vector at the cell-centers
// delatation=(\nabla\,,variable)
(
	int	nc,//number of cells
	struct Cell	*cell,
	struct Face	*face,
	double	*Vol,//cell volumers
	double	*V,//face-vector variable which gradient is to be computed
	double	*D //scalar cell-variable representing the delatation of V
);
void	delatCellVec2NodeScl
//Computes delatation of a face-vector at the nodes
// delatation=(\nabla\,,variable)
(	
	int	nc,//number of cells
	int	nn,//number of nodes
	int	nbc,//number of boundary cells
	struct Cell	*cell,
	struct Face	*face,
	struct Node	*node,
	struct BoundaryVertex	*Bv,
	double	*Vol,//node-based volume
	//double	*U,//cell-velocity for boundary convection
	double	*V,//face-vector variable which gradient is to be computed
	double	*D //scalar node-variable representing the delatation of V
);
void	delatCellVec2CellScl
//Computes delatation of a cell-vector at the cell centers
// delatation=(\nabla\,,variable)
(	
	int	nc,//number of cells
	int	nn,//number of nodes
	int	nbc,//number of boundary cells
	struct Cell	*cell,
	struct Face	*face,
	struct Node	*node,	
	double	*Y,//face-center coordinates
	double	*Z,//cell-center coordinates
	double	*Volc,//cell-based volume
	double	*U,//cell-velocity for boundary convection
	double	*V,//cell-vector variable which delatation is to be computed
	double	*W,//face-velocity at the boundary
	double	*D //scalar cell-variable representing the delatation of V
);
void	mapCellVec2NodeSclBnd
//Map cell-center vector boundary node scalar
//(used to compute the boundary source term in NS)
(
	int	nbf,//number of boundary faces
	int	nbn,//number of boundary nodes
	struct Cell	*cell,
	struct Face	*face,
	struct Node	*node,	
	struct BoundaryVertex	*Bv,
	//double	*Volc,//cell-based volume
	double	*Barea,//Boundary node-based volume (temporal variable)
	double	*V,//cell-vector
	double	*B //boundary node scalar representing the normal projection of V
);
void	interVectorCellFaceCentral
//Centeral-interpolates cell-centered vector to the faces
(
	int	nf, //number of faces
	int	nbf,//number of boundary faces
	struct Cell	*cell,
	struct Face	*face,
	double	*U, //cell-center vector
	double	*V  //face-interpolated variable
);
void	exterVectorCellFaceCentral
//Centeral-interpolates cell-centered vector to the faces
(
	int	nbf,//number of boundary faces
	struct Cell	*cell,
	struct Face	*face,
	double	*U, //cell-center vector
	double	*W, //boundary-face vector
	double	*V  //extrapolated face vector
);
void	mapVecCell2Face
//Centeral-interpolates cell-centered vector to the faces
(
	int	nf, //number of faces
	int	nbf,//number of boundary faces
	struct Cell	*cell,
	struct Face	*face,
	double	*U, //cell-center vector
//	double	*W, //boundary-face vector
	double	*V  //face-interpolated variable
);
void	mapVecCell2FaceUpwind
//Centeral-interpolates cell-centered vector to the faces
(
	int	nf, //number of faces
	int	nbf,//number of boundary faces
	struct Cell	*cell,
	struct Face	*face,
	double	*U, //cell-center vector
	double	*W, //boundary-face vector
	double	*V  //face-interpolated variable
);
//	void	interVectorCellFaceUp
//	//Upwind interpolates cell-centered vector to the faces
//	(
//		int	nc, //number of cells
//		int	nf, //number of faces
//		int	nbf,//number of boundary faces
//		struct Cell	*cell,
//		struct Face	*face,
//		double	*V, //cell-center vector
//		double	*U  //face-interpolated variable
//	);
void	conVecFace2Cell
(	//Convection of face variables to the cell-center
	int	nc,
	struct Cell	*cell,
	struct Face	*face,
	double	*Vol,//cell-volumes 
	double	*V,//face-variable
	double	*U,//face-velocity
	double	*C //cell-centered convection term
);
void	conCellVecUp
(	//Convection of face variables to the cell-center 
	//with upwind interpolation of face-varialbes
	int	nc,
	struct Cell	*cell,
	struct Face	*face,
	double	*Vol,//cell-volumes 
	double	*V,//cell-centered variable
	double	*U,//cell-centered velocity
	double	*W,//boundary velocity
	double	*C //cell-centered convection term
);
void	normgradCellVec2FaceScl
(	//Face-normal derivative of cell-center vector
	// (\nabla\,,vector)
	int	nf,
	int	nbf,
	struct Cell	*cell,
	struct Face	*face,
	double	*Y,//face-center coordinates
	double	*Z,//cell-center coordinates
	double	*V,//cell-center vector
	//double	*U,//boundar-face vector for V
	//RETURNS:
	double	*G //face-normal derivative of V (scalar: rank=0)
);
void	difCellVec
//Computes diffusion of a cell-centered vector
(
	int	nc, //number of cells
	int	nbc,//number of boundary cells
	struct Cell	*cell,
	struct Face	*face,
	double	*Vol,//cell-volume
	double	*Y,//face-center coordinates
	double	*Z,//cell-center coordinates
	double	*V,//cell-centered vector variable
	double	*W,//face-vector at the boundary
	double	*D //cell-centered diffusion term
);
void	addvec
(	//A=(A+B)
	int	n,
	double	*A,
	double	*B
);
void	subvec
(	//A=(A+B)*d
	int	n,
	double	*A,
	double	*B
);
void	addvec
(	//A=(A+B)*d
	int	n,
	double	*A,
	double	*B,
	double	d
);
void	subvec
(	//A=(A-B)*d
	int	n,
	double	*A,
	double	*B,
	double	d
);
void	div(int n, double *A, double *B);
void	zero(int n, double *A);
void	zerovec(int n, double *A);
