/*
	Number of cells:
	Core: 3*6=18
	Skirt (first layer): 2*6+2*6+6=30
	...
	plane layer: 
*/

//Number of Nodes in the Ring number n:
#define	NNR(n)	(n>=0?3*(n+1)*(n)+1:0)
//Number of Trapeziods in the Ring number n:
#define	NTD(n)	6*(n)*(n)
//Number of Cells in the Ring number n:
#define	NCR(n)	18*(n)*(n)
//Index of the cell
#define	NB(itrapez,itetra)	(itrapez<0||itrapez>=nt?boundary:3*(itrapez)+(itetra))

class	Cylinder
{
	unsigned int	bc[3];
	int	//number of
		nc, //cells
		nn, //nodes
		np,	//planes in the cylinder
		ns,	//slices in the cylinder=np-1
			//slice=a disk region between two cross-sectional planes
		nt,	//trapezoids in the cylinder
		ni,	//circles in the cylinder=nr+1
		nr,	//rings in the cylinder=nd-1
		nnp,	//nodes on the plane
		ntr,	//trapezoids in the ring
		nts,	//trapezoids in the slice 
		ncs,	//cells in the slice
		nx[2];//nx[0]: number of rings; nx[1]: number of slices
	double	R,L;
	public:
	~Cylinder();
//?	void	ring(int i, int ip, double *X, Cell *C);
	void	ring(int i, int ip, DNode **I, DCell **C);
	void	init(int	nx[], double dd[]);
//?	void	create(double *X, Node *n, Cell *c);
	void	create(DNode *n, DCell *c);
//?/	int	setBoundaryVertexes(Node *n);
//?/	int	setBoundaryVertexes(DNode *n);
//?/	void	connect(struct Cell *C);
//?/	void	connect(DCell *C);
	int	getnodes();
	int	getcells();
	void	setbc(int b[]);
	void	setgrid(int n[]);
};
