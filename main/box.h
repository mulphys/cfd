
class	Box
{
	public:
	unsigned int	bc[6];
	int	nx[DIM];
	void	create(double xmin[], double xmax[], double *X, struct Node *n);
	void	create(double xmin[], double xmax[], struct DNode *n);
	int	setBoundaryVertexes(struct Node *n);
	void	connect(struct Cell *c);
	void	connect(struct DNode *n, struct DCell *c);
};

