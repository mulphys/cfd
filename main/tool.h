
enum Shapes
{	half_plane=0,
	inside_sphere,
	outside_sphere,
	inside_ellipsoid,
	outside_ellipsoid,
	maxshapes
};
struct	Frame
{	ElementStatus	type;//used to set elment types
	double
		radius,
		refinement,
		x[DIM];
	Frame	*next;//set of frmanodes forms a ring
};
class	Tool
{	int
		ivolume,//index to the volume varialbe 
		ixold,//index to old x-coordinate values
		nbn;//number of boundary vertexes
	double	surface_roughness,
				refinement;
	ElementStatus	etype;//used to set elment types
	public:
	enum ToolType
	{	coral=0,// cell-growth model
		crystal,// crystall-growth model
		maxtooltypes
	}	type;
	enum Shapes	shape;
	enum States
	{	passive,
		add_cells,
		remove_cells,
		push_cells
	}	state;
	struct
	{	unsigned int	active:	1;//1: tool is active 
		unsigned int	manual:	1;//1: tool is operated manually
	           //0: tool is attached to the frame
///		unsigned int	construct:	1; //1 - constructing grid
///		unsigned int	destroy:	1; //1 - destroing grid
		unsigned int	setnodes:	1; //1 - assigning node types
		unsigned int	setcells:	1; //1 - assigning cell types
		unsigned int	setboundary:	1; //1 - assigning boundary conditions 
	}	mode;
	struct	NodeCellList
	{	DNode	*node;
		DCellList	*neib_cell_root;
		NodeCellList	*next;
	}	*node_cell_root;
	double
		*rootpos,
		center[DIM],
		tip[DIM];//position of the tip: 
		         // used in destructive mode
	char	*configfile;
	Sphere	sphere;
	DNodeList 
		*first_relax_dnode,//all nodes inside the relaxation radius
		                   // used in volume relaxation scheme
		*first_tool_boundary_dnode,//all boundary nodes inside the tool radius:
		                      // used to advance the tool boundary by
									// computing repulsive force from the boundary
									// nodes, inside tool nodes and the tool
		*first_tool_dnode;//all internal nodes inside the tool radius
		                  // used to compute the force on boundary nodes
	DCell *root_cell;//cell containing a critical point (center or tip)
	DCellList
		*first_tool_boundary_dcell,//boundary cells inside the tool radius:
		                      // used in cell splitting scheme
		*first_relax_dcell;//all cells inside the relaxation radius:
		                   // used in volume relaxation scheme
	Frame	*first_framenode,*last_framenode,*current_framenode;
	Tool	*next;//a set of tools should form a ring
	Tool();
	Tool(char *name);
	~Tool();
	void	load();
	void	addnode(double *x, int type, double radius, double refinement);
	void	setype(ElementStatus type);
	void	setdefaults();
	void	pos(double x[]);
	void	getpos(double *x);
	void	move(double dx[]);
	void	setpos(double x, double y, double z);
	void	setradius(double r);
	void	setroughness(double r);
	void	setrefinement(double r);
	void	setVolInd(int ivol);
	void	setXoldInd(int ix);
	void	setboundary();
	double	getradius();
	double	getinsideradius();
	double	getrelaxradius();
	double	getlengthscale();
	double	force_neib(double distance, double interaction_radius);
	double	force_tool(double distance, double interaction_radius);
	int	createToolNodeList(Domain *D);
	int	createLists(Domain *D);
	int	updateLists(Domain *D);
	void	cleanSurface();
	int	growCoral(Domain *domain);
	int	growCrystal(Domain *domain);
	void	deleteLists();
	void	raid(DCell *cell);
	void	relax(Domain *D, double relaxation_factor);
	void	push();
	void	advanceBoundaryNodes();
	void	destroy(Domain *D);
	///template <class Element>
	///void	deleteElement(Element *&element);
	void	insert(DCell	*newcell);
	template	<class List>
	void	deleteList(List *&element);
	void	unlinkNodeListPointer(DNode *kill, DNodeList *root);
	void	unlinkCellListPointer(DCell *kill, DCellList *root);
};
class	Bio:Tool
{

};
class	Crystal:Tool
{

};

