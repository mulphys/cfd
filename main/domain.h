//add a Ring class here
template <class Object>
class Ring;     // Incomplete declaration.

template <class Object>
class RingNode
{
	Object   element;
	RingNode *next, *prev;

	RingNode(const Object& theElement):element(theElement)
	{
		next = NULL;
		prev = NULL;
  	}

	friend class Ring<Object>;
};


template <class Object>
class Ring
{
private:
  	RingNode<Object> *header, *current;
  	int n;

public:
	Ring();
	~Ring( );

	bool isEmpty( ) const;

	//return the pointer pointing to the first, last or current
	RingNode<Object>* firstNode() const;
	RingNode<Object>* lastNode() const;
	RingNode<Object>* currentNode() const;

	//return the element value of the first , last or current
	Object firstElement() const;
	Object lastElement() const;
	Object currentElement() const;
	//return the next element of current, and move current to the next
	Object nextElement();

	//insert x before current node
	void insertBefore( const Object & x);
	
	//insert x after current node
	void insertAfter( const Object & x);

	//remove the current node
	void removeNode();

	//remove the first occurence of x
	void removeNode(const Object & x );

	//current advanced to next 
	void advance();
	
	//number of nodes of the ring
	int size() const;
	
	//make the ring empty
	void clear();

};

/**implementation of Ring class

 * Constructor, construct a empty ring.
 */
template <class Object>
Ring<Object>::Ring( )
{
	header = NULL;
	current = header;
	n = 0;
}

/**
 * Destructor
 */
template <class Object>
Ring<Object>::~Ring( )
{
	clear( );
	delete header;
	delete current;
}

/**
 * Test if the list is logically empty.
 * return true if empty, false otherwise.
 */
template <class Object>
bool Ring<Object>::isEmpty( ) const
{
	return n == 0;
}

/**
 * return the size of the ring
 */
template <class Object>
int Ring<Object>::size( ) const
{
	return n;
}

/**
 * Make the list logically empty.
 */
template <class Object>
void Ring<Object>::clear( )
{
	while(!isEmpty()) removeNode();
	delete header;
	delete current;
}

template <class Object>
RingNode<Object>* Ring<Object>::firstNode( ) const
{
	return header;
}

template <class Object>
RingNode<Object>* Ring<Object>::lastNode( ) const
{
	if(!isEmpty())
	{
		return header->prev;
	}
	else
	{
		return NULL;
	}
}

template <class Object>
RingNode<Object>* Ring<Object>::currentNode( ) const
{
	if(!isEmpty())
	{
		return current;
	}
	else
	{
		return NULL;
	}
}

template <class Object>
Object Ring<Object>::firstElement( ) const
{
	if(!isEmpty())
	{
		return header->element;  
	}
	else
	{
		throw "The ring is empty";
	}
}

template <class Object>
Object Ring<Object>::lastElement( ) const
{

	if(!isEmpty())
	{
		return header->prev->element;
	}
	else
	{
		throw "The ring is empty";
	}
}

template <class Object>
Object Ring<Object>::currentElement( ) const
{
	if(!isEmpty())
	{
		return current->element;
	}
	else
	{
		throw "The ring is empty";
	}
}

template <class Object>
Object Ring<Object>::nextElement( )
{
	if(!isEmpty())
	{
		current = current->next;
		return current->element;
	}
	else
	{
		throw "The ring is empty";
	}
}

template <class Object>
void Ring<Object>::advance( )
{
	if(!isEmpty())	current = current->next;
}

/**
 * Insert item x before current node
 */
template <class Object>
void Ring<Object>::insertBefore( const Object & x)
{
	RingNode<Object>* node = new RingNode<Object>(x);
	if(isEmpty())
	{
		header = node;
		header->prev = node;
		header ->next = node;
	}
	else
	{
		node->next = current;
		node->prev = current->prev;
		current->prev->next = node;
		current->prev=node;
	}
	/*point current to the new node */
	current = node;
	/*increase the size*/
	n++;
}

/*insert x after the current node */
template <class Object>
void Ring<Object>::insertAfter(const Object & x)
{
	RingNode<Object>* node = new RingNode<Object>(x);
	if(isEmpty())
	{
		header = node;
		header->prev = node;
		header ->next = node;
	}
	else
	{
		node->next = current->next;
		node->prev = current;
		current->next->prev = node;
		current->next=node;
	}
	/*point current to the new node */
	current = node;
	/*increase the size*/
	n++;
}

/**
 * Remove the current node
 */
template <class Object>
void Ring<Object>::removeNode()
{
	RingNode<Object> *tmp;
	if(size()>1)
	{
		tmp = current;
		current->prev->next = current->next;
		current->next->prev = current->prev;

		if(current==header) header = header->next;
		current = current->next;
		delete tmp;
		n--;
	}
	if(size()==1)
	{
		tmp = current;
		current = header = NULL;
		n--;
		delete tmp;
	}

}

/**
 * Remove the first occurrence of an item x.
 */
template <class Object>
void Ring<Object>::removeNode(const Object & x )
{

}

/**end of implementation of Ring class *****/

/***********************/

#ifdef WITH_DOVE
struct Dove
{
	int	
		n, //number of overlapping nodes of the parallel domain
		*I;//node-indexes of the prallel doman
	double	*X;//node-coordinates of the parallel domain
	DCell	**C;//cells of this domain containing X
	Domain	*dom;//parallel domain
	Dove	*next;
};
#endif
class	Domain
{
	int
		model,//number of the model used
		nvar,//number of variables
		nbv,//number of boundary vertexes
		nbf,//number of boundary faces
		nbc,//number of boundary cells
		np, //number of active points (particles)
		mp, //the max. number of points = Ne[points]
		level,//current nesting level
		fixedlimits,
		dispvar,//variable to display
		dispvarcomp;//component of the variable to display

	int	n_nodes; //number of nodes, added by zhang

	double
		xmin[DIM],xmax[DIM];
	char	output_type_name[MAXLINLEN];
	struct	Display
	{	int
			lighting;
		float
			rgbcolor[3];//0-red, 1-green, 2-blue'
	}	disp[maxelements];
	public:
	Tool *tool;
	char
		name[MAXLINLEN];
	int
#ifdef WITH_DOVE
		iproc,//process assigned to this domain
#endif
		type,
		//Me[maxtypes],//Ne[i]=maximum number of grid elements of type i=0:DIM
		Ne[maxelements];//Ne[i]=number of active grid elements of type i=0:DIM
	Variable
		*coordinates,//coordinates (vectors)
		*volumes,//volumes of the elements
		*variable;//a collection of variables
	Point	*origin,*first,*last;//a loop of particles
	DNode	*dnode_root;
	DCell	*dcell_root;
	BNodeList	*bnode_root;//list of boundary nodes
	BFaceList	*bface_root;//list of boundary faces
	DCellList	*bcell_root;//list of boundary cells 
#ifdef WITH_DOVE
	double	*Xb;//array of boundary coordinates

	/*replace *dove_root with *dove_ring */
	Ring<Dove*> *dove_ring;
	/*************************/

	DNode	**bnode;//bnode[0:nbv]: boundary node pointer array
	int	countBoundaryNodes();
	int indexBoundaryNodes();
	void storeBoundaryCoordinates();
	void setBoundaryFaceTypes();

	void ConnectThisProcDomains(int iproc);

	void sendBoundary(Domain *domi, Domain *domj);
	void receiveBoundary(Domain *domi, Domain *domj);

	void sendOverlap();
	void receiveOverlap();

	void sendVecArray(int n, double *X, Domain *dom);
	void receiveArraySize(Domain *dom);
	int receiveVecArraySize(Domain *D);
	void receiveVecArray(int n, double *X, Domain *D);

	void addDoveI(int nbn, double *X, Domain *D);
	void addDoveX(int nbn, double *X, Domain *D);

	void putVar(int ivar, Dove *dove);
	void sendVar(int ivar, Dove *dove);
	int receiveVar(int ivar, Dove *dove);

	void updateDove(int ivar);
	void updateDove();

	void sendDove(int ivar);
	void receiveDove(int ivar);

	void getArray(int idom, int n, double *X);
#endif //END WITH_DOVE

	BoundaryVertex	*Bv;
	void	(*init)(Domain *D);
	void	(*step)(double dt, Domain *D);
	Domain();
	~Domain();
	void setDomain();
	void setModel();
	void setGeom(char file[]);
	void setGeomCFG(char file[]);
	void setGeomXML(char file[]);
	void setCellVolumes(int ivar);
	void setNodeVolumes
	(
		int	ivarvoln,
		int	ivarvolc
	);
	void setNodeVolumes(int ivarvoln);
	void saveGeom(char *filename);
	void loadGeom(char *filename);
	void loadGeomBIN(FILE *inp);
	void loadGeomASC(FILE *inp);
	void loadGeomBZIP(FILE *inp);
	void (*defvar)
	(//Mod-function
		int	ivar,
		char	*name,
		int	&type,
		int	&rank,
		Domain	*dom
	);
	void	(*initvar)
	(//Mod-function
		int	ivar,
		char	*name,
		int	type,
		int	size,
		int	rank,
		double	*val,
		Domain	*dom
	);
	void	addDCell
	(	int	varbufsize,
	 	int	ivert,
		int	jvert,
	 	DNode	*&newvert,
		DCell	*&cell,
		int	iprev_neib_cell,
		int	ioldclone_neib_next,//local index of the next neighbor of the old clone
		DCell	*&oldclone,//neighbor cell's twin cell
		DCell	*&newclone
	);
	template <class Element>
	void insert(Element *&p, Element *&q);
	template <class Element>
	void insert_before(Element *&p, Element *&q);
	template <class Element>
	void alloc(int type, Element *&el);
	template <class Element>
	void	alloc
	(	int	ne, int type,
		Element	*&element_root
	);
	template <class Element>
	void	allocElementVars
	(	int	bufsize,
		Element	*element_root
	);
	int getnvar();
	void allocvar(int type, Variable *var);
	void allocgeovar(Variable *var);
	void setBoundaryVertexes();
	void setBoundaryCellFlags();
	void cleanInternalFaces();
	void createBoundaryNodeList(BNodeList *&root);
	void createBoundaryFaceList(BFaceList *&root);
	void updateBoundaryFaceList(BFaceList *&root);
	void createBoundaryCellList(DCellList *&root);
	void connectBoundaryCellsToFaces(BFaceList *&root);
	void setBC();
	int	Steps(int iter);
	void Advance(double new_time);
	void buildBrick(int *nx, int *bc, double *xmin, double *xmax);
	void brickConnect();
	void setFaces();
	void indexNeibNodes();
	void setVolumes();
	void getBrickBC(char *p, char *s, char *filename);
	void sort();
	void swapCells(int i, int j);
	void setDisplay();
	void	setDisplayVariable(int ivar);
	int	getDisplayVariable();
	void	setDisplayVariableComp(int ivar);
	int	getDisplayVariableComp();
	void	setElementColor(int element_type, float rgbcolor[]);
	void	getElementColor(int element_type, float rgbcolor[]);
	void	getGridLimits
	(
		double	xmin[],
		double	xmax[]
	);
	int	getNoVariables();
	int	getNoVariables(Elements vartype);
	int	getNoVariables(Elements vartype, int rank);
	int	getVarBufSize(Elements vartype);
	int	getNoBoundaryVertexes();
	int	getNoBoundaryFaces();
	int	getNoBoundaryCells();
	int	getLighting(int type);
	template <class Anything>
	void	swap(Anything *&a, Anything *&b);
	void	printVariables();
	//Point functions
	int	ihostcell(double *x, DCell *c);
	DCell	*hostcell(double *x, DCell *c);
	ElementStatus	hostcell
	(	double *z, //z[DIM] coordinates of a point 
		DCell *cell,//initial cell
		DCell	*&newcell
	);
	ElementStatus	hostcell//returns either bc or a face number
	(
		double *z, //z[DIM] coordinates of a point 
		DCell *cell,//initial cell
		DCell	*&newcell,
		int	*iface//cell face number
	);
	ElementStatus	hostcell//returns either bc or a face number
	(
		double *z, //z[DIM] coordinates of a point 
		DCell *cell,//initial cell
		DCell	*&newcell,
		double	*norm, //norm[DIM]: face-normal vector
		int	*iface//cell face number
	);
	DCell	*findcell(double *x);
	DCell	*findcell(double *x, DCell *c);
	ElementStatus	findcell
	(	double *x,//x[DIM] coordinates of a point 
		DCell *cell,//current cell
		DCell	*&newcell //new cell
	);
	ElementStatus	findcell//returns both cell and face number
	(	double *x,//x[DIM] coordinates of a point 
		DCell *cell,//current cell
		DCell	*&newcell, //new cell
		int	*iface //-1 if inside the cell, 0:Nf if outside the cell
	);
	ElementStatus	findcell
	(	double *x,//x[DIM] coordinates of a point 
		DCell *cell,//current cell
		DCell	*&newcell, //new cell
		double	*norm //norm[DIM]: boundary-face normal vector
	);
	void	initp();
	void	addp();
	int	putp
	(//Put point at a certain posisition and initialize
	 // a single scalar variable
		//int	iscl,//index to the scalar variable
		//double	scl,//initialization value for the scalar
		double	*x //coordinates
	);
	void	setsclp(int ip, int iscl, double scl);
	void	setvecp
	(//Initialize a single vector variable
		int	ip,  //index to the particle
		int	ivec,//index to the vector variable
		double	*vec//initialization value for the vector
	);
	void	injectParticlesAtFrame(int ninj);
	void	injectParticles(int ninj);
	int	injp
	(//Inject a point at a certain posisition 
	 // with a given velocity and initialize
	 // a single scalar variable
		int	nvar,
		int	*ivar,//index to the variables to be initialized
		double	*var,//variables initialization array
		double	*x //coordinates
	);
	int	delp(Point *&p);
	int	getNoPoints();
	int	getMaxNoPoints();
	void	getPoints
	(
		int &m,
		int	&n,
		struct Point *&origin0,
		struct Point *&first0,
		struct Point *&last0,
		double	*&X
	);
	void	deleteDNode(DNode *&dnode);
	template	<class List>
	void	deleteRing(List *&element);
	template	<class List>
	void	deleteList(List *&head);
	template <class Element>
	void	deleteElement(Element *&element, Element *&root);
	template <class Element>
	void	deleteGrid(Element *&root);
	void	relax
	(
		int	ixold,
	 	int	ivoln,
	 	int	ivolc,
		double	relaxation_factor
	);
	void	checkGrid
	(	DNode	*node_root,
		DCell *cell_root
	);
	//DIFFERENTIAL ALGEBRA:
	void	VijUj   (int iVij , int iUj , int iVijUj  );
	void	dVdXi   (int ivoln, int	iV	, int idVdXi  );
	void	dVidXi  (int ivoln, int	iVi , int idVidXi );
	void	dVidXj  (int ivoln, int	iVi , int idVidXj );
	void	dVijdXj (int ivoln, int	iVij, int	idVijdXj);
	void	dViUjdXj(int ivoln, int	iVi , int iUj, int	idViUjdXj);
	void	dVidXjXj(int ivoln,	int	iVi ,	int	idVidXjXj);
	//PARSING
	int	parseStringNodes(int nnodes, char *string);
	int	parseStringCells(int nnodes, int ncells, char *string);
	//OUTPUT
	void	setOutputDefault();
	void	setOutputEnsight();
	void	setOutputTecplot();
	void	setOutputUser   ();
	void	outputData(char *filename);
	void	outputDefault(char *filename);
	void	outputEnsight(char *filename);
	void	outputTecplot(char *filename);
	void	outputUser(char *filename);
	void	saveDispVarDefault(int ivar, char *filename);
	void	saveGeomDefault(char *filename);
	void	saveGeomEnsight(char *filename);
	void	saveGeomTecplot(char *filename);
	void	saveGeomUser(char *filename);
	
	#ifdef WITH_MPI
	void	sendDispVar2Proc0(int ivar);
	#endif
};

extern int	ndomains;
extern Domain	*domain_root;

