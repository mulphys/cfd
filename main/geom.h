#define PI	3.14159265358979323846
#define TWOPI	6.283185307179586232

#ifndef double
#define double	double
#endif

#ifndef DIM
#define DIM	3
#endif
#ifndef	Nv
/* number of vertices of the cell */
#define	Nv	4
#define	Nv1	3
#endif
#ifndef	Nf
/* number of faces of the cell */
#define	Nf	4
#define	Nf1	3
/* number of face-vertices of the face */
#define	Nfv	3
#define	Nfv1	2
#endif

enum	GridType
{
	structured_grid=0,
	block_structured_grid,
	unstructured_tetrahedral_grid,
	unstructured_tetrahexa_grid,
	unstructured_mixed_grid,
	max_grid_types
};
enum	DomainType
{
	generic=0,
	box,
	cylinder,
	user,
	dynamic,
	max_domain_types
};
enum  ElementStatus
{
   maxElementStatus=-14, //if >= maxElementStatus -> parallel domain
   dead,       //-13
   dove,       //-12 <-- this is the new line
   internal,   //-11
   insidetool, //-10
   presoutlet, //-9
   presinlet,  //-8
   outlet,     //-7
   inlet,      //-6
   boundary,   //-5
   bounce,     //-4
   slip,       //-3
   noslip,     //-2
   composite  //-1
};
enum	Elements
{
	points=0,
	nodes,
	edges,
	faces,
	cells,
	boundary_nodes,
	boundary_faces,
	boundary_cells,
	maxelements
};
struct	Loop
{
	struct Loop	*next,*prev;
};
struct	Edge
{
	int	node[2];
	struct Edge *next,*prev;
};
struct	Face
{	int
		vert[Nfv],
		cell[2],//global indexes to neighbor cells
		        // for the boundary face: 
		        // cell[0] points to the internal cell
		        // cell[1] specifies the boundary conditions
		cf;//local indexes to this face from neighbor cells
	double
		area,
		norm[DIM];//normal vector directed toward cell[1]
};
/****************
*             this vertex--- >  o      
*                              / \     
*                             /   \    
*             this cell--->  /     \   
*                           o-------o  
*                            \     /   
*                             \   /    
*             neighbor cell    \ /     
*           opposite vertex-->  o      
*/
struct NodeState
{
	unsigned int	insidetool:	1;
	unsigned int	toolsurface:	1;
	unsigned int	boundary:	1;
	unsigned int	relax:	1;
	unsigned int	fixed:	1;	
};
struct	DNode
{//Dynamic node
	int	index;//used in file IO
	struct NodeState	state;
	int	type; //<0 => ElementStatus, >=0 => index to boundary vertex array
	struct Bond	*bond;//Neighbor-pointer list
	double
		x[DIM],//Coordinates
		*var;  //Variable storage
	struct DNode	
		*next,*prev;
	DNode();
	DNode(int mvar);
	DNode(int mvar, double x[]);
	~DNode();
};
struct	Bond
{
	DNode	*node;
 	Bond	*next;
	double	x[DIM];//initial bond positions (used in lattice model)
};
struct	DEdge
{
	DNode	*node[2];
	DEdge *next,*prev;
};
struct	DCell
{//Dynamic cell: can be created or destroyed
	struct
	{
		unsigned int	boundary:	1;
		unsigned int	relax:	1;
	}	state;
	unsigned int	index; //used in file IO
	ElementStatus	facetype[Nf]; //<0: Boundary conditions on the faces,
	                      //>=0: index of the parallel domain
	DNode	*vert[Nv];
	void	*neib[Nf];//neighbors or boundary elements
	DCell	*next,*prev;//Used in dynamic cell creation/destruction
	double	*var;
	DCell(void);
	DCell(int varbufsize);
	DCell(DCell *c);
	DCell(DCell *c, int varbufsize);
	~DCell();
};
struct	DNodeList
{
	DNode	*node;
	DNodeList	*next;
};
struct	BNodeList
{//Boundary node list
	DNode	*node;
//	double	*var;//nodal variables
	BNodeList	*next;
};
struct	DCellList //Cell pointer list
{
//	int	status;//For boundary cells:
	          //-1: no refinement, 
	          //>0: type/6=number of the first vertex of the edge 
	          // to be refined, type%6=second vertex of the edge to be refined
	DCell	*cell;
	DCellList	*next;
};
struct	BFaceList //Boundary cell pointer list
{
	int	iface;//boundary face index from the cell: iface=0:Nf-1
	double
		x[DIM],//face center coordinates
		area,norm[DIM];//area and its normal vector
	ElementStatus	type;
	DCell	*cell;//ponter to the cell
	            //there can be up to 4 DCells pointing
	            //to the same DCell
	BFaceList	*next;
	double	*var;
	BFaceList();
	BFaceList(int mvar);
	~BFaceList();
};
struct	BoundaryVertex
{
	unsigned int	i;/* index of this node for a static domain */
	int	b; //>=maxElementStatus => number of the parallel domain
	double	norm[DIM];//boundary normal vector
};
struct	Point
{
	void	*cell;//host-cell 
	Point	*next,*prev;
};
struct	Sphere
{
	double	radius;
};
extern	char	*elementype[];
extern	char	*spantype[];
extern double	nullvector[];
extern unsigned int	mask[Nv]; //binary mask used to retrieve face indexes from C[].f
//Computes the edge-vectors of a tetrahedron
#define EDGES(e,verts)	{double	*x=verts[Nv1]->x;      \
			for (int iv=0; iv<Nv; iv++)                    \
			{ double	*y=verts[iv]->x;                       \
				for (int k=0; k<DIM; k++)e[iv][k]=y[k]-x[k]; \
				x=y;                                         \
			}}
//Computes the area-vector of a tetrahedron's face (actually twice the area)
#define	AREAVEC(iv,e,a) { double         \
					*e0=e[(iv+1)%Nv],             \
					*e1=e[(iv+2)%Nv],             \
					*e2=e[(iv+3)%Nv];             \
				VECP(a,e1,e2);                 \
				if(SCLP(e0,a)<0.0)             \
					for(int i=0;i<DIM;i++)       \
						a[i]*=-0.5;                \
				else                           \
					for(int i=0;i<DIM;i++)       \
						a[i]*= 0.5;                \
                      }
#define	AREA2(iv,e,a) { int            \
					iv0=(iv+1)%Nv,               \
					iv1=(iv+2)%Nv,               \
					iv2=(iv+3)%Nv;               \
				VECP(a,e[iv1],e[iv2]);         \
				if(SCLP(e[iv0],a)<0.0)         \
					for(int i=0;i<DIM;i++)       \
						a[i]=-a[i];                \
                      }
#define	FACECENTER(ivert,verts,xf) for(int i=0;i<DIM;i++) \
				{	double	xcenter=0.0;                              \
					for(int jvert=1;jvert<Nv;jvert++)               \
						xcenter+=verts[(ivert+jvert)%Nv]->x[i];       \
					xf[i]=xcenter/(double)Nv1;                        \
				}
