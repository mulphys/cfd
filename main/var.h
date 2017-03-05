class	Variable
{
	struct
	{
		unsigned int	output	:1;
		unsigned int	initialize	:1;
	}	status;
	public:
	int
		type,//type of the element the variable is defined for
		     //node,edge,face,cell; 
		size,//size=Ne[type] for grid var
			   //size=determined from config file for local var's
		level,//current nesting level where the variable was defined
		dimension,//
		loc;//location of the variable in the element var-buffer
	char
		name[MAXLINLEN],
		*inpfile,
		*outfile;
	double
		*initval,
		*min,*max,//min[0:DIM^rank]=variable limi values
		*val;//val[DIM^irank*(0:Ne[type])+0:DIM^irank]=tensor components
	void	setinit(int dimension);
	void	delinit();
	void	init(int nelements);
	void	reset();
	void	computeLimits();
	void	getLimits(int icomp, double &min, double &max);
};

