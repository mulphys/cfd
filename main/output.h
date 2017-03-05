namespace	Output
{
	extern	char	ioconfigfile[];
	struct OutputStatus
	{
		unsigned int	active: 1;
	};
	extern enum	OutputTypes	outype;
	extern char	*outypename[];
	//output.cc
	void initOutput();
	void setOutput();
	void readGeom(char *filename);
	void writeGeom(char *filename);
	void outputBoundary();
///		void setOutputDefault();
///		void setOutputUser();
///		void setOutputTecplot();
///		void setOutputEnsight();
///		void outputEnsight(unsigned int);
///		void outputTecplot(unsigned int);
	void toggle();
	void xwd(char *windowname);
};
