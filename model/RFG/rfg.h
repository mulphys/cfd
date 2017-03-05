#ifndef	LINELENGTH
#define LINELENGTH	510
#endif
#ifndef	INT
#define	INT	int
#endif
#ifndef	double
#define	double	double
#endif
#ifndef	DIM
#define	DIM	3
#endif

extern	int
	ne; /* Number of terms in the series Eq.(15),celik.bib:\cite{LiAhetalJAS94} */

extern	double
	fe,
	*Omega, /* Eq.15,celik.bib:\cite{LiAhetalJAS94} */
	*U1,*U2, /* velocity vectors (Eqs.15-17),celik.bib:\cite{LiAhetalJAS94} */
	*K; /* wave vectors (Eqs.15-17),celik.bib:\cite{LiAhetalJAS94} */
