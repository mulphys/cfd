#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "rfg.h"
#include "vecalg.h"

#ifndef SMALL
#define SMALL	1.e-30
#endif

int	ne = 1; /* Number of terms in the series Eq.(15),celik.bib:\cite{LiAhetalJAS94} */

double
	fe = 1.41421,
	*Omega, /* Eq.15,celik.bib:\cite{LiAhetalJAS94} */
	*U1,*U2, /* velocity vectors (Eqs.15-17),celik.bib:\cite{LiAhetalJAS94} */
	*K; /* wave vectors (Eqs.15-17),celik.bib:\cite{LiAhetalJAS94} */

void    Allocate(double **A, int n)
{
	if ((*A = (double *) malloc(sizeof(**A)*n)) == NULL)
	{
		fprintf(stderr,"CAN'T ALLOCATE MEMORY\n");
		exit(1);
	}
}
void	genspec_
(
	int	*Ne
)
/*
 *	Generate spectral expansion coefficients 
 */
{
	extern	double	gauss_(); /* get a Gaussian random variable */
	extern	void	seed_(), gaussn_(double *, double, int); /* get an array of Gaussian random numbers */
	int	i,ie;
	double pi=4.*atan(1.0);
	ne=*Ne;
	if (ne<=0) return;
	fe=sqrt(2./(double)ne);
	Allocate(&Omega,ne);
	Allocate(&K ,ne*DIM);
	Allocate(&U1,ne*DIM);
	Allocate(&U2,ne*DIM);
	seed_(); /* initialize the random generator */
//	seed0(666); //DEBUG
	gaussn_(K,.5,ne*DIM);
	for (ie=0; ie<ne; ie++)
	{	int	j=ie*DIM;
		double	a,V1[DIM],V2[DIM], /* random vectors xi and zeta (Eq.16) */
			*u1=U1+j,
			*u2=U2+j,
			*k=K+j;
		*k*=2.0*pi;
		Omega[ie] = 2.0*pi*gauss_(); 
		gaussn_(V1, 1., DIM);
		gaussn_(V2, 1., DIM);
		VECP(u1,V1,k); /* Eq.16\cite{LiAhetalJAS94} */
		VECP(u2,V2,k);
//		for (i=0; i<DIM; i++)
//			k[i]/=d[i]>SMALL?turb_time*d[i]:TL[i];
	}
}
void	delspec_()
{
	free(Omega);
	free(K);
	free(U1);
	free(U2);
}
void	genvel_
(
//	INPUT:
	double	*t,  // time
	double	*x,  // coordinates
	double	*TT, // Turbulent Time: scalar
	double	*TL, // Turbulent Length: vector: TL[0:2]
	double	*UU, // velocity correlations: UU,UV,VV,UW,VW,WW
// OUTPUT:
	double	*v   // velocities
)
{
#ifdef FORTRAN
	extern	void	diag
	(
	  double *A, /* velocity correlations */  
	  double *D  /* diagonal vector after diagonalization of A */
	);
#endif
	int	i,ie;
	double	a,c,s,
		d[DIM],
		turb_time=*TT;

	if (ne<=0)  return;
#ifdef FORTRAN
	diag(UU,d);
#endif
//	ftrans(wn); /* principal axes of velocity correlation tensor 
//	               are alligned with turbulent length-scales */
	for (i=0; i<DIM; i++)
	{	v[i]=0.0;
		d[i]=sqrt(fabs(d[i]));
	}
	for (ie=0; ie<ne; ie++)
	{	int n=ie*DIM;
		double	k[DIM],
		      *u1=U1+n,*u2=U2+n;
		for (i=0; i<DIM; i++)
			k[i]=K[n+i]/TL[i];
//			k[i]=K[n+i]/(d[i]>SMALL?turb_time*d[i]:TL[i]);
//		for (i=0; i<DIM; i++)k[i]=K[n+i]*wn[i];//K-anisotropy
		a=SCLP(k,x)+Omega[ie]**t/turb_time;
		c=cos(a); s=sin(a);
		for (i=0; i<DIM; i++) 
			v[i]+=u1[i]*c+u2[i]*s;
	}
	for (i=0; i<DIM; i++)v[i]*=fe*d[i];//V-anisotropy
#ifdef FORTRAN
	btrans(v);
#endif
}

