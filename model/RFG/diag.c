/*
	Diagonalization of a real symmetric matrix 

		A=\{a_{ij}\}

	Source: Numerical Recipes in C

	Equations: 

	1. eq.(11.0.13)

	A = Transpose[Z].A'.Z
=> 
a_{ij}=z_{ii'}a'_{i'j'}z_{j'j}
% Checking the inverse transform:
z_{ki}a_{ij}=z_{ki}z_{ii'}a'_{i'j'}z_{j'j}
            =\delta_{ki'}a'_{i'j'}z_{j'j}
            =a'_{kj'}z_{j'j}
z_{ki}a_{ij}z_{jl}=a'_{kj'}z_{j'j}z_{jl}
z_{ki}a_{ij}z_{jl}=a'_{kj'}\delta_{j'l}
                  =a'_{kl}
=>a'_{kl}=z_{ki}a_{ij}z_{jl}
=> the inverse is consistent

A regular tensor relation:

a_{ij}=z_{ii'}z_{jj'}a'_{i'j'}
% Checking the inverse transform:
z_{ki}a_{ij}=z_{ki}z_{ii'}z_{jj'}a'_{i'j'}
z_{ki}a_{ij}=\delta_{ki'}z_{jj'}a'_{i'j'}
            =z_{jj'}a'_{kj'}
z_{lj}z_{ki}a_{ij}=z_{lj}z_{jj'}a'_{kj'}
                  =\delta_{lj'}a'_{kj'}=a'_{kl}
=> a'_{kl}=z_{ki}z_{lj}a_{ij}
=> the inverse is consistent

*/
#include	<stdio.h>
#include	<math.h>

//#include "mainrfg.h"
#include "rfg.h"

#define	M	3
#define	N	M
#define	I(i,j)	(i)*N+(j)

double	Z[M*N];

void	diag
(
	double	*UU, /* IN: correlations uu,uv,vv,uw,vw,ww */
	double	*D  /* OUT: diagonal elements of A after diagonalization */
)
{
	INT	i,j,k,
				m=M,n=N,ierr;
	double	A[M*N],W[N],F[N];
	k=0;
	for (i=0; i< n; i++)
	for (j=0; j<=i; j++)
	{	A[I(j,i)]=UU[k++];
		A[I(i,j)]=A[I(j,i)];
	}
	tred2_(&m,&n,A,W,F,Z);
	tql2_(&m,&n,W,F,Z,&ierr);
	if (ierr!=0)
	{
		fprintf(stderr,"\nWARNING: Diagonalization failed for the correlation matrix:\n");
		for (i=0; i<n; i++)
		{
			for (j=0; j<n; j++)
				fprintf(stderr,"\t%g",A[I(i,j)]);
			fprintf(stderr,"\n");fflush(stderr);
		}
		return;
	}
	for (i=0; i<n; i++)
	{	int	ii,jj;
		double	s=0.;
		j=i;
		for (ii=0; ii<n; ii++)
		for (jj=0; jj<n; jj++)
			s+=Z[I(i,ii)]*Z[I(j,jj)]*A[I(ii,jj)];
		D[i]=s;
	}
}
void	ftrans
(
	double	*v
)
{	int	i,j;
	double	u[N];
	for (i=0; i<N; i++)
	{	double	vi=0.;
		for (j=0; j<N; j++)
			vi+=Z[I(i,j)]*v[j];
		u[i]=vi;
	}
	for (i=0; i<N; i++)
		v[i]=u[i];
}
void	btrans
(
	double	*v
)
{	int	i,j;
	double	u[N];
	for (i=0; i<N; i++)
	{	double	vi=0.;
		for (j=0; j<N; j++)
			vi+=Z[I(j,i)]*v[j];
		u[i]=vi;
	}
	for (i=0; i<N; i++)
		v[i]=u[i];
}
