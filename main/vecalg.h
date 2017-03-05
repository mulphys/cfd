/*
	Vector algebra
*/

//	#ifndef	double
//	#define	double	double
//	#endif
//	#ifndef	DIM
//	#ifdef	DM
//	#define	DIM	DM
//	#else
//	#define	DIM	3
//	#endif
//	#endif

//extern int	E[DIM][DIM][DIM]; /* assymmetric tensor used in vector product */

#ifndef ZERO
#define	ZERO(n,x)	for (int i=0; i<(n); i++)x[i]=0.0
#endif

#define ADD(n,x,y)	for (int i=0; i<(n); i++) (x)[i] += (y)[i]
#define SUB(n,x,y)	for (int i=0; i<(n); i++) (x)[i] -= (y)[i]
#define ADDC(n,x,y,d)	for (int i=0; i<(n); i++) (x)[i] += (y)[i]*(d)
#define SUBC(n,x,y,d)	for (int i=0; i<(n); i++) (x)[i] -= (y)[i]*(d)
#define MUL(n,x,y)	for (int i=0; i<(n); i++) (x)[i] *= (y)[i]
#define MULC(n,x,c)	for (int i=0; i<(n); i++) (x)[i] *= (c)
#define DIV(n,x,y)	for (int i=0; i<(n); i++) (x)[i] /= (y)[i]
#define MULVEC(x,s)	for (int i=0; i<DIM; i++) (x)[i] *= (s)
#define ADDVEC(x,y)	for (int i=0; i<DIM; i++) (x)[i] += (y)[i]
#define COPYVEC(x,y)	for (int i=0; i<DIM; i++) (x)[i] = (y)[i]
#define ZEROVEC(x)	for (int i=0; i<DIM; i++) (x)[i] = 0.
#define ZERO3(x)	(x)[0] = (x)[1] = (x)[2] = 0.0
#define COPY3(x,y)	(x)[0]=(y)[0];(x)[1]=(y)[1];(x)[2]=(y)[2];
#define SCLP(A,B)	(A[0]*B[0]+A[1]*B[1]+A[2]*B[2])
#define VECP(A,B,C)	A[0]=B[1]*C[2]-B[2]*C[1];\
							A[1]=B[2]*C[0]-B[0]*C[2];\
							A[2]=B[0]*C[1]-B[1]*C[0];
#define SCLP1(A,B)	*(A+0)*B[0]+*(A+1)*B[1]+*(A+2)*B[2]
#define LENGTH(A)	sqrt(SCLP(A,A))

#define DET(a)  - a[0][2]*a[1][1]*a[2][0] + a[0][1]*a[1][2]*a[2][0] \
                + a[0][2]*a[1][0]*a[2][1] - a[0][0]*a[1][2]*a[2][1] \
                - a[0][1]*a[1][0]*a[2][2] + a[0][0]*a[1][1]*a[2][2]

#define	INV( d, a, b, x)	x[0]=                              \
         -(  b[2]*a[0][2]*a[1][1] - b[2]*a[0][1]*a[1][2]     \
           - b[1]*a[0][2]*a[2][1] + b[0]*a[1][2]*a[2][1]     \
           + b[1]*a[0][1]*a[2][2] - b[0]*a[1][1]*a[2][2])/d; \
   x[1]= -(- b[2]*a[0][2]*a[1][0] + b[2]*a[0][0]*a[1][2]     \
           + b[1]*a[0][2]*a[2][0] - b[0]*a[1][2]*a[2][0]     \
           - b[1]*a[0][0]*a[2][2] + b[0]*a[1][0]*a[2][2])/d; \
   x[2]=  (- b[2]*a[0][1]*a[1][0] + b[2]*a[0][0]*a[1][1]     \
           + b[1]*a[0][1]*a[2][0] - b[0]*a[1][1]*a[2][0]     \
           - b[1]*a[0][0]*a[2][1] + b[0]*a[1][0]*a[2][1])/d

//	#define DET(d) d[0][0]*d[1][1]*d[2][2] \
//								+d[0][1]*d[1][2]*d[2][0] \
//								+d[1][0]*d[2][1]*d[0][2] \
//								-d[2][0]*d[1][1]*d[0][2] \
//								-d[0][0]*d[2][1]*d[1][2] \
//								-d[1][0]*d[0][1]*d[2][2] 
//	#define INV(dd,d,b,g) g[0]= \
//				 (b[0]*d[1][1]*d[2][2]+d[0][1]*d[1][2]*b[2]+b[1]*d[2][1]*d[0][2] \
//				 -b[2]*d[1][1]*d[0][2]-b[0]*d[2][1]*d[1][2]-b[1]*d[0][1]*d[2][2])/dd; \
//		g[1]=(d[0][0]*b[1]*d[2][2]+b[0]*d[1][2]*d[2][0]+d[1][0]*b[2]*d[0][2] \
//				 -d[2][0]*b[1]*d[0][2]-d[0][0]*b[2]*d[1][2]-d[1][0]*b[0]*d[2][2])/dd; \
//		g[2]=(d[0][0]*d[1][1]*b[2]+d[0][1]*b[1]*d[2][0]+d[1][0]*d[2][1]*b[0] \
//				 -d[2][0]*d[1][1]*b[0]-d[0][0]*d[2][1]*b[1]-d[1][0]*d[0][1]*b[2])/dd

//	void	inivec();
//	void	vecp
//	/*
//		Vector product A = [B,C]
//	*/
//	(
//		double *A,
//		double *B,
//		double *C
//	);
//	double	average
//	(
//		double	*A
//	);
//	double	sclp
//	/* 	Scalar product: a = (B,C)
//	*/
//	(
//		double *A,
//		double *B
//	);
