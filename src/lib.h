
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define pi 3.14159265
//#define RAND_MAX 0xffffffff

#define TINY 1.0e-20;
void ludcmp(a,n,indx,d);
#undef TINY
void lubksb(a,n,indx,b);
double matrix_logdet(X, n);
/* Y=inv(X), return d=log(det(X)) */ 
double matrix_inverse(X, Y, n);
/* Y=inv(X), return d=log(det(X)) */
int matrix_inverse_diag(X, Y, diag, n);
double matrix_trace(A,p);
int matrix_sum(A,B,C,n,p);
/* Matrix: A: n by p; B: p by m;  C: n by m */
int matrix_multiply(A,B,C,n,p,m);
int matrix_vector_prod(A,b,d,n,p);
double vector_matrix_vector(a,X,b,m,n);
void copy_vector(a,b,p)         ;                                                                                                                     
void copy_matrix(a,b,n,p);
int choldc(a, n, D);
/* calculate log(Gamma(s))  */
double loggamma(xx);
/* calculate log(k!) */
double logpum(k);

/* generate the random variable form Gamma(a,b) */
double Rgamma(a,b);

/* Generate a random variable from Beta(1,k), where
  the first parameter is 1, the second parameter is b */
double Rbeta(b);

/* Generate deviates from Dirichlet(a1,a2,\ldots,a_k) */
int RDirichlet(w,a,k);

double gasdev();

double Rgasdev(mean,variance);

int RNORM(x,mu,Sigma,p);

int Rwishart(B,df,Sigma,p);

/* calculated the log-density of  z~gamma(a,b) */
double dloggamma(x,a,b);

double dloggauss(z,mean,variance);

double DLOGGAUSS(z,mean,variance,p);

double Dlogwishart(D,df,Sigma,p);

int uniform_direction(d, n);

int dmaxclass(z,n);
                                                                                                                        
int imaxclass(z,n);

int binary_trans(k,l,d);
                                                                                                                                         
double logsum(a,b);

double maxvector(x,n);

double minvector(x,n);
                                                                                                                                        
double sample_variance(x,n);

/* Return the value ln[Gamma(xx)] for xx>0 */
double gammln(xx);

#define ITMAX 100
#define EPS 3.0e-7
void gser(gamser,a,x,gln);

#undef ITMAX
#undef EPS

#define ITMAX 100
#define EPS 3.0e-7
void gcf(gammcf,a,x,gln);
#undef ITMAX
#undef EPS


double gammp(a,x);

/* Return the CDF of the standard normal distribution */
double Gaussp(x);

/* return Gamma'(z)/Gamma(z)   */
/* Refer to "mathematics handbook pp.287" */
double diGamma(z);

double correlation(z1,z2,p);

int permut_sample(sam,n);
int random_order(x,n);

/* Generate a subset sample of size M from the set 1:N */
int subset_sample(x,M,N);

void indexx(n,arrin,indx);

float* vector();
float** matrix();
float** convert_matrix();
double* dvector();
double** dmatrix();
int* ivector();
int** imatrix();
float** submatrix();
void free_vector();
void free_dvector();
void free_ivector();
void free_matrix();
void free_dmatrix();
void free_imatrix();
void free_submatrix();
void free_convert_matrix();
void nrerror();
