// #include <stdio.h>
// #include <stdlib.h>
// #include <math.h>

#include <R.h>
#include "RBF.h"

#define ZERO 1e-10
#define TOL_INVERSE ZERO
#define min(a,b) ((a)<(b)?(a):(b))

double psi_huber_w( double r, double k) {
double sal, aux1, aux2;
aux1 = fabs(r);
aux2 = k/aux1;
// sal = min(aux2, 1.0);
if( aux2 < 1 ) { sal = aux2; }
else { sal = 1;}
return( sal );
}


double psi_tukey_w( double r, double k) {
double sal=0.0, aux1;
aux1 = fabs(r/k);
if( aux1 < 1.0 ) {
	sal = (1+aux1) * (1+aux1) * (1-aux1) * (1-aux1);
};
return( sal );
}

// MAL
// double psi_tukey_w( double r, double k) {
// double sal, aux1, aux2;
// aux1 = fabs(r/k);
// aux2 = (1+aux1) * (1+aux1) * (1-aux1) * (1-aux1);
// if( aux2 < 1 ) { sal = aux2; }
// else { sal = 0;}
// return( sal );
// }




//31
//Nucleo de Epanechnikov

double kepan( double b ) {
double a, aux;
a = 0.75 * (1 - (b*b) );
aux = fabs(b);
if( aux < 1 ) { return(a); }
else { return(0); }
}
//41
//Suma de coordenadas de un vector


double suma_vec(double *x, int n){
register int i;
double sal=0;
for(i=0;i<n;i++){
	sal = sal + x[i];
}
return(sal);
}

double norma_dos(double *x, int n){
register int i;
double sal=0;
for(i=0;i<n;i++){
	sal = sal + (x[i])*(x[i]);
}
return(sqrt(sal));
}

double l2dist(double *x, double *y, int n){
register int i;
double sal=0;
for(i=0;i<n;i++){
	sal = sal + (x[i] - y[i]) * (x[i] - y[i]);
}
return(sqrt(sal));
}




void kernel_huber_pos(double *punto, double *x, int *nrow, double *y, double *muhat_initial,
 double *ventanas, double *eps, double *sigmahat, double *prob,
 double *kh, int *maxit, double *salida) {

double kepan(double);
double psi_huber_w(double, double);
double suma_vec(double*, int);
register int i;
int n = *nrow, it=0;
double corte, muhat, muold, khh;
double *pesos, *aux1, *aux2, *res;

pesos = (double *) malloc( n * sizeof(double));
aux1 = (double *) malloc( n * sizeof(double));
aux2 = (double *) malloc( n * sizeof(double));
res = (double *) malloc( n * sizeof(double));
corte = 10 * (*eps);
khh = *kh;
muhat = *muhat_initial;
for(i=0;i<n;i++) {
     pesos[i] = kepan( (x[i] - *punto) / *ventanas )  / prob[i];
}

while( (corte > (*eps)) && (it< (*maxit)) ){
	for(i=0;i<n;i++){
		res[i] = ( y[i] - muhat) / *sigmahat;
		aux1[i] = pesos[i] * psi_huber_w(res[i],khh) * y[i];
		aux2[i] = pesos[i] * psi_huber_w(res[i],khh);
	};
	muold = muhat;
	muhat = suma_vec(aux1,n) / suma_vec(aux2,n);
	corte = fabs(muold-muhat) / ( fabs(muold) + *eps );
	it = it + 1;
};
*salida = muhat;
free(pesos); free(aux1); free(aux2); free(res);
return;
}


void kernel_tukey_pos(double *punto, double *x, int *nrow, double *y, double *muhat_initial,
 double *ventanas, double *eps, double *sigmahat, double *prob,
 double *kt, int *maxit, double *salida) {

double kepan(double);
double psi_tukey_w(double, double);
double suma_vec(double*, int);
register int i;
int n = *nrow, it=0;
double corte, muhat, muold, ktt;
double *pesos, *aux1, *aux2, *res;

pesos = (double *) malloc( n * sizeof(double));
aux1 = (double *) malloc( n * sizeof(double));
aux2 = (double *) malloc( n * sizeof(double));
res = (double *) malloc( n * sizeof(double));
corte = 10 * (*eps);
ktt = *kt;
muhat = *muhat_initial;
for(i=0;i<n;i++) {
     pesos[i] = kepan( (x[i] - *punto) / *ventanas )  / prob[i];
}

while( (corte > (*eps)) && (it< (*maxit)) ){
	for(i=0;i<n;i++){
		res[i] = ( y[i] - muhat) / *sigmahat;
		aux1[i] = pesos[i] * psi_tukey_w(res[i],ktt) * y[i];
		aux2[i] = pesos[i] * psi_tukey_w(res[i],ktt);
	};
	muold = muhat;
	muhat = suma_vec(aux1,n) / suma_vec(aux2,n);
	corte = fabs(muold-muhat) / (fabs(muold) + *eps );
	it = it + 1;
};
*salida = muhat;
free(pesos); free(aux1); free(aux2); free(res);
return;
}




void kernel_huber_lin(double *punto, double *x, int *nrow, double *y,
 double *z, int *degree, double *beta_initial,
 double *ventanas, double *eps, double *sigmahat, double *prob,
 double *kh, int *maxit, double *salida) {

double kepan(double);
double psi_huber_w(double, double);
double suma_vec(double*, int);
void reset_mat(double**, int, int);
void reset_vec(double*, int);
void sum_mat(double **a, double **b, double **c, int n, int m);
void vec_vecprime(double **a, double *v1, double *v2, int n);
void scalar_vec(double *a, double b, double *c, int n);
void sum_vec(double *a, double *b, double *c, int n);
double vecprime_vec(double *a, double *b, int n);
void mat_vec(double **a, double *b, double *c, int n, int m);
double norma_dos(double *x, int n);
double l2dist(double *x, double *y, int n);
int lu(double **a,int *P, double *x);
register int i, j;
int n = *nrow, it=0, q = *degree + 1;
double corte, khh; // muhat, muold, khh;
double *pesos, *aux1, *aux2, *res, **zz, *beta, *beta_old;
double **tmp1, **tmp2, *tmp3, *tmp4, **tmp5;

zz = (double **) malloc( n * sizeof(double *));
for(i=0;i<n;i++) zz[i] = (double *) malloc ( q * sizeof(double) );

pesos = (double *) malloc( n * sizeof(double));
aux1 = (double *) malloc( n * sizeof(double));
aux2 = (double *) malloc( n * sizeof(double));
res = (double *) malloc( n * sizeof(double));
beta = (double *) malloc( q * sizeof(double));
beta_old = (double *) malloc( q * sizeof(double));

tmp1 = (double **) malloc( q * sizeof(double *));
tmp2 = (double **) malloc( q * sizeof(double *));
tmp5 = (double **) malloc( q * sizeof(double *));
tmp3 = (double *) malloc( q * sizeof(double));
tmp4 = (double *) malloc( q * sizeof(double));

for(j=0;j<q;j++) {
	tmp1[j] = (double *) malloc( q * sizeof(double) );
	tmp2[j] = (double *) malloc( (q + 1) * sizeof(double) );
	tmp5[j] = (double *) malloc( q * sizeof(double) );
}

corte = 10 * (*eps);
khh = *kh;

for(i=0;i<n; i++) {
     pesos[i] = kepan( (x[i] - *punto) / *ventanas )  / prob[i];
     for(j=0; j<q; j++) zz[i][j] = *(z + j*n + i);
};
reset_mat(tmp1, q, q);
reset_vec(tmp3, q);
reset_mat(tmp5, q, q);
reset_vec(beta, q);
reset_vec(beta_old, q);
for(j=0; j<q; j++) beta[j] = beta_initial[j];
while( (corte > (*eps)) && (it< (*maxit)) ){
	reset_mat(tmp2, q, q+1);
	reset_vec(tmp4, q);
	for(j=0; j<q; j++) beta_old[j] = beta[j];
	for(i=0;i<n;i++){
		res[i] = ( y[i] - vecprime_vec(zz[i], beta, q) ) / *sigmahat;
		aux1[i] = pesos[i] * psi_huber_w(res[i],khh) * y[i];
		aux2[i] = pesos[i] * psi_huber_w(res[i],khh);
		scalar_vec(zz[i], aux2[i], tmp3, q);
		vec_vecprime(tmp1, tmp3, zz[i], q);
		sum_mat(tmp1, tmp2, tmp2, q, q);
		scalar_vec(zz[i], aux1[i], tmp3, q);
		sum_vec(tmp3, tmp4, tmp4, q);
	};
	for(j=0;j<q;j++) tmp2[j][q] = tmp4[j];
	if( lu(tmp2, &q, beta) == 1) {
		it = *maxit; // System became singular
		for(j=0; j<q; j++) beta[j] = NA_REAL;
	};
	// if( lu(tmp2, &q, beta) == 1) it = *maxit; // System became singular
	corte = l2dist(beta, beta_old, q) / ( norma_dos(beta_old, q) + *eps);
	it = it + 1;
};
for(j=0; j<q; j++) salida[j] = beta[j];
free(pesos); free(aux1); free(aux2); free(res);
for(i=0; i<n; i++) free(zz[i]); free(zz);
for(j=0;j<q;j++) {
	free( tmp1[j] ); free( tmp2[j] ); free( tmp5[j] );
}
free(tmp1); free(tmp2); free(tmp3); free(tmp4); free(tmp5);
free(beta); free(beta_old);
return;
}


void kernel_tukey_lin(double *punto, double *x, int *nrow, double *y,
 double *z, int *degree, double *beta_initial,
 double *ventanas, double *eps, double *sigmahat, double *prob,
 double *kt, int *maxit, double *salida) {

double kepan(double);
double psi_tukey_w(double, double);
double suma_vec(double*, int);
void reset_mat(double**, int, int);
void reset_vec(double*, int);
void sum_mat(double **a, double **b, double **c, int n, int m);
void vec_vecprime(double **a, double *v1, double *v2, int n);
void scalar_vec(double *a, double b, double *c, int n);
void sum_vec(double *a, double *b, double *c, int n);
double vecprime_vec(double *a, double *b, int n);
void mat_vec(double **a, double *b, double *c, int n, int m);
double norma_dos(double *x, int n);
double l2dist(double *x, double *y, int n);
int lu(double **a,int *P, double *x);

register int i, j;
int n = *nrow, it=0, q = *degree + 1;
double corte, ktt; // muhat, muold, ktt;
double *pesos, *aux1, *aux2, *res, **zz, *beta, *beta_old;
double **tmp1, **tmp2, *tmp3, *tmp4, **tmp5;

zz = (double **) malloc( n * sizeof(double *));
for(i=0;i<n;i++) zz[i] = (double *) malloc ( q * sizeof(double) );

pesos = (double *) malloc( n * sizeof(double));
aux1 = (double *) malloc( n * sizeof(double));
aux2 = (double *) malloc( n * sizeof(double));
res = (double *) malloc( n * sizeof(double));
beta = (double *) malloc( q * sizeof(double));
beta_old = (double *) malloc( q * sizeof(double));

tmp1 = (double **) malloc( q * sizeof(double *));
tmp2 = (double **) malloc( q * sizeof(double *));
tmp5 = (double **) malloc( q * sizeof(double *));
tmp3 = (double *) malloc( q * sizeof(double));
tmp4 = (double *) malloc( q * sizeof(double));

for(j=0;j<q;j++) {
	tmp1[j] = (double *) malloc( q * sizeof(double) );
	tmp2[j] = (double *) malloc( (q + 1) * sizeof(double) );
	tmp5[j] = (double *) malloc( q * sizeof(double) );
}

corte = 10 * (*eps);
ktt = *kt;

for(i=0;i<n; i++) {
     pesos[i] = kepan( (x[i] - *punto) / *ventanas )  / prob[i];
     for(j=0; j<q; j++) zz[i][j] = *(z + j*n + i);
};
reset_mat(tmp1, q, q);
reset_vec(tmp3, q);
reset_mat(tmp5, q, q);
reset_vec(beta, q);
reset_vec(beta_old, q);
for(j=0; j<q; j++) beta[j] = beta_initial[j];
while( (corte > (*eps)) && (it< (*maxit)) ){
	reset_mat(tmp2, q, q+1);
	reset_vec(tmp4, q);
	for(j=0; j<q; j++) beta_old[j] = beta[j];
	for(i=0;i<n;i++){
		res[i] = ( y[i] - vecprime_vec(zz[i], beta, q) ) / *sigmahat;
		aux1[i] = pesos[i] * psi_tukey_w(res[i],ktt) * y[i];
		aux2[i] = pesos[i] * psi_tukey_w(res[i],ktt);
		scalar_vec(zz[i], aux2[i], tmp3, q);
		vec_vecprime(tmp1, tmp3, zz[i], q);
		sum_mat(tmp1, tmp2, tmp2, q, q);
		scalar_vec(zz[i], aux1[i], tmp3, q);
		sum_vec(tmp3, tmp4, tmp4, q);
	};
	for(j=0;j<q;j++) tmp2[j][q] = tmp4[j];
	if( lu(tmp2, &q, beta) == 1) {
		it = *maxit; // System became singular
		for(j=0; j<q; j++) beta[j] = NA_REAL;
	};
	// if( lu(tmp2, &q, beta) == 1) it = *maxit; // System became singular
	corte = l2dist(beta, beta_old, q) / (norma_dos(beta_old, q) + *eps);
	it = it + 1;
};
for(j=0; j<q; j++) salida[j] = beta[j];
free(pesos); free(aux1); free(aux2); free(res);
for(i=0; i<n; i++) free(zz[i]); free(zz);
for(j=0;j<q;j++) {
	free( tmp1[j] ); free( tmp2[j] ); free( tmp5[j] );
}
free(tmp1); free(tmp2); free(tmp3); free(tmp4); free(tmp5);
free(beta); free(beta_old);
return;
}



void lu_R(double *a, double *b, int *kk, double *y) {
int i,j, k = *kk;
double **m1;
int lu(double **a,int *P, double *x);
m1 = (double **) malloc( k * sizeof(double *) );
for(j=0;j<k;j++)
	m1[j] = (double *) malloc( (k+1) * sizeof(double) );
for(i=0;i<k;i++) {
	for(j=0;j<k;j++)
		m1[i][j] = *(a + j*k + i);
	m1[i][k] = b[i];
};
lu(m1, kk, y);
for(j=0;j<k;j++) free(m1[j]);
free(m1);
return;
}





int lu(double **a,int *P, double *x)
{
int *pp,p;
register int i,j,k;
double *kk,s;
p = *P;
if ((pp = (int *) malloc(p*sizeof(int)))==NULL)
	{ printf("\nNot enough memory in LU\n");
	  exit(1); }
/* pp vector storing the permutations */
for(j=0;j<p;j++)   /* cols */
{ pp[j]=j;
  for(i=j;i<p;i++)   /* filas */
	if ( fabs( a[i][j] ) > fabs( a[pp[j]][j] ) )
		pp[j]=i;
  if ( pp[j] != j )       /* permuto las filas cambiando los punt */
	{ kk=a[j];
	  a[j]=a[pp[j]];
	  a[pp[j]]=kk;
	};
  /* salida si el sistema resulta singular (det=0)
   * se detecta si el pivote (j,j) es cero  */
/*  if ( a[j][j] == 0 ) {   free(pp);
				return(1);
				}; */
    if ( fabs(a[j][j]) < TOL_INVERSE ) {   free(pp);
				return(1);
				};
  for(k=(j+1);k<p;k++)
	a[k][j] = a[k][j] / a[j][j];
  for(k=(j+1);k<p;k++)
	for(i=(j+1);i<p;i++)
		a[k][i] = a[k][i] - a[k][j] * a[j][i];

};    /* cierra el for de j */
for(i=0;i<p;i++)
	{ s=0.0;
	  for(j=0;j<i;j++)
	    s += a[i][j] * x[j];
	    x[i] = a[i][p] - s;          /* y[i]=a[i][p] */
	};
for(i=(p-1);i>=0;i--)
	{ s=0;
	  for(j=(i+1);j<p;j++)
	    s += a[i][j] * x[j];
	  x[i] = (x[i] - s) / a[i][i];
	  };
free(pp);
return(0);
}


void sum_mat(double **a, double **b, double **c, int n, int m)
{
register int i,j;
for(i=0;i<n;i++)
	for(j=0;j<m;j++)
		c[i][j] = a[i][j] + b[i][j];
}

void vec_vecprime(double **a, double *v1, double *v2, int n)
{
register int i,j;
for(i=0;i<n;i++)
	for(j=0;j<n;j++)
		a[i][j] = v1[i] * v2[j];/* could take advantage of symmetry */
}

void scalar_mat(double **a, double b, double **c, int n, int m)
{
register int i,j;
for(i=0;i<n;i++)
        for(j=0;j<m;j++)
	c[i][j]  = b * a[i][j];
}

void scalar_vec(double *a, double b, double *c, int n)
{
register int i;
for(i=0;i<n;i++)
	c[i]  = b * a[i];
}

double vecprime_vec(double *a, double *b, int n)
{
register int i;
double s = 0.0;
for(i=0;i<n;i++) s += a[i] * b[i];
return(s);
}

void sum_vec(double *a, double *b, double *c, int n)
{
register int i;
for(i=0;i<n;i++) c[i] = a[i] + b[i];
}

void dif_vec(double *a, double *b, double *c, int n)
{
register int i;
for(i=0;i<n;i++) c[i] = a[i] - b[i];
}

void dif_mat(double **a, double **b, double **c, int n, int m)
{
register int i,j;
for(i=0;i<n;i++)
	for(j=0;j<m;j++) c[i][j] = a[i][j] - b[i][j];
}

void mat_vec(double **a, double *b, double *c, int n, int m)
{
register int i,j;
for(i=0;i<n;i++)
	for(c[i]=0,j=0;j<m;j++) c[i] += a[i][j] * b[j];
}

void mat_mat(double **a, double **b, double **c, int n,
		int m, int l)
{
register int i,j,k;
for(i=0;i<n;i++)
	for(j=0;j<l;j++) {
	c[i][j] = 0;
	for(k=0;k<m;k++) c[i][j] += a[i][k] * b[k][j];
	};
}

/* void disp_vec(double *a, int n)
{
register int i;
Rprintf("\n");
for(i=0;i<n; i++) Rprintf("%lf ",a[i]);
Rprintf("\n");
}

void disp_mat(double **a, int n, int m)
{
register int i,j;
for(i=0;i<n;i++) {
Rprintf("\n");
for(j=0;j<m;j++) Rprintf("%10.8f ",a[i][j]);
};
Rprintf("\n");
}
*/

int inverse(double **a, double **b, int n)
{
int lu(double **, int *, double *);
void mat_vec(double **, double *, double *, int, int);
register int i,j,k;
double **c, *e;
c = (double **) malloc( n * sizeof(double *));
e = (double *) malloc( n * sizeof(double));
for(i=0;i<n;i++) c[i] = (double *) malloc ( (n+1) * sizeof(double) );
for(i=0;i<n;i++) {   /* i-th column */

for(j=0;j<n;j++)
	for(k=0;k<n;k++) c[j][k] = a[j][k];
for(j=0;j<i;j++) c[j][n] = 0.0;
c[i][n] = 1.0;
for(j=i+1;j<n;j++) c[j][n] = 0.0;
if( lu(c,&n,e) == 1) {
	for(i=0;i<n;i++) free(c[i]);
	free(c);free(e);
	return(1);
	};
for(j=0;j<n;j++) b[j][i] = e[j] ;
};
for(i=0;i<n;i++) free(c[i]);
free(c);free(e);
return(0);
}

void reset_mat(double **a, int n, int m)
{
register int i,j;
for(i=0;i<n;i++)
	for(j=0;j<m;j++)
		a[i][j] = 0.0;
}

void reset_vec(double *a, int n)
{
register int i;
for(i=0;i<n;i++) a[i] = 0.0;
}


double median(double *x, int n)
{
double kthplace(double *,int,int);
double *aux,t;
register int i;
if ( (aux = (double *) malloc (n*sizeof(double)) )==NULL)
	{printf("\nNot enought memory in median\n"); exit(1); };
for(i=0;i<n;i++) aux[i]=x[i];
if ( (n/2) == (double) n / 2 )
	t = ( kthplace(aux,n,n/2) + kthplace(aux,n,n/2+1) ) / 2 ;
else	t = kthplace(aux,n, n/2+1 ) ;
free(aux);
return(t);
}

double median_abs(double *x, int n)
{
double kthplace(double *,int,int);
double *aux,t;
register int i;
if ( (aux = (double *) malloc (n*sizeof(double)) )==NULL )
	{ printf("\nNot enought memory in med_abs\n");exit(1);};
for(i=0;i<n;i++) aux[i]=fabs(x[i]);
if ( (n/2) == (double) n / 2 )
	t = ( kthplace(aux,n,n/2) + kthplace(aux,n,n/2+1) ) / 2 ;
else 	t = kthplace(aux,n, n/2+1 ) ;
free(aux);
return(t);
}


double kthplace(double *a, int n, int k)
{
int jnc,j;
int l,lr;
double ax,w;
k--;
l=0;
lr=n-1;
while (l<lr)
	{ ax=a[k];
	  jnc=l;
	  j=lr;
	  while (jnc<=j)
		{ while (a[jnc] < ax) jnc++;
		  while (a[j] > ax) j--;
		  if (jnc <= j)
			{ w=a[jnc];
			  a[jnc]=a[j];
			  a[j]=w;
			  jnc++;
			  j--;
			};
		};
	  if (j<k) l=jnc;
	if (k<jnc) lr=j;
	};
return(a[k]);
}


void huber_pos(int *nrow, double *y, double *muhat_initial,
 double *eps, double *sigmahat, double *prob,
 double *kh, int *maxit, double *salida) {

double psi_huber_w(double, double);
double suma_vec(double*, int);
register int i;
int n = *nrow, it=0;
double corte, muhat, muold, khh;
double *aux1, *aux2, *res;

aux1 = (double *) malloc( n * sizeof(double));
aux2 = (double *) malloc( n * sizeof(double));
res = (double *) malloc( n * sizeof(double));
corte = 10 * (*eps);
khh = *kh;
muhat = *muhat_initial;

while( (corte > (*eps)) && (it< (*maxit)) ){
	for(i=0;i<n;i++){
		res[i] = ( y[i] - muhat) / *sigmahat;
		aux1[i] = psi_huber_w(res[i],khh) * y[i] / prob[i];
		aux2[i] = psi_huber_w(res[i],khh) / prob[i];
	};
	muold = muhat;
	muhat = suma_vec(aux1,n) / suma_vec(aux2,n);
	corte = fabs(muold-muhat) / (fabs(muold) + *eps);
	it = it + 1;
};
*salida = muhat;
free(aux1); free(aux2); free(res);
return;
}

void tukey_pos(int *nrow, double *y, double *muhat_initial,
 double *eps, double *sigmahat, double *prob,
 double *kt, int *maxit, double *salida) {

double psi_tukey_w(double, double);
double suma_vec(double*, int);
register int i;
int n = *nrow, it=0;
double corte, muhat, muold, ktt;
double *aux1, *aux2, *res;

aux1 = (double *) malloc( n * sizeof(double));
aux2 = (double *) malloc( n * sizeof(double));
res = (double *) malloc( n * sizeof(double));
corte = 10 * (*eps);
ktt = *kt;
muhat = *muhat_initial;

while( (corte > (*eps)) && (it< (*maxit)) ){
	for(i=0;i<n;i++){
		res[i] = ( y[i] - muhat) / *sigmahat;
		aux1[i] = psi_tukey_w(res[i],ktt) * y[i] / prob[i];
		aux2[i] = psi_tukey_w(res[i],ktt) / prob[i];
	};
	muold = muhat;
	muhat = suma_vec(aux1,n) / suma_vec(aux2,n);
	corte = fabs(muold-muhat) / (fabs(muold) + *eps);
	it = it + 1;
};
*salida = muhat;
free(aux1); free(aux2); free(res);
return;
}

