#include <R.h>
#include <math.h>
#include <gsl/gsl_math.h>

// utilities
double *create_real_vector(int size){
  double *vec = (double *) Calloc(size,double);
  return vec;
}

void destroy_real_vector(double *vector, int size){
  Free(vector);
}

double **create_real_matrix(int rows, int cols){
  int i;
  double **mat = (double **) Calloc(rows,double*);
  for(i = 0; i < rows; i++) mat[i] = (double *) Calloc(cols,double);
  return mat;
}

double **create_and_copy_real_matrix(int rows, int cols, double *data){
  int i,j;
  double **mat = (double **) Calloc(rows,double*);
    mat[i] = (double *) Calloc(cols,double);
    for(i = 0; i < rows; i++){
      for(j=0;j<cols;++j) mat[i][j] = data[j * rows + i];
  }
  return mat;
}

void real_matrix_copy(double **mat, int d1, int d2, double *data){
  int i, j;
  for( i=0; i<d1; ++i)
    for( j=0; j<d2; ++j)
      data[i+d1*j] = mat[i][j];
}

void destroy_real_matrix(double **matrix, int rows, int cols){
  int i;
  for(i = 0; i < rows; i++){ Free(matrix[i]); }
  Free(matrix);
}

double **create_real_array3d(int d1, int d2, int d3){
  int i,j;
  double ***array = (double **) Calloc(d1,double**);
  for (i=0; i < d1; ++i){
    array[i] = Calloc(d2,double*);
    for(j=0; j<d2;++j){
      array[i][j] = Calloc(d3,double);
    }
  }
  return array;
}

void real_array3d_copy(double ***array, int d1, int d2, int d3, double *data){
  int i, j, k;
  for( i=0; i<d1; ++i)
    for( j=0; j<d2; ++j)
      for( k=0; k<d3; ++k )
        data[i+d1*j+d1*d2*k] = array[i][j][k];
}

void destroy_real_array3d(double ***array, int d1, int d2, int d3){
  int i, j;
  for (i=0; i < d1; ++i) {
    for (j=0; j < d2; ++j) Free(array[i][j]);
    Free(array[i]);
  }
  Free(array);
}

void chol(double **L, double **M, int n) {
  int i, j, k;
  double s;
  for(i = 0; i < n; i++) {
    for(j = 0; j < (i+1); j++) {
      s = 0;
      for( k = 0; k < j; k++) s += L[i][k] * L[j][k];
      L[i][j] = (i == j) ? sqrt(M[i][i] - s) : (1.0 / L[j][j] * (M[i][j] - s));
    }
  }
}

// update of S = a*S + b*x*x'
void rank_one_chol_up(double **L_new, double **L_old, double *x, int n, double a, double b, double *work) {
  int i,j;
  double r, c, s;

  memcpy(work,x,sizeof(double)*n);
  for(i=0;i<n;++i) work[i] = sqrt(b)*work[i]; // this could be avoided maybe?

  for(i=0;i<n;++i){
    r = sqrt(a*L_old[i][i]*L_old[i][i] + work[i]*work[i]);
    c = r /( sqrt(a)*L_old[i][i] );
    s = (work[i]) / (sqrt(a)*L_old[i][i]);
    L_new[i][i] = r;
    for(j=i+1;j<n;++j){
      L_new[j][i] = (sqrt(a)*L_old[j][i] + s*work[j]) / c;
      work[j] = c*work[j] - s*L_new[j][i];
    }
  }
}

// Update of the form S = a*S + WW*
void rank_r_chol_up(double **L_new, double **L_old, double **W, int n, double a, double b, double **work, int r) {
  int i,j,p;
  double alpha_bar;
  double *alpha, *d, *gamma;

  work = create_and_copy_real_matrix(n,n,W);

  for( i = 1; i < r; ++i){
    alpha[i] = 1;
  }
  for( j = 1; j < n; ++j){
    for( i = 1; i < r; ++i){
      alpha_bar = alpha[i] + work[j][i]*work[j][i]/d[j];
      d[j] = d[j] * alpha_bar;
      gamma[i] = work[j][i]/d[j];
      d[j] = d[j]/alpha[i];
      alpha[i] = alpha_bar;
    }
    for( p = j; p < n; ++p){
      for( i = 1; i < r; ++i){
        work[p][i] = work[p][i] - work[j][i] * L_old[p][j];
        L_new[p][j] = L_old[p][j] + gamma[i] * work[p][i];
      }
    }
  }
}

void fwdinv(double **L_inv, double **L, int n){
  int i,j,k;
  for( k=0; k<n; ++k ){
    for( i=0; i<n; ++i ){
      L_inv[i][k] = (i==k?1:0);
      for( j=0; j<i; ++j ) L_inv[i][k] -= L[i][j] * L_inv[j][k];
      L_inv[i][k] = L_inv[i][k]/L[i][i];
    }
  }
}

void matvec(double *y, double **A, double *x, int n){
  int i,j;
  for(i=0;i<n;++i){
    y[i] = A[i][0]*x[0];
    for(j=1;j<n;++j){
      y[i] += A[i][j]*x[j];
    }
  }
}


// Bivariate DCC(1,1) Model Filter
void bidcc_filter(int *status, double *rho, double* eps, double *loglik, double *param, double *_y, int *T, int n){

  int t, i;
  double logden;
  double alpha,beta;
  double **Q, **y;
  double rho_bar;
  *loglik = 0;

  // sanity check
  if( !finite(param[0]) || !finite(param[1]) ){
    *loglik = -HUGE_VAL;
    return;
  }

  alpha = param[0];
  beta  = param[1];

  // check constraints
  if( alpha <= 1e-5 || beta < 0 || (alpha+beta)>1 ){
    *loglik = -HUGE_VAL;
    return;
  }

  // allocate
  Q = create_real_matrix(*T,3);
  y = create_and_copy_real_matrix(*T,2,_y);



  // loop
  *loglik = 0;
  for( t=1; t<*T; ++t ){
    Q[t][0] = (1-alpha-beta)         + alpha*y[t-1][0]*y[t-1][0] + beta*Q[t-1][0];
    Q[t][1] = (1-alpha-beta)         + alpha*y[t-1][1]*y[t-1][1] + beta*Q[t-1][1];
    Q[t][2] = rho_bar*(1-alpha-beta) + alpha*y[t-1][0]*y[t-1][1] + beta*Q[t-1][2];
    rho[t]  = Q[t][2]/sqrt(Q[t][0]*Q[t][1]);

    logden  = -0.5*log(2*PI) - 0.5*log(1-rho[t]*rho[t]) - 0.5*(y[t][0]*y[t][0]+y[t][1]*y[t][1]-2*y[t][0]*y[t][1]*rho[t])/ (1.0-rho[t]*rho[t]);

    if( finite(logden) ){
      *loglik += logden;
    }
    else{
      Rprintf("problem at time %d\n",t);
    }

  }

  // safeguard
  if( !isfinite(*loglik) ){
    *loglik = -HUGE_VAL;
  }

  // cleanup
  destroy_real_matrix(Q,*T,3);
  destroy_real_matrix(y,*T,2);
}

/*
// n-variate DCC(1,1) Model Filter
void dcc_filter(int *status, double *rho, double *eps, double *loglik, double *param, double *_y, int *T, int *N, double **rho_bar){

  int t;
  double logden;
  double alpha,beta;
  double **Q, **y;
  *loglik = 0;

  // sanity check
  if( !finite(param[0]) || !finite(param[1]) ){
    *loglik = -HUGE_VAL;
    return;
  }

  alpha = param[0];
  beta  = param[1];

  // check constraints
  if( alpha <= 1e-5 || beta < 0 || (alpha+beta)>1 ){
    *loglik = -HUGE_VAL;
    return;
  }

  // allocate
  Q = create_real_array3d(*T,*N,*N);
  y = create_and_copy_real_matrix(*T,*N,_y);

  // init
  Q[0] = **rho_bar

  // loop
  *loglik = 0;
  for( t=1; t<*T; ++t ){
    Q[t][0] = (1-alpha-beta)         + alpha*y[t-1][0]*y[t-1][0] + beta*Q[t-1][0];
    Q[t][1] = (1-alpha-beta)         + alpha*y[t-1][1]*y[t-1][1] + beta*Q[t-1][1];
    Q[t][2] = rho_bar*(1-alpha-beta) + alpha*y[t-1][0]*y[t-1][1] + beta*Q[t-1][2];
    rho[t]  = Q[t][2]/sqrt(Q[t][0]*Q[t][1]);

    logden  = -0.5*log(2*PI) - 0.5*log(1-rho[t]*rho[t]) - 0.5*(y[t][0]*y[t][0]+y[t][1]*y[t][1]-2*y[t][0]*y[t][1]*rho[t])/ (1.0-rho[t]*rho[t]);

    if( finite(logden) ){
      *loglik += logden;
    }
    else{
      Rprintf("problem at time %d\n",t);
    }

  }

  // safeguard
  if( !isfinite(*loglik) ){
    *loglik = -HUGE_VAL;
  }

  // cleanup
  destroy_real_matrix(Q,*T,3);
  destroy_real_matrix(y,*T,2);
}
*/
