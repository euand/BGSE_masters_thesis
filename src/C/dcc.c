
#include <R.h>
#include <math.h>

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
  for(i = 0; i < rows; i++){
    mat[i] = (double *) Calloc(cols,double);
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

/*
// update of S = a*S + b*x*x'
void chol_up(double **L_new, double **L_old, double *x, int n, double a, double b, double *work) {
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
*/

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

// filters
// Gaussian GARCH(1,1) Filter
void garch_filter(int *status, double *sigma2, double* eps, double *loglik, double *param, double *y, int *T){

  double logden;
  double omega, alpha, beta;
  int t;

  omega = param[0];
  alpha = param[1];
  beta  = param[2];

  // check constraints
  if( alpha <= 1e-6 || beta < 0 || omega<=0 || (alpha+beta)>1 ){
    *loglik = -HUGE_VAL;
    return;
  }

  // init
  sigma2[0]=0;
  for( t=0; t<10; ++t ){ sigma2[0] += y[t]*y[t]; }
  sigma2[0] /= 10;
  eps[0] = y[0]/sqrt( sigma2[0] );
  *loglik = 0;

  // loop
  for( t=1 ; t<*T ; ++t ){
    sigma2[t] = omega + alpha * y[t-1]*y[t-1] + beta * sigma2[t-1];
    eps[t]    = y[t]/sqrt( sigma2[t] );

    logden    = -0.5 *log(2*PI) -0.5*log( sigma2[t] ) -0.5*(y[t]*y[t])/sigma2[t];
    *loglik   += logden;
  }

  // safeguard
  if( !isfinite(*loglik) ){
    *loglik = -HUGE_VAL;
  }

}

// Gaussian GARCH(1,1) numerical OPG
void garch_opg_num(int *status, double *OPG, double *param, double *y, int *T, double *eps){

  double logden_plus, logden_minus;
  double param_plus[4], param_minus[4];
  double dlogden[4];
  double omega, alpha, beta;
  int t,i,j,p;

  /*
  // init
  sigma2[0]=0;
  for( t=0; t<10; ++t ){
  sigma2[0] += y[t]*y[t];
  }
  sigma2[0] /= 10;

  // loop
  for( t=1; t<*T; ++t ){
  for( p=0; p<4; ++p ){

  memcpy(param_plus ,param,4*sizeof(double));
  param_minus[p] += eps[p];
  omega = param_plus[0]; alpha = param_plus[1]; beta  = param_plus[2];
  sigma2[t] = omega + alpha * y[t-1]*y[t-1] + beta * sigma2[t-1];
  logden_plus = -0.5 *log(2*PI) -0.5*log( sigma2[t] ) -0.5*(y[t]*y[t])/sigma2[t];

  memcpy(param_minus,param,4*sizeof(double));
  param_plus[p] += eps[p];
  omega = param_minus[0]; alpha = param_minus[1]; beta  = param_minus[2];
  sigma2[t] = omega + alpha * y[t-1]*y[t-1] + beta * sigma2[t-1];
  logden_minus = -0.5 *log(2*PI) -0.5*log( sigma2[t] ) -0.5*(y[t]*y[t])/sigma2[t];

  dlogden[p] = (logden_plus-logden_minus)/(2*eps[p]);

  }

  for( i=0; i<4; ++i ){
  for( j=0; j<=i; ++j){
  dlogden[i]*dlogden[j];
  {
  }
  }
  */

}

// Gaussian GARCH(1,1) numerical Hessian
void garch_hessian_num(int *status, double *OPG, double *param, double *y, int *T, double *eps){

}

// TARCH(1,1) Model Filter
void tarch_filter(int *status, double *sigma2, double* eps, double *loglik, double *param, double *y, int *T){

  double logden;
  double omega, alpha, gamma, beta;
  int t;

  omega = param[0];
  alpha = param[1];
  gamma = param[2];
  beta  = param[3];

  // check constraints
  if( alpha <= 0 || beta < 0 || omega<0 || (alpha+beta)>1 ){
    *loglik = -HUGE_VAL;
    return;
  }

  // init
  sigma2[0]=0;
  for( t=0; t<10; ++t ){ sigma2[0] += y[t]*y[t]; }
  sigma2[0] /= 10;
  eps[0] = y[0]/sqrt( sigma2[0] );

  // loop
  *loglik = 0;
  for( t=1 ; t<*T ; ++t ){
    sigma2[t] = omega + alpha * y[t-1]*y[t-1] + gamma * y[t-1]*y[t-1]*((double)(y[t-1]<0)) + beta * sigma2[t-1];
    eps[t]    = y[t]/sqrt( sigma2[t] );

    logden    = -0.5 *log(2*PI) -0.5*log(sigma2[t]) -0.5*(y[t]*y[t])/sigma2[t];
    *loglik   += logden;
  }
}

// Bivariate DCC(1,1) Model Filter
void bidcc_filter(int *status, double *rho, double* eps, double *loglik, double *param, double *_y, int *T){

  int t;
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

  // init
  rho_bar = 0;
  for( t=0; t<*T; ++t ){ rho_bar += y[t][0]*y[t][1]; }
  rho_bar /= *T;

  Q[0][0] = 1;
  Q[0][1] = 1;
  Q[0][2] = rho_bar;
  rho[0]  = rho_bar;

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

// MEWMA Model Filter
void mewma_filter(int *status, double *_s, double* _eps, double *loglik, double *param, double *_y, int *T, int *N){

  int t,i,j;
  double logden;
  double lambda;
  double ***S, **y, **eps;
  double **Sig;
  double *work1, **work2;
  double rho_bar;
  *loglik = 0;

  lambda = param[0];

  // check constraints
  if( lambda <= 1e-5 || lambda>1 ){
    *loglik = -HUGE_VAL;
    return;
  }

  // allocate
  S     = create_real_array3d(*T,*N,*N);
  y     = create_and_copy_real_matrix(*T,*N,_y);
  eps   = create_and_copy_real_matrix(*T,*N,_eps);
  work1 = create_real_vector(*N);
  work2 = create_real_matrix(*N,*N);

  // init
  for( i=0; i<*N; ++i){
    for( j=0; j<=i; ++j ){
      work2[i][j]=0;
      for( t=0; t<*T; ++t ){ work2[i][j] += y[t][i]*y[t][j]; }
      work2[i][j] /= *T;
      work2[j][i] = work2[i][j];
    }
  }
  chol(S[0],work2,*N);
  *loglik = 0;

  // loop
  for( t=1; t<*T; ++t ){

    //chol_up(S[t],S[t-1],y[t-1],*N,lambda,1.0-lambda,work1);
    fwdinv(work2,S[t],*N);
    matvec(eps[t],work2,y[t],*N);

    logden = -0.5*(*N)*log(2*PI);
    for(i=0;i<*N;++i) logden += -log( S[t][i][i] )-0.5*eps[t][i]*eps[t][i];

    *loglik += logden;
  }

  // safeguard
  if( !isfinite(*loglik) ){
    *loglik = -HUGE_VAL;
  }

  // copy results
  real_array3d_copy(S,*T,*N,*N,_s);
  real_matrix_copy(eps,*T,*N,_eps);

  // cleanup
  destroy_real_array3d(S,*T,*N,*N);
  destroy_real_matrix(y,*T,*N);
  destroy_real_matrix(eps,*T,*N);
  destroy_real_vector(work1,*N);
  destroy_real_matrix(work2,*N,*N);

}

// update of S = a*xx' + b*S
void chol_up(double **L, double *x, int n, double a, double b, double *work) {
  int i,j;
  double r, c, s;

  memcpy(work,x,sizeof(double)*n);
  for(i=0;i<n;++i) work[i] = sqrt(a/b)*work[i]; // this could be avoided maybe?

  for(i=0;i<n;++i){
    r = sqrt(L[i][i]*L[i][i] + work[i]*work[i]);
    c = r/L[i][i];
    s = work[i] / L[i][i];
    L[i][i] = r;
    for(j=i+1;j<n;++j){
      L[j][i] = (L[j][i] + s*work[j]) / c;
      work[j] = c*work[j] - s*L[j][i];
    }
  }
  for(i=0;i<n;++i){
    for(j=0;j<=i;++j){
      L[i][j] *= sqrt(b);
    }
  }
}

void solve(double **L, double *b, double *y, double *x, int n){
  int i,j,k;

  y[0] = b[0] / L[0][0];
  for( i=0; i<n; ++i ){
    y[i] = b[i];
    for(j=0; j<i; ++j) y[i] -= y[j] * L[i][j];
    y[i] = y[i] / L[i][i];
  }
  x[n-1] = y[n-1] / L[n-1][n-1];
  for( i=(n-2); i>-1; --i ){
    x[i] = y[i];
    for(j=(n-1); j>i; --j) x[i] -= x[j] * L[j][i];
    x[i] = x[i] / L[i][i];
  }
}

// n-variate DCC Model Filter
void dcc_filter(int *status, double *_L, double *_eps, double *loglik, double *param, int *T, int *N){

  int i,j,t;
  double logden;
  double alpha, beta, L_tilde_denom;
  double *y, *x;
  double *work1;
  double **L, **L_new, **eps, **L_tilde;

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
  work1 = create_real_vector(*N);
  x = create_real_vector(*N);
  y = create_real_vector(*N);
  L = create_and_copy_real_matrix(*N,*N, _L);
  L_tilde = create_and_copy_real_matrix(*N,*N, _L);
  eps   = create_and_copy_real_matrix(*T, *N, _eps);

  // First part of log-likelihood: log determinant of Rt
  logden = 0;
  L_tilde_denom = 0;
  // Calculate Q^(-1/2) * L
  for(i = 0; i < *N; ++i){
    for(j = 0; j <= i; ++j) L_tilde_denom += L[i][j] * L[i][j];
    for(j = 0; j <= i; ++j) L_tilde[i][j] = L[i][j] / sqrt(L_tilde_denom);
    L_tilde_denom = 0;
  }

  // Add log determinant
  for(i = 0; i < *N; ++i) *loglik += log(L_tilde[i][i]);
  // Add likelihoods from each time step
  //*loglik += logden;
  for(t=1; t<*T; ++t ){

    logden = 0;
    // Update L
    // update of S = a*S + b*x*x'
    //void chol_up(double **L, double *x, int n, double a, double b, double *work) {
    //work1 = eps[t];
    chol_up(L, eps[t], *N, alpha, beta, work1);

    // Calculate Q^(-1/2) * L
    for(i = 0; i < *N; ++i){
      for(j = 0; j <= i; ++j) L_tilde_denom += L[i][j] * L[i][j];
      for(j = 0; j <= i; ++j) L_tilde[i][j] = L[i][j] / sqrt(L_tilde_denom);
      L_tilde_denom = 0;
    }

    for(i = 0; i < *N; ++i) *loglik += log(L_tilde[i][i]);

    //if( finite(logden) ){
    //  *loglik += logden;
    //}
    //else{
    //  Rprintf("problem at time %d\n",t);
    //}

    solve(L_tilde, eps[t], y, x, *N);

    logden = 0;

    for(i = 0; i < *N; ++i) *loglik += y[i] * eps[t][i];


  }

  // safeguard
  if( !isfinite(*loglik) ){
    *loglik = -HUGE_VAL;
  }

  // cleanup
  //destroy_real_matrix(L,*T,*T);
  //destroy_real_matrix(eps,*T,*N);
}
//
