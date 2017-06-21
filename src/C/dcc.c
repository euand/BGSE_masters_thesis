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

// update of S = a*xx' + b*S
void chol_up(double **L, double *x, int n, double a, double b, double *work) {
  int i,j;
  double r, c, s;

  memcpy(work,x,sizeof(double)*n);
  for(i=0;i<n;++i) work[i] = sqrt(a/b)*work[i]; // this could be avoided maybe?

  for(i=0;i<n;++i){
    r = sqrt(L[i][i]*L[i][i] + work[i]*work[i]);
    c = r/L[i][i];
    s = (work[i]) / L[i][i];
    L[i][i] = r;
    for(j=i+1;j<n;++j){
      L[j][i] = (L[j][i] + s*work[j]) / c;
      work[j] = c*work[j] - s*L[j][i];
    }
  }
  for(i=0;i<n;++i){
    for(j=0;j<i;++j){
      L[i][j] *= sqrt(b);
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

double solve(double **L, double *b, double *y, double *x, int n){
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
  return *x;
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

// n-variate DCC Model Filter
void dcc_filter(int *status, double **omega, double **eps, double *loglik, double *param, double *_y, int *T, int *N){

  int i,t;
  double logden1, logden2;
  double alpha,beta;
  double *y, *x;
  double *work1;
  double **work2, **Q;

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
  y = create_and_copy_real_matrix(*T,*N,_y);
  work1 = create_real_vector(*N);
  work2 = create_real_matrix(*N,*N);
  x = create_real_vector(*N);
  *loglik = 0;

  // init
  chol(Q, work2, *N);

  // First part of log-likelihood: log determinant of Rt
  logden1 = 0;

  for(i = 0; i < *N; ++i){
    logden1 += log(work2[i][i]);
  }

  // loop
  *loglik = 0;
  for(t=1; t<*T; ++t ){

    if( finite(logden1) ){
      *loglik += logden1;
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
}
