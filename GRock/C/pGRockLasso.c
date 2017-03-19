/*----------------------------------------------------------
 * Solve a distributed lasso problem, i.e.,
 *
 *   minimize lambda ||x||_1 + 0.5 * ||Ax - b||_2^2
 *
 * The implementation uses MPI for distributed communication
 * and the GNU Scientific Library (GSL) for math.
 *
 * Compile: make
 * run code: mpiexec -n 1 ./pGRockLasso $dataDirectory$
 *
 * Copyright (c)     Zhimin Peng @ Math, UCLA
 * Created Date:     01/11/2013
 * Modified Date:    05/02/2013
 *                   05/02/2014, 
 *                   add dynamic update P, 
 *                   fixed some bugs
 * --------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "mmio.h"
#include <mpi.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort_vector.h>

struct value {
  int    ID;
  double data;
};

typedef int (*compfn)(const void*, const void*);
int compare(struct value *, struct value *);
void shrink(gsl_vector *v, double sigma);
double objective(gsl_vector *x, double lambda, gsl_vector *z, int N);
void pointwise(gsl_vector *a, gsl_vector *b, double n);
void abs_vector(gsl_vector *a, double n);
int main(int argc, char **argv) {

  const int MAX_ITER  = 20;
  const double TOL = 1e-12;
  
  int rank;
  int size;
  int P = 8; // number of blocks to update P <= size

  /* -----------------------------------
     mode controls the selection schemes, 
       mode =0, fixed P
       mode =1, dynamic update P
     ----------------------------------*/
  int mode=1; // number of processors used to update each time
  double lambda = 0.1;
  srand (time(NULL));
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Determine current running process
  MPI_Comm_size(MPI_COMM_WORLD, &size); // Total number of processes
  
  // data directory (you need to change the path to your own data directory)
  char* dataCenterDir = "../Data/Gaussian";
  char* big_dir;
  if(argc==2)
    big_dir = argv[1];
  else
    big_dir = "big1";

  /* Read in local data */
  
  FILE *f, *test;
  int m, n, j;
  int row, col;
  double entry, startTime, endTime;
  double total_start_time, total_end_time;
  /*
   * Subsystem n will look for files called An.dat and bn.dat
   * in the current directory; these are its local data and do not need to be
   * visible to any other processes. Note that
   * m and n here refer to the dimensions of the *local* coefficient matrix.
   */
  
  /* ------------
     Read in A 
     ------------*/
  if(rank ==0){
    printf("=============================\n");
    printf("|    Start to load data!     |\n");
    printf("=============================\n");
  }
  char s[100];
  sprintf(s, "%s/%s/A%d.dat",dataCenterDir,big_dir, rank + 1);
  printf("[%d] reading %s\n", rank, s);
  f = fopen(s, "r");
  if (f == NULL) {
    printf("[%d] ERROR: %s does not exist, exiting.\n", rank, s);
    exit(EXIT_FAILURE);
  }
  mm_read_mtx_array_size(f, &m, &n);
  gsl_matrix *A = gsl_matrix_calloc(m, n);
  for (int i = 0; i < m*n; i++) {
    row = i % m;
    col = floor(i/m);
    fscanf(f, "%lf", &entry);
    gsl_matrix_set(A, row, col, entry);
  }
  fclose(f);
  
  /* ------------
      Read in b 
     -------------*/
  sprintf(s, "%s/%s/b.dat", dataCenterDir, big_dir);
  printf("[%d] reading %s\n", rank, s);
  f = fopen(s, "r");
  if (f == NULL) {
    printf("[%d] ERROR: %s does not exist, exiting.\n", rank, s);
    exit(EXIT_FAILURE);
  }
  mm_read_mtx_array_size(f, &m, &n);
  gsl_vector *b = gsl_vector_calloc(m);
  for (int i = 0; i < m; i++) {
    fscanf(f, "%lf", &entry);
    gsl_vector_set(b, i, entry);
  }
  fclose(f);
  
  /* ------------
     Read in xs 
     ------------*/
  sprintf(s, "%s/%s/xs%d.dat", dataCenterDir, big_dir, rank + 1);
  printf("[%d] reading %s\n", rank, s);
  f = fopen(s, "r");
  if (f == NULL) {
    printf("[%d] ERROR: %s does not exist, exiting.\n", rank, s);
    exit(EXIT_FAILURE);
  }
  mm_read_mtx_array_size(f, &m, &n);
  gsl_vector *xs = gsl_vector_calloc(m);
  
  for (int i = 0; i < m; i++) {
    fscanf(f, "%lf", &entry);
    gsl_vector_set(xs, i, entry);
  }
  fclose(f);
  
  m = A->size1;
  n = A->size2;
  MPI_Barrier(MPI_COMM_WORLD);
  
  /*----------------------------------------
   * These are all variables related to GRock
   ----------------------------------------*/
  
  struct value table[size];
  gsl_vector *x        = gsl_vector_calloc(n);
  gsl_vector *As       = gsl_vector_calloc(n);
  gsl_vector *invAs    = gsl_vector_calloc(n);
  gsl_vector *local_b  = gsl_vector_calloc(m);
  gsl_vector *beta     = gsl_vector_calloc(n);
  gsl_vector *tmp      = gsl_vector_calloc(n);
  gsl_vector *d        = gsl_vector_calloc(n);
  gsl_vector *absd     = gsl_vector_calloc(n);
  gsl_vector *oldx     = gsl_vector_calloc(n);
  gsl_vector *tmpx     = gsl_vector_calloc(n);
  gsl_vector *z        = gsl_vector_calloc(m);
  gsl_vector *tmpz     = gsl_vector_calloc(m);
  gsl_vector *Ax       = gsl_vector_calloc(m);
  gsl_vector *Atmpx    = gsl_vector_calloc(m);
  gsl_vector *xdiff    = gsl_vector_calloc(n);
  gsl_permutation *idx = gsl_permutation_calloc(n);
  double send[1]; 
  double recv[1]; 
  double err;

  int num_upd = (int)(n*0.08);
  double sigma = 0.01;

  double xs_local_nrm[1], xs_nrm[1];
  double local_old_obj, global_old_obj, local_new_obj, global_new_obj;
  //calculate the 2 norm of xs
  xs_local_nrm[0] = gsl_blas_dnrm2(xs);
  xs_local_nrm[0] *=xs_local_nrm[0];
  MPI_Allreduce(xs_local_nrm, xs_nrm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  xs_nrm[0] = sqrt(xs_nrm[0]);
  
  // evaluate the two norm of the columns of A
  for(j=0;j<n;j++){
    gsl_vector_view column = gsl_matrix_column(A, j);
    double d;
    d = gsl_blas_dnrm2(&column.vector);
    gsl_vector_set(As, j, d*d);
    gsl_vector_set(invAs, j, 1./(d*d));
  }
  
  if (rank == 0) {
    printf("=============================\n");
    printf("|GRock start to solve Lasso!|\n");
    printf("|---------------------------|\n");
    printf("|lambda=%1.2f, m=%d, n=%d  |\n", lambda, m, n*size);
    if(mode==1) printf("| Mode: dynamic update P.   |\n");
    else  printf("|   Mode: fixed update P    |\n");
    printf("=============================\n");
    printf("%3s %8s %8s %5s\n", "iter", "rel_err", "obj", "P");
    startTime = MPI_Wtime();
    sprintf(s, "results/test%d.m", size);
    test = fopen(s, "w");
    fprintf(test,"res = [ \n");
  }
  
  /* Main BCD loop */
  total_start_time = MPI_Wtime();
  int iter = 0;
  while (iter < MAX_ITER) {
    startTime = MPI_Wtime();

    /*---------- restore the old x ------------*/
    gsl_vector_memcpy(oldx, x);
    
    /*------- calculate local_b = b - sum_{j \neq i} Aj*xj--------- */ 
    gsl_blas_dgemv(CblasNoTrans, 1, A, x, 0, Ax); // Ax = A * x
    MPI_Allreduce(Ax->data, z->data,  m, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    gsl_vector_sub(z, b); // z = Ax - b
    gsl_vector_memcpy(local_b, Ax);
    gsl_vector_sub(local_b, z);
    
    /* -------calculate beta ------------------*/
    gsl_blas_dgemv(CblasTrans, -1, A, z, 0, beta); // beta = A'(b - Ax) + ||A.s||^2 * xs
    gsl_vector_memcpy(tmp, As);    
    pointwise(tmp, x, n);
    gsl_vector_add(beta, tmp);
    shrink(beta, lambda);
    // x = 1/|xs|^2 * shrink(beta, lambda)
    gsl_vector_memcpy(x, beta);
    pointwise(x, invAs, n); 
  
    /* ------calcuate proposed decrease -------- */
    gsl_vector_memcpy(d,x);
    gsl_vector_sub(d, oldx);
    if(mode ==1){
      gsl_vector_memcpy(absd, d);
      abs_vector(absd, n);
      // sort the local array d
      gsl_vector_scale(absd, -1.0);
      gsl_sort_vector_index(idx, absd);

      //    printf("|d(0)| = %lf, |d(1)| = %lf \n", gsl_vector_get(absd,0), gsl_vector_get(absd, 3));
      // calculate current objective value;
      local_old_obj = objective(oldx, lambda, z, size);
      MPI_Allreduce(&local_old_obj, &global_old_obj, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      num_upd = fmin(num_upd+1, (int)(0.1*n));    
      gsl_vector_memcpy(tmpx, oldx);
      int upd_idx;
      double local_delta = 0, delta=0.0;
      for(int i=0; i<num_upd; i++){
	upd_idx = gsl_permutation_get(idx, i);
	//      printf("%d\n", upd_idx);
	gsl_vector_set(tmpx, upd_idx, gsl_vector_get(x, upd_idx));
	local_delta += gsl_vector_get(d, upd_idx) * gsl_vector_get(d, upd_idx);
      }
      MPI_Allreduce(&local_delta, &delta,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);    
      gsl_blas_dgemv(CblasNoTrans, 1, A, tmpx, 0, Atmpx); // Ax = A * x
      MPI_Allreduce(Atmpx->data, tmpz->data,  m, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      gsl_vector_sub(tmpz, b); // z = Ax - b
    
      local_new_obj = objective(tmpx, lambda, tmpz, size);
      MPI_Allreduce(&local_new_obj, &global_new_obj, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      while(global_new_obj - global_old_obj> -sigma * delta){
	num_upd = fmax(num_upd-1, 1);
	for(int i=0; i<num_upd; i++){
	  upd_idx = gsl_permutation_get(idx, i);
	  gsl_vector_set(tmpx, upd_idx, gsl_vector_get(x, upd_idx));
	  local_delta += gsl_vector_get(d, upd_idx) * gsl_vector_get(d, upd_idx);
	}
	MPI_Allreduce(&delta, &local_delta,  1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);    
	gsl_blas_dgemv(CblasNoTrans, 1, A, tmpx, 0, Atmpx); // Ax = A * x
	MPI_Allreduce(Atmpx->data, tmpz->data,  m, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	gsl_vector_sub(tmpz, b); // z = Ax - b
	
	local_new_obj = objective(tmpx, lambda, tmpz, size);
	MPI_Allreduce(&local_new_obj, &global_new_obj, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
	if(num_upd==1)
	  break;
      }

      gsl_vector_memcpy(x, tmpx);
    }  

    if(mode==0){
      CBLAS_INDEX_t id = gsl_blas_idamax(d);
      double *store = (double*)calloc(size, sizeof(double));
      double foo[1];
      foo[0] = gsl_vector_get(d,id);
      MPI_Allgather(foo, 1, MPI_DOUBLE, store, 1, MPI_DOUBLE, MPI_COMM_WORLD);
      for(int i=0;i<size;i++){
	table[i].ID   = i;
	table[i].data = fabs(store[i]);
      }
      // quick sort to decide which block to update
      qsort((void *) & table, size, sizeof(struct value), (compfn)compare );
      gsl_vector_memcpy(x, oldx);
      
      if(size>P){
	for(int i=0;i<P;i++){
	  if(rank == table[i].ID)
	    gsl_vector_set(x, id, gsl_vector_get(oldx, id) + gsl_vector_get(d, id));
	}
      }else
	gsl_vector_set(x, id, gsl_vector_get(oldx, id) + gsl_vector_get(d, id));
      local_new_obj = objective(x, lambda, z, size);
      MPI_Allreduce(&local_new_obj, &global_new_obj, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    
    /*------------------------------
      calculate the relative error
      ------------------------------*/
    gsl_vector_memcpy(xdiff,xs);
    gsl_vector_sub(xdiff, x);
    err = gsl_blas_dnrm2(xdiff);
    send[0] = err*err;
    MPI_Allreduce(send, recv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    recv[0] = sqrt(recv[0])/xs_nrm[0];
 
    endTime = MPI_Wtime();
    if(mode==1) P = num_upd*size;
    if (rank == 0) {
      if(iter%5 == 0)
	printf("%3d %10.2e %10.4f %3d\n", iter,
	       recv[0],  global_new_obj, P);
      fprintf(test, "%e \n",recv[0]);
    }

    /* termination check */
    if(recv[0] < TOL){
      break;
    }
    iter++;
  }
  total_end_time = MPI_Wtime();  
  /* Have the master write out the results to disk */
  if (rank == 0) {
    printf("=============================\n");
    printf("|    GRock solved Lasso!    |\n");
    printf("|---------------------------|\n");
    printf("|Summary:                   |\n");
    printf("|   # of iteration: %d      |\n", iter);
    printf("|   relative error: %4.2e|\n", recv[0]);
    printf("|  objective value: %4.2f    |\n", global_new_obj);
    printf("|             time: %4.1es|\n", total_end_time - total_start_time);
    printf("=============================\n");
    
    fprintf(test,"] \n");
    fprintf(test,"semilogy(1:length(res),res); \n");
    fprintf(test,"xlabel('# of iteration'); ylabel('||x - xs||');\n");
    fclose(test);
    f = fopen("results/solution.dat", "w");
    fprintf(f,"x = [ \n");
    gsl_vector_fprintf(f, x, "%lf");
    fprintf(f,"] \n");
    fclose(f);
    endTime = MPI_Wtime();
  }
  
  MPI_Finalize(); /* Shut down the MPI execution environment */
  
  /* Clear memory */
  gsl_matrix_free(A);
  gsl_vector_free(b);
  gsl_vector_free(x);
  gsl_vector_free(z);
  gsl_vector_free(xdiff);
  gsl_vector_free(Ax);
  gsl_vector_free(As);
  gsl_vector_free(invAs);
  gsl_vector_free(tmpx);
  gsl_vector_free(oldx);
  gsl_vector_free(local_b);
  gsl_vector_free(beta);
  gsl_vector_free(tmpz);
  gsl_vector_free(absd);
  gsl_vector_free(Atmpx);
  gsl_permutation_free(idx);

  return 0;
}


/* ----------- evaluate the objective function --------------*/
double objective(gsl_vector *x, double lambda, gsl_vector *z, int N) {
  double obj = 0;
  double temp =0.0;
  temp = gsl_blas_dnrm2(z);
  temp = temp*temp/(double)(2.0*N);
  double foo;
  foo = gsl_blas_dasum(x);
  //  double recv;
  //  MPI_Allreduce(&foo, &recv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  obj = lambda*foo + temp;
  return obj;
}

/*----------- shrinkage function ---------------------- */
void shrink(gsl_vector *v, double sigma) {
  double vi;
  for (int i = 0; i < v->size; i++) {
    vi = gsl_vector_get(v, i);
    if (vi > sigma)       { gsl_vector_set(v, i, vi-sigma); }
    else if (vi < -sigma) { gsl_vector_set(v, i, vi+sigma); }
    else              { gsl_vector_set(v, i, 0); }
  }
}


/* ----------- point wise product function ----------------- */
void pointwise(gsl_vector *a, gsl_vector *b, double n){
  for(int i=0; i<n; i++)
    gsl_vector_set(a, i, gsl_vector_get(a, i) * gsl_vector_get(b, i));
}

void abs_vector(gsl_vector *a, double n){
  for(int i=0; i<n; i++)
    gsl_vector_set(a, i, fabs(gsl_vector_get(a, i)));
}


int compare(struct value *tab1, struct value *tab2){

  if(tab1->data < tab2->data)
    return 1;
  else if (tab1->data > tab2->data)
    return -1;
  else
    return 0;
}
