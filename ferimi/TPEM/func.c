#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "tpem.h"

#ifdef NO_LAPACK
#include "diag.h"
#endif


/* Global Variables */
double beta, miu, a, b;

double orig_func( double x)
{
    
    //return log( 1+exp( -beta*(a*x+b-miu) ) );
    return 1;
}

double kernel(size_t cutoff, size_t i)
{
	double y;

	y = 1.0/(cutoff+1.0)*( (cutoff-i+1.0)*cos(M_PI*i/(cutoff+1.0))+sin(M_PI*i/(cutoff+1.0))*tan(M_PI/2-M_PI/(cutoff+1.0)) );
	return y;
}

void tpem_statedensity_( int *num_x1, int *cut, int *rank, int *nonzeros, double *val1, double *val2, int *col, int *rowp, double *stateden)
{
    size_t num_x = *num_x1;
    size_t i, j, cutoff= *cut;

    /* para in this func*/
	double eps = 0.00001;
	double *moment = tpem_malloc (cutoff * sizeof *moment);
	/* para in this func*/

    tpem_sparse *matrix = new_tpem_sparse(*rank, *nonzeros);
    for (i = 0; i < *nonzeros; i++){
        matrix->colind[i] = col[i];
        matrix->values[i] = tpem_cast(val1[i], val2[i]);
    }
    for (i = 0; i <= *rank; i++){
        matrix->rowptr[i] = rowp[i];
    }
    tpem_init();

	tpem_calculate_moment_tpem(matrix, cutoff, moment, eps);
	
	for (j = 0; j < cutoff; j++){
	    printf("%f",moment[j]);
	    printf("\n");
	}
	
	for (i = 0; i < num_x; i++){ 
		stateden[i] = kernel(cutoff,0)*moment[0]/M_PI/sqrt(1-(-1+2*(i+0.5)/num_x)*(-1+2*(i+0.5)/num_x));
		for (j = 1; j < cutoff; j++){ 
			stateden[i] = stateden[i]+2.0*kernel(cutoff,j)*moment[j]*tpem_chebyshev(j,(-1+2*(i+0.5)/num_x))/M_PI/sqrt(1-(-1+2*(i+0.5)/num_x)*(-1+2*(i+0.5)/num_x));
		}
		/*
	    printf("%f",stateden[i]);
	    printf("\n");
		*/
	}

	free(moment);
    tpem_sparse_free(matrix);
    tpem_finalize();

}

void tpem_logweight_( int *cut, int *rank, int *nonzeros, double *val1, double *val2, int *col, int *rowp, double *beta1, double *miu1, double *a1, double *b1, double *weight )
{
    beta = *beta1;
    miu = *miu1;
    a = *a1;
    b = *b1;
    size_t i, cutoff=*cut;

    /* para in this func*/
	double	eps = 0.00001;
	double	*coeffs = tpem_malloc (cutoff * sizeof *coeffs);
	double	*moment	= tpem_malloc (cutoff * sizeof *moment);
    /* para in this func*/

    tpem_sparse *matrix = new_tpem_sparse(*rank, *nonzeros);
    for (i = 0; i < *nonzeros; i++){
        matrix->colind[i] = col[i];
        matrix->values[i] = tpem_cast(val1[i], val2[i]);
    }
    for (i = 0; i <= *rank; i++){
        matrix->rowptr[i] = rowp[i];
    }
    tpem_init();
	
#ifndef NO_GSL	
	tpem_calculate_coeffs (cutoff, coeffs, orig_func);
#else
	tpem_calculate_coeffs_alt (cutoff, coeffs, orig_func);
#endif
	tpem_calculate_moment_tpem(matrix, cutoff, moment,eps);
	for (i = 1; i < cutoff; i++){ 
		moment[i] = moment[i]*kernel(cutoff,i);
	}
	*weight = tpem_expansion (cutoff, moment, coeffs);
	free (coeffs);
	free (moment);

    tpem_sparse_free(matrix);
    tpem_finalize();
}
