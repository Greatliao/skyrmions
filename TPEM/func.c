#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "tpem.h"

#ifdef NO_LAPACK
#include "diag.h"
#endif


/* Global Variables */
double beta1, miu1, a1, b1;

double orig_func( double x)
{
    
    return log( 1+exp( -beta1*(a1*x+b1-miu1) ) );
}

static double tpem_apply (tpem_sparse *matrix, tpem_func_ptr funcptr,
		size_t cutoff)
{
	double	eps = 0.00001;
	double	*coeffs = tpem_malloc (cutoff * sizeof *coeffs);
	double	*moment	= tpem_malloc (cutoff * sizeof *moment);
	double	ret;
	
#ifndef NO_GSL	
	tpem_calculate_coeffs (cutoff, coeffs, funcptr);
#else
	tpem_calculate_coeffs_alt (cutoff, coeffs, funcptr);
#endif
	tpem_calculate_moment_tpem(matrix, cutoff, moment,eps);
	ret = tpem_expansion (cutoff, moment, coeffs);
	free (coeffs);
	free (moment);
	return ret;
}

double ctofortran_( int *rank, int *nonzeros, double *val1, double *val2, int *col, int *rowp, double *beta, double *miu, double *a, double *b )
{
    beta1 = *beta;
    miu1 = *miu;
    a1 = *a;
    b1 = *b;
    size_t i, cutoff=30;
    double ln_weight;

    tpem_sparse *matrix = new_tpem_sparse(*rank, *nonzeros);
    for (i = 0; i < *nonzeros; i++){
        matrix->colind[i] = col[i];
        matrix->values[i] = tpem_cast(val1[i], val2[i]);
    }
    for (i = 0; i <= *rank; i++){
        matrix->rowptr[i] = rowp[i];
    }
    tpem_init();
    ln_weight = tpem_apply (matrix, orig_func, cutoff);
    tpem_sparse_free(matrix);
    tpem_finalize();
    return ln_weight;

}


