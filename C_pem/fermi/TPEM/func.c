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
    
    return log( 1+exp( -beta*(a*x+b-miu) ) );
	
}

/*****************************************************************************************************************************/

double orig_func_particle_density( double x)
{
    
    return 1.0 / ( 1.0+exp( beta*(a*x+b-miu) ) );
	
}

/*****************************************************************************************************************************/

double kernel(size_t cutoff, size_t i)
{
	double y;

	y = 1.0/(cutoff+1.0)*( (cutoff-i+1.0)*cos(M_PI*i/(cutoff+1.0))+sin(M_PI*i/(cutoff+1.0))*tan(M_PI/2-M_PI/(cutoff+1.0)) );
	return y;
}

/*##############################################################################################################################*/

void pem_moments_( int *cut, int *rank, int *nonzeros, double *val1, double *val2, int *col, int *rowp, double *moments)
{
    size_t i, j, cutoff= *cut;

    tpem_sparse *matrix = new_tpem_sparse(*rank, *nonzeros);
    for (i = 0; i < *nonzeros; i++){
        matrix->colind[i] = col[i];
        matrix->values[i] = tpem_cast(val1[i], val2[i]);
    }
    for (i = 0; i <= *rank; i++){
        matrix->rowptr[i] = rowp[i];
    }

    tpem_init();
	tpem_calculate_moment_pem(matrix, cutoff, moments);

    tpem_sparse_free(matrix);
    tpem_finalize();

}


/*****************************************************************************************************************************/

void pem_statedensity_( int *num_x1, int *cut, int *rank, int *nonzeros, double *val1, double *val2, int *col, int *rowp, double *stateden)
{
    size_t i, j, cutoff= *cut, num_x = *num_x1;
	double x;
	double *moment = tpem_malloc (cutoff * sizeof *moment);

    tpem_sparse *matrix = new_tpem_sparse(*rank, *nonzeros);
    for (i = 0; i < *nonzeros; i++){
        matrix->colind[i] = col[i];
        matrix->values[i] = tpem_cast(val1[i], val2[i]);
    }
    for (i = 0; i <= *rank; i++){
        matrix->rowptr[i] = rowp[i];
    }

    tpem_init();
	tpem_calculate_moment_pem(matrix, cutoff, moment);

    for (j = 0; j < num_x; j++){

		x = -1+2.0 * j/num_x;

	    stateden[j] = kernel(cutoff,0)*moment[0];
	    for (i = 1; i < cutoff; i++){ 

		    stateden[j] += 2.0*tpem_chebyshev(i,x)*kernel(cutoff,i)*moment[i];

	    }
        stateden[j] /=  *rank/M_PI/sqrt(1-x*x);
    }
	free (moment);

    tpem_sparse_free(matrix);
    tpem_finalize();

}

/*****************************************************************************************************************************/

void pem_statedensity_integral_( int *num_x1, int *interval1, int *cut, int *rank, int *nonzeros, double *val1, double *val2, int *col, int *rowp, double *stateden)
{
    size_t i, j, cutoff= *cut, num_x = *num_x1, interval = *interval1;
	double x;
	double *moment = tpem_malloc (cutoff * sizeof *moment);
	double	*coeffs = tpem_malloc (cutoff * sizeof *coeffs);
	/*
	FILE *fp;
	fp = fopen("coeffs.out","w");
	*/
    tpem_sparse *matrix = new_tpem_sparse(*rank, *nonzeros);
    for (i = 0; i < *nonzeros; i++){
        matrix->colind[i] = col[i];
        matrix->values[i] = tpem_cast(val1[i], val2[i]);
    }
    for (i = 0; i <= *rank; i++){
        matrix->rowptr[i] = rowp[i];
    }

    tpem_init();
	tpem_calculate_moment_pem(matrix, cutoff, moment);

	for (i = 0; i < cutoff; i++){ 
		moment[i] = moment[i]*kernel(cutoff,i)/ *rank;
		/*
		printf("%f",coeffs[i]);
	   	 printf("\n");
		*/
	}

    for (j = 0; j < num_x; j++){
		tpem_calculate_coeffs_dos (cutoff, coeffs, num_x, interval, j, j+1);
		/*
		fprintf(fp,"%zu\n",j);	
		for (i = 0; i < cutoff; i++){
			if (fp != NULL)
				fprintf(fp,"%f\n",coeffs[i]);		
		} 
		*/
		stateden[j] = tpem_expansion (cutoff, moment, coeffs);
	}
	free (coeffs);
	free (moment);

    tpem_sparse_free(matrix);
    tpem_finalize();

}


/*****************************************************************************************************************************/

void pem_test_( int *cut, int *rank, int *nonzeros, double *val1, double *val2, int *col, int *rowp, double *miun_tpem, double *miun_pem, double *miun_ev, double *eigv)
{
    size_t i, j, cutoff= *cut;

    /* para in this func*/
	double	eps = 0.00001;
	double *moment_tpem = tpem_malloc (cutoff * sizeof *moment_tpem);
	double *moment_pem = tpem_malloc (cutoff * sizeof *moment_pem);
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
	tpem_calculate_moment_tpem(matrix, cutoff, moment_tpem, eps);
	tpem_calculate_moment_pem(matrix, cutoff, moment_pem);
	for (i = 0; i < cutoff; i++){ 

	    miun_tpem[i] = moment_tpem[i];
		miun_pem[i] = moment_pem[i];
		miun_ev[i] = 0.0;
		for (j = 0; j < *rank; j++){
			miun_ev[i] += tpem_chebyshev(i,eigv[j]);
	    }
	}
    
	free (moment_tpem);
	free (moment_pem);
    tpem_sparse_free(matrix);
    tpem_finalize();
}

/*****************************************************************************************************************************/

void pem_ev_dos_( int *cut, int *num_x1, int *rank, double *eigv, double *stateden)
{
    size_t i, j, cutoff= *cut, num_x = *num_x1;
	double *moment = tpem_malloc (cutoff * sizeof *moment);
	double x;

	for (i = 0; i < cutoff; i++){ 
		moment[i] = 0.0;
		for (j = 0; j < *rank; j++){
			moment[i] += tpem_chebyshev(i,eigv[j]);
		}
	}

    for (j = 0; j <= num_x; j++){
		x = -1+2.0 * j/num_x;

	    stateden[j] = kernel(cutoff,0)*moment[0];
	    for (i = 1; i < cutoff; i++){ 

		    stateden[j] += 2.0*tpem_chebyshev(i,x)*kernel(cutoff,i)*moment[i];
   
	    }
        stateden[j] /=  *rank/M_PI/sqrt(1-x*x);
	}
	free (moment);
}

/*****************************************************************************************************************************/

void pem_logweight_( int *cut, int *rank, int *nonzeros, double *val1, double *val2, int *col, int *rowp, double *beta1, double *miu1, double *a1, double *b1, double *weight )
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
	//tpem_calculate_moment_tpem(matrix, cutoff, moment,eps);
	tpem_calculate_moment_pem(matrix, cutoff, moment);
	for (i = 0; i < cutoff; i++){ 

		moment[i] = moment[i]*kernel(cutoff,i);
		
		//printf("%f\n",coeffs[i]);
		
	}
	*weight = tpem_expansion (cutoff, moment, coeffs);
	free (coeffs);
	free (moment);

    tpem_sparse_free(matrix);
    tpem_finalize();
}

/*****************************************************************************************************************************/

void pem_particle_density_( int *cut, int *rank, int *nonzeros, double *val1, double *val2, int *col, int *rowp, double *beta1, double *miu1, double *a1, double *b1, double *ParticleDensity )
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
	tpem_calculate_coeffs (cutoff, coeffs, orig_func_particle_density);
#else
	tpem_calculate_coeffs_alt (cutoff, coeffs, orig_func_particle_density);
#endif
	//tpem_calculate_moment_tpem(matrix, cutoff, moment,eps);
	tpem_calculate_moment_pem(matrix, cutoff, moment);
			
	for (i = 0; i < cutoff; i++){ 

		moment[i] = moment[i]*kernel(cutoff,i)/ *rank;

	}
	*ParticleDensity = tpem_expansion (cutoff, moment, coeffs);
	free (coeffs);
	free (moment);

    tpem_sparse_free(matrix);
    tpem_finalize();
}
