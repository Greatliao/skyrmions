#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "tpem.h"

#ifdef NO_LAPACK
#include "diag.h"
#endif

/* Random number generator */
double my_random (void)
{
	static int next = 1;
	
	next = 16807 * (next % 127773) - 2836 * (next / 127773);
	if (next <= 0) next += 2147483647;
	return ((double) next) / 2147483647.0;
}


/* Initialize the matrix using nearest neighbor hoppings and random potentials */
tpem_sparse *new_tpem_sparse_random (size_t rank, double range)
{
	tpem_sparse *this = new_tpem_sparse (rank, 3 * rank);
	double scale = 2.0 + fabs (range), epsilon;
	tpem_t	hopping = tpem_cast (-1.0 / scale, 0.0);
	size_t	i, j;

	range *= 2.0 / scale;
	for (i = j = 0; i < rank; i++) {
		this->rowptr[i] = j;
		this->colind[j] = i;
		epsilon = range * (my_random () - 0.5);
		this->values[j] = tpem_cast (epsilon, 0.0);
		j++;
		this->colind[j] = (i + 1) % rank;
		this->values[j] = hopping;
		j++;
		this->colind[j] = (i + rank - 1) % rank;
		this->values[j] = hopping;
		j++;
	}
	return this;
}

/* Calculates an observable defined by the function funcptr using the TPEM for
   the model given by the Hamiltonian matrix *matrix */
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

#ifndef NO_LAPACK
/* Calculates an observable defined by the function funcptr using exact diagonalization
    for the model given by the Hamiltonian matrix *matrix */
static double diag_apply (tpem_sparse *matrix, double (*funcptr)(double))
{
	char	jobz	= 'n';
	char	uplo	= 'u';
	int	n	= matrix->rank, info, i;
	size_t	k;
#ifdef	COMPLEX
	int	lwork 	= 2 * n - 1;
	double	*rwork	= tpem_malloc ((3 * n - 2) * sizeof *rwork);
	extern	void	zheev_ (char *, char *, int *, tpem_t *, int *,
			double *, tpem_t *, int *, double *, int *);
#else
	int	lwork	= 3 * n - 1;
	extern	void	dsyev_ (char *, char *, int *, tpem_t *, int *,
			double *, tpem_t *, int *, int *);
#endif
	tpem_t	*a	= tpem_malloc (n * n * sizeof *a);
	double	*w	= tpem_malloc (n * sizeof *w), ret;
	tpem_t	*work	= tpem_malloc (lwork * sizeof *work);

	for (i = 0; i < n * n; i++)
		a[i] = tpem_cast (0.0, 0.0);
	for (i = 0; i < n; i++)
		for (k = matrix->rowptr[i]; k < matrix->rowptr[i + 1]; k++)
			a[i * n + matrix->colind[k]] = matrix->values[k];
#ifdef	COMPLEX
	zheev_ (&jobz, &uplo, &n, a, &n, w, work, &lwork, rwork, &info);
	free (rwork);
#else
	dsyev_ (&jobz, &uplo, &n, a, &n, w, work, &lwork, &info);
#endif
	ret = 0.0;
	for (i = 0; i < n; i++) {
		ret += funcptr (w[i]);
		/* fprintf(stderr,"Eigenvalues %f\n",w[i]); */
	}	
	free (a);
	free (w);
	free (work);
	return ret;
}
#else
static double diag_apply (tpem_sparse *matrix, double (*funcptr)(double))
{
	
	int	n	= matrix->rank, i;
	size_t	k;
	ccomplex *a;
	double ret;
	double *w;
	
	a = malloc(n*n*sizeof(ccomplex));
	w = malloc(n*sizeof(double));
	for (i = 0; i < n * n; i++)
		a[i].re = a[i].im = 0.0;
	for (i = 0; i < n; i++) {
		for (k = matrix->rowptr[i]; k < matrix->rowptr[i + 1]; k++) {
			a[i * n + matrix->colind[k]].re = tpem_real(matrix->values[k]);
			a[i * n + matrix->colind[k]].im = tpem_imag(matrix->values[k]);
		}
	}
	diag(n, a, w);
	
	ret = 0.0;
	for (i = 0; i < n; i++) {
		ret += funcptr (w[i]);
		/* fprintf(stderr,"Eigenvalues %f\n",w[i]); */
	}	
	free(a);
	free(w);
	
	return ret;
}
#endif

/* Calculates the difference in the observable defined by the function funcptr 
     using exact diagonalization
    for the matrices matrix0 and matrix1 */
static double tpem_apply_diff (tpem_sparse *matrix0, tpem_sparse *matrix1,
		tpem_func_ptr funcptr, size_t cutoff)
{
	double  eps_prod = 0.00001, eps_trace = 0.001;
	double  *coeffs = tpem_malloc (cutoff * sizeof *coeffs);
	double  *moment = tpem_malloc (cutoff * sizeof *moment);
	size_t	nsupp = 2, support[2];
	double  ret;

	support[0] = 0;
	support[1] = matrix0->rank / 2 - 1;
#ifndef NO_GSL	
	tpem_calculate_coeffs (cutoff, coeffs, funcptr);
#else
	tpem_calculate_coeffs_alt (cutoff, coeffs, funcptr);
#endif
	tpem_calculate_moment_diff_tpem (matrix0, matrix1, cutoff, moment,
			nsupp, support, eps_trace, eps_prod); 
	/* tpem_calculate_moment_diff_tpem (matrix0, matrix1, cutoff, moment); */
	ret = tpem_expansion (cutoff, moment, coeffs);
	free (coeffs);
	free (moment);
	return ret;
}


/* Emulates an energy function */
double orig_func1(double x)
{
	return 5.0 * x *  (1.0 - tanh (10.0 * x));
}

/* Emulates a density function */
double orig_func2(double x)
{
	return 0.5 * (1.0 - tanh (10.0 * x));
}

/* Emulates an action function */
double orig_func3(double x)
{
	return log (1.0 + exp (-10.0 * x));
}


int main (void)
{
	tpem_sparse *matrix0 = new_tpem_sparse_random (400, 10.0);
	tpem_sparse *matrix1 = new_tpem_sparse_random (400, 10.0);
	double	tpem0, tpem1, tpemd, naive0, naive1;
	size_t	cutoff;
	
	tpem_init();
	
	puts ("********************************************************");
	puts ("****** TESTING TRUNCATED POLYNOMIAL EXPANSION **********");
	puts ("********************************************************");
	puts ("");
	puts ("");
	puts ("This testing program calculates model properties in two ways:");
	puts ("(i)  Using standard diagonalization");
	puts ("(ii) Using the truncated polynomial expansion method");
	puts ("");
	puts ("All tests are done for a nearest neighbour interaction with");
	puts ("random (diagonal) potentials.");
	puts ("");
	puts ("");
	
	tpem_sparse_copy (matrix0, matrix1);
	matrix1->values[matrix1->rowptr[0]] =
		tpem_cast (2.4 * (my_random () - 0.5), 0.0);
	matrix1->values[matrix1->rowptr[matrix1->rank / 2 - 1]] =
		tpem_cast (2.4 * (my_random () - 0.5), 0.0);
	
	puts ("-------------------------------------------------------------");
	puts ("TEST 1: MEAN VALUE FOR THE FUNCTION:                         ");
	puts ("        N(x) =  0.5 * (1.0 - tanh (10.0 * x))");
	puts ("");
	naive0 = diag_apply (matrix0, orig_func2);
	printf ("** Using diagonalization <N>=%f\n", naive0);
	puts ("** Using TPEM <N>=(cutoff--> infinity) "
			"lim<N_cutoff> where <N_cutoff> is");
	puts ("cutoff\t<N_cutoff>\t\%Error(compared to diag.)");
	for (cutoff = 10; cutoff <= 40; cutoff++) {
		tpem0 = tpem_apply (matrix0, orig_func2, cutoff);
		printf ("%d\t%.12f\t%f%%\n", cutoff, tpem0,
				100.0 * fabs (1.0 - tpem0 / naive0));
		fflush(stdout);
	}
	
	puts ("-------------------------------------------------------------");
	puts ("TEST 2: MEAN VALUE FOR THE FUNCTION:                         ");
	puts ("        E(x) = 5.0 * x * (1.0 - tanh (10.0 * x))");
	puts ("");
	naive0 = diag_apply (matrix0, orig_func1);
	printf ("** Using diagonalization <E>=%f\n", naive0);
	puts ("** Using TPEM <E>=(cutoff--> infinity) "
			"lim<E_cutoff> where <E_cutoff> is");
	puts ("cutoff\t<E_cutoff>\t\%Error(compared to diag.)");
	for (cutoff = 10; cutoff <= 40; cutoff++) {
		tpem0 = tpem_apply (matrix0, orig_func1, cutoff);
		printf ("%d\t%.12f\t%f\t%f\n", cutoff, tpem0,naive0,
				100.0 * fabs (1.0 - tpem0 / naive0));
		fflush(stdout);
	}
	

	
	puts ("-------------------------------------------------------------");
	puts ("TEST 3: MEAN VALUE AND DIFFERENCE FOR THE FUNCTION:          ");
	puts ("        S(x) =  log (1.0 + exp (-20.0 * x))");
	puts ("");
	naive0 = diag_apply (matrix0, orig_func3);
	naive1 = diag_apply (matrix1, orig_func3);
	printf ("** Using diagonalization <S[matrix0]>= %.12f\n", naive0);
	printf ("** Using diagonalization <S[matrix1]>= %.12f\n", naive1);
	printf ("** Using diagonalization <S[matrix1]>-<S[matrix0]>=%.12f\n",
			naive0 - naive1);
	puts ("** Using TPEM <S>=(cutoff--> infinity) lim<S_cutoff>");
	puts ("cutoff\tDelta_S_cutoff\tS_cutoff[diff]\tError (to diag.)");
	for (cutoff = 10; cutoff <= 40; cutoff++) {
		tpem0 = tpem_apply (matrix0, orig_func3, cutoff);
		tpem1 = tpem_apply (matrix1, orig_func3, cutoff);
		tpemd =0.0;
		tpemd = tpem_apply_diff (matrix0, matrix1, orig_func3, cutoff);
		printf ("%d\t%.12f\t%.12f\t%f%%\n", cutoff,
				tpem0 - tpem1, tpemd,
				100.0 * fabs (1.0 - tpemd / (naive0-naive1)));
		fflush(stdout);
	}
	
	puts ("-------------------------------------------------------------");
   	
	tpem_sparse_free (matrix0);
	tpem_sparse_free (matrix1);
	tpem_finalize();
	return 0;
}
