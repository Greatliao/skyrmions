#ifndef TPEM_H
#define TPEM_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#ifndef NO_GSL
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#endif

#ifdef	__cplusplus
extern	"C" {
#endif




#define	COMPLEX	1


#ifdef	COMPLEX
typedef	struct	{ double real, imag; } tpem_t;
#else
typedef	double	tpem_t;
#endif	/* COMPLEX */


void	*tpem_malloc	(size_t);
void	tpem_print	(tpem_t a, FILE *b); 
tpem_t	tpem_cast	(double, double);
double	tpem_real	(tpem_t);
double tpem_imag (tpem_t z);
double tpem_chebyshev(int m,double x);


/* tpem_sparse.c
 * 
 * struct tpem_sparse implements a sparse matrix in compressed row storage
 * format.
 * The CRS format puts the subsequent nonzero elements of the matrix rows
 * in contiguous memory locations. We create 3 vectors: one for tpem_t numbers
 * (*values) and the other two for integers (*colind, *rowptr).
 *
 * The *values vector stores the values of the nonzero elements of the matrix,
 * as they are traversed in a row-wise fashion.
 *
 * The *colind vector stores the column indices of the elements of the *values
 * vector. That is, if
 *	values[k] = a[i][j] then colind[k] = j
 *
 * The *rowptr vector stores the locations in the *values vector that start
 * a row, that is
 *	values[k] = a[i][j] if rowptr[i] <= i < rowptr[i + 1]
 *
 * By convention, we define rowptr[n] = #(nz), the number of nonzeros in the
 * matrix. The storage savings of this approach is significant. Instead of
 * storing n^2 elements, we need only 2#(nz) + n + 1 storage locations.
 *
 * As an example, consider the nonsymmetric matrix defined by
 *	[10  0  0  0 -2  0 ]
 *	[ 3  9  0  0  0  3 ]
 * A =	[ 0  7  8  7  0  0 ]
 *	[ 3  0  8  7  5  0 ]
 *	[ 0  8  0  9  9 13 ]
 *	[ 0  4  0  0  2 -1 ]
 * The CRS format for this matrix is then specified by the arrays given below
 * values = [10 -2  3  9  3  7  8  7  3 ... 9 13  4  2 -1 ]
 * colind = [ 0  4  0  1  5  1  2  3  0 ... 4  5  1  4  5 ]
 * rowptr = [ 0  2  5  8 12 16 19 ]
 *
 * Here's an example illustrating matrix vector product:
 *
 * void tpem_sparse_mult (sparse *matrix, tpem_t *dest, tpem_t *src)
 * {
 *	size_t row, col, k;
 *	tpem_t tmp;
 *
 *	for (row = 0; row < matrix->rank; row++) {
 *		sum = 0.0;
 *		for (k = matrix->rowptr[row]; k < matrix->rowptr[row+1]; k++) {
 *			col  = matrix->colind[k];
 *			sum += matrix->values[k] * src[col];
 *		}
 *		dest[row] = sum;
 *	}
 * }
 *
 * new_tpem_sparse (rank, size) is a constructor
 *	rank:	number of rows (we don't consider degenerate matrices)
 *	size:	#(nz), the total number of nonzero elements
 *
 * tpem_sparse_element (this, row, col)
 *	this:	pointer to the matrix (think pathetic OO)
 *	row:	row index
 *	col:	column index
 *	returns a pointer to a[row][col] or freaks
 * 
 */
typedef	struct {
	size_t	rank;		/* number of rows */
	size_t	*rowptr;
	size_t	*colind;
	tpem_t	*values;
} tpem_sparse;


tpem_sparse	*new_tpem_sparse	(size_t, size_t);
void		tpem_sparse_free	(tpem_sparse *);
void		tpem_sparse_copy	(tpem_sparse *, tpem_sparse *);
tpem_t		*tpem_sparse_element	(tpem_sparse *, size_t, size_t);
/* void		tpem_sparse_print	(tpem_sparse *, FILE *); */

void		tpem_matrix_free	(tpem_t **);
tpem_t		**new_tpem_matrix	(size_t, size_t);


/* subspace.c
 *
 * struct tpem_subspace is opaque to the user, don't worry about it. But if
 * you must know, the boolean *flags vector tells whether the i^th state has
 * been used so far, think of the macro
 *	#define	isused(i)	(flags[i] != 0)
 * *stack and *top implement a stack of all states that have been used since
 * the last call to tpem_subspace_reset. The rest should be self-explanatory.
 *
 */
typedef	struct	{
	size_t	size;
	int	*flags;
	size_t	*stack, *top;
} tpem_subspace;


tpem_subspace	*new_tpem_subspace	(size_t);
void		tpem_subspace_free	(tpem_subspace *);
void		tpem_subspace_reset	(tpem_subspace *);
void		tpem_subspace_push	(tpem_subspace *, size_t);
/* void		tpem_subspace_print	(tpem_subspace *, FILE *);
*/

/* tpem.c */
typedef	double	(*tpem_func_ptr)	(double);
void	tpem_off_diagonal_element_tpem	(tpem_sparse *, size_t, tpem_t **,
					size_t, double);
void	tpem_calculate_moment_tpem		(tpem_sparse *, size_t, double *,
					double);
void	tpem_calculate_moment_pem		(tpem_sparse *, size_t, double *);
void	tpem_calculate_moment_diff_tpem	(tpem_sparse *, tpem_sparse *,
					size_t, double *,
					size_t, size_t *, double, double);
void	tpem_calculate_moment_diff_pem	(tpem_sparse *, tpem_sparse *,
					size_t, double *);
void	tpem_calculate_coeffs	(size_t, double *,tpem_func_ptr f);
void	tpem_calculate_coeffs_alt	(size_t, double *,tpem_func_ptr f);
void 	tpem_calculate_coeffs_dos 	(size_t , double *, size_t , size_t , size_t, size_t );
double	tpem_expansion			(size_t, double *, double *);
void	tpem_finalize			(void);
void 	tpem_init			(void);

/*
 * 	CALL TREE:
 * 	----------
 *
 *	tpem_calculate_moment:
 *		\tpem_diagonal_element:
 *			\tpem_sparse_product
 *
 *	tpem_calculate_moment_diff:
 *		|tpem_subspace_for_trace:
 *		|	\tpem_diagonal_element:
 *		|		\tpem_sparse_product
 *		\tpem_diagonal_element:
 *			\tpem_sparse_product
 *
 *	tpem_finalize:
 *		|tpem_diagonal_element
 *		|tpem_calculate_moment
 *		\tpem_calculate_moment_diff
 */


#ifdef	__cplusplus
}
#endif


#endif
