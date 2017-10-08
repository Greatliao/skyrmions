

#include "tpem.h"


void *tpem_malloc (size_t n)
{	/* malloc and report if error */
	void	*p = malloc (n);

	if (p == NULL) {	/* print error message and exit */
		fprintf (stderr, "malloc of %u bytes failed: ", n);
		fprintf (stderr, "%s\n", strerror (errno));
		exit (2);	/* conventional value for failed execution */
	}
	return p;
}


void tpem_sparse_free (tpem_sparse *this)
{
	if (this && this->rowptr) free (this->rowptr), this->rowptr = NULL;
	if (this && this->colind) free (this->colind), this->colind = NULL;
	if (this && this->values) free (this->values), this->values = NULL;
	if (this) free (this), this = NULL;
}


tpem_sparse *new_tpem_sparse (size_t rank, size_t size)
{
	tpem_sparse *this = tpem_malloc (sizeof *this);

	this->rank = rank;
	this->rowptr = tpem_malloc ((rank + 1) * sizeof *this->rowptr);
	this->colind = tpem_malloc (size * sizeof *this->colind);
	this->values = tpem_malloc (size * sizeof *this->values);
	this->rowptr[rank] = size;
	return this;
}


void tpem_sparse_copy (tpem_sparse *this, tpem_sparse *other)
{
	size_t i, n;

	n = other->rank = this->rank;
	for (i = 0; i < n; i++)
		other->rowptr[i] = this->rowptr[i];
	n = other->rowptr[n] = this->rowptr[n];
	for (i = 0; i < n; i++) {
		other->colind[i] = this->colind[i];
		other->values[i] = this->values[i];
	}
}


tpem_t *tpem_sparse_element (tpem_sparse *this, size_t row, size_t col)
{
	size_t i;

	for (i = this->rowptr[row]; i < this->rowptr[row + 1]; i++)
		if (this->colind[i] == col)
			return this->values + i;
	fprintf (stderr, "tpem_sparse_element: no element at [%d][%d]\n",
			row, col);
	fprintf (stderr, "tpem_sparse_element: debug: rowprtr[row]=%d, rowptr[row+1]=%d\n",
			this->rowptr[row], this->rowptr[row+1]);
	return NULL;
}


void tpem_sparse_print (tpem_sparse *this, FILE *fp)
{
	size_t i, k;

	for (i = 0; i < this->rank; i++) {
		for (k = this->rowptr[i]; k < this->rowptr[i + 1]; k++) {
			fprintf (fp, "[%u][%u] = ", i, this->colind[k]);
			tpem_print (this->values[k], fp);
			fputc ('\t', fp);
		}
		fputc ('\n', fp);
	}
}


void tpem_matrix_free (tpem_t **this)
{
	if (this && *this) free (*this), *this = NULL;
	if (this) free (this), this = NULL;
}


tpem_t **new_tpem_matrix (size_t row, size_t col)
{
	size_t i;
	tpem_t **this = tpem_malloc (row * sizeof *this);

	*this = tpem_malloc (row * col * sizeof **this);
	for (i = 1; i < row; i++)
		this[i] = this[i - 1] + col;
	return this;
}
