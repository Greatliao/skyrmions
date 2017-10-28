#include <stdio.h>
#include <stdlib.h>
#include "tpem.h"


void tpem_subspace_free (tpem_subspace *this)
{
	if (this && this->stack) free (this->stack), this->stack = NULL;
	if (this && this->flags) free (this->flags), this->flags = NULL;
	if (this) free (this), this = NULL;
}


void tpem_subspace_reset (tpem_subspace *this)
{
	size_t i;

	for (i = 0; i < this->size; i++)
		this->flags[i] = 0;
	this->top = this->stack;
}


tpem_subspace *new_tpem_subspace (size_t size)
{
	tpem_subspace *this = tpem_malloc (sizeof *this);

	this->size = size;
	this->flags = tpem_malloc (size * sizeof *this->flags);
	this->stack = tpem_malloc (size * sizeof *this->stack);
	tpem_subspace_reset (this);
	return this;
}


void tpem_subspace_push (tpem_subspace *this, size_t state)
{
	if (this->flags[state] == 0) {
		this->flags[state] = 1;
		*(this->top)++ = state;
	}
}


void tpem_subspace_print (tpem_subspace *this, FILE *fp)
{
	size_t *p;

	for (p = this->stack; p < this->top; p++)
		fprintf (fp, "state[%u] = %u\n", p - this->stack, *p);
	fprintf (fp, "tpem_subspace_print done\n");
}
