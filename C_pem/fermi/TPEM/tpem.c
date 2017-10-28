#include "tpem.h"

struct my_f_params { int m;  tpem_func_ptr funk;};

void my_handler (const char * reason, const char * file, int line, int gsl_errno);

void tpem_print (tpem_t z, FILE *fp)
{
#ifdef	COMPLEX
	fprintf (fp, "(%f,%f)", z.real, z.imag);
#else
	fprintf (fp, "%f", z);
#endif
}


tpem_t tpem_cast (double real, double imag)
{
#ifdef	COMPLEX
	tpem_t	z;

	z.real = real;
	z.imag = imag;
	return	z;
#else
	return real;
#endif
}


double tpem_real (tpem_t z)
{
#ifdef	COMPLEX
	return	z.real;
#else
	return	z;
#endif
}

double tpem_imag (tpem_t z)
{
#ifdef	COMPLEX
	return	z.imag;
#else
	return	z;
#endif
}

tpem_t tpem_conj (tpem_t z)
{
#ifdef	COMPLEX
	tpem_t	w;

	w.real = z.real;
	w.imag = -z.imag;
	return	w;
#else
	return	z;
#endif
}

static tpem_t tpem_mult (tpem_t x, tpem_t y)
{
#ifdef	COMPLEX
	tpem_t	z;

	z.real = x.real * y.real - x.imag * y.imag;
	z.imag = x.real * y.imag + x.imag * y.real;
	return	z;
#else
	return	x * y;
#endif
}

tpem_t tpem_conj_mult (tpem_t x, tpem_t y)
{
	return tpem_mult (tpem_conj (x),y);
}

static double tpem_norm (tpem_t z)
{
#ifdef	COMPLEX
	return	z.real * z.real + z.imag * z.imag;
#else
	return	z * z;
#endif
}


static tpem_t tpem_add (tpem_t x, tpem_t y)
{
#ifdef	COMPLEX
	tpem_t	z;

	z.real = x.real + y.real;
	z.imag = x.imag + y.imag;
	return	z;
#else
	return	x + y;
#endif
}


static tpem_t tpem_sub (tpem_t x, tpem_t y)
{
#ifdef	COMPLEX
	tpem_t	z;

	z.real = x.real - y.real;
	z.imag = x.imag - y.imag;
	return	z;
#else
	return	x - y;
#endif
}
	

static void tpem_sparse_product_tpem (tpem_sparse *matrix, tpem_subspace *info,
		tpem_t *dest, tpem_t *src, double eps)
{	
	size_t	*oldtop, *p, i, j, k;
	tpem_t	t;
	double	u;

	for (i = 0; i < matrix->rank; i++)
		dest[i] = tpem_cast (0.0, 0.0);

	oldtop = info->top;
	/* loop over states that have been used so far */
	for (p = info->stack; p < oldtop; p++) {
		j = *p;
		/* loop over nonzero elements of j^th row */
		for (k = matrix->rowptr[j]; k < matrix->rowptr[j + 1]; k++) {
			i = matrix->colind[k];
			t = tpem_mult (src[j], matrix->values[k]);
			u = tpem_norm (t);
			if (info->flags[i]) {
				dest[i] = tpem_add (dest[i], t);
			} else if (u > eps) {
				dest[i] = t;
				tpem_subspace_push (info, i);
			}
		}
	}
}

static void tpem_sparse_product_pem (tpem_sparse *matrix,tpem_t *dest, tpem_t *src)
{	
	size_t	i, j, k;
	tpem_t	t;
	
	for (i = 0; i < matrix->rank; i++)
		dest[i] = tpem_cast (0.0, 0.0);

	
	/* loop over all rows */
	for (j=0;j<matrix->rank;j++) {
		/* loop over nonzero elements of j^th row */	
		for (k = matrix->rowptr[j]; k < matrix->rowptr[j + 1]; k++) {
			i = matrix->colind[k];
			t = tpem_mult (src[i], matrix->values[k]);
			dest[j] = tpem_add (dest[j], t);
		}
	}
}


static tpem_subspace *tpem_diagonal_element_tpem (tpem_sparse *matrix, size_t n,
		double *moment, size_t ket, double eps)
{	
	static	int	firstcall = 1;
	static	tpem_t	*tmp, *jm0, *jm1;
	static	tpem_subspace	*info;
	size_t	i, m, *p;
	tpem_t	sum1, sum2, keep;
	
	if (matrix == NULL || moment == NULL) {
		if (tmp) free (tmp), tmp = NULL;
		if (jm0) free (jm0), jm0 = NULL;
		if (jm1) free (jm1), jm1 = NULL;
		tpem_subspace_free (info);
		firstcall = 1;
		return NULL;
	}
	if (firstcall) {
		firstcall = 0;
		tmp = tpem_malloc (matrix->rank * sizeof *tmp);
		jm0 = tpem_malloc (matrix->rank * sizeof *jm0);
		jm1 = tpem_malloc (matrix->rank * sizeof *jm1);
		info = new_tpem_subspace (matrix->rank);
	} else if (info->size != matrix->rank) {
		tpem_finalize ();
		tpem_diagonal_element_tpem (matrix, n, moment, ket, eps);
	}
	
	for (i = 0; i < matrix->rank; i++)
		tmp[i] = jm0[i] = jm1[i] = tpem_cast (0.0, 0.0);
	jm0[ket] = tpem_cast (1.0, 0.0);	/* set |j,0> */
	
	tpem_subspace_reset (info);
	tpem_subspace_push (info, ket);

	/* calculate |j,1> = X|j,0> */
	tpem_sparse_product_tpem (matrix, info, jm1, jm0, eps);
	
	sum1 = tpem_conj_mult (jm0[ket], jm1[ket]);
	moment[1] += tpem_real (sum1);
	
	sum2 = tpem_cast (0.0, 0.0);
	for (p = info->stack; p < info->top; p++)
		sum2 = tpem_add (sum2, tpem_conj_mult (jm1[*p], jm1[*p]));
	moment[2] += tpem_real (sum2);

	/* calculate |j,m> = 2X|j,m-1> - |j,m-2>
	 *
	 * begin (m=2) pass	jm0 = |j,0>	jm1 = |j,1>
	 * end   (m=2) pass	jm0 = |j,1>	jm1 = |j,2>
	 * ...
	 * begin (m=k) pass	jm0 = |j,k-2>	jm1 = |j,k-1>
	 * end   (m=k) pass	jm0 = |j,k-1>	jm1 = |j,k>
	 */
	for (m = 2; m < n / 2; m++) {
		/* calculate |tmp> = X|jm1> */
		tpem_sparse_product_tpem (matrix, info, tmp, jm1, eps);
		sum1 = sum2 = tpem_cast (0.0, 0.0);
		for (p = info->stack; p < info->top; p++) {
			i = *p;
			keep = tpem_sub (tpem_add (tmp[i], tmp[i]), jm0[i]);
			/* for moment[2 * m    ] */
			sum1 = tpem_add (sum1, tpem_conj_mult (keep, keep));
			/* for moment[2 * m - 1] */
			sum2 = tpem_add (sum2, tpem_conj_mult (keep, jm1[i]));
			/* set |j,m-1> and |j,m-2> for next iteration */
			jm0[i] = jm1[i];
			jm1[i] = keep;
		}
		moment[m + m    ] += tpem_real (sum1);
		moment[m + m - 1] += tpem_real (sum2);
	}
	return info;
}

static void  tpem_diagonal_element_pem (tpem_sparse *matrix, size_t n,
		double *moment, size_t ket)
{	
	static	int	firstcall = 1;
	static	tpem_t	*tmp, *jm0, *jm1;
	size_t	i, j,m;
	tpem_t	sum1, sum2, keep;
	
	if (matrix == NULL || moment == NULL) {
		if (tmp) free (tmp), tmp = NULL;
		if (jm0) free (jm0), jm0 = NULL;
		if (jm1) free (jm1), jm1 = NULL;
		firstcall = 1;
		return;
	}
	if (firstcall) {
		firstcall = 0;
		tmp = tpem_malloc (matrix->rank * sizeof *tmp);
		jm0 = tpem_malloc (matrix->rank * sizeof *jm0);
		jm1 = tpem_malloc (matrix->rank * sizeof *jm1);
	} 
	
	for (i = 0; i < matrix->rank; i++)
		tmp[i] = jm0[i] = jm1[i] = tpem_cast (0.0, 0.0);
	jm0[ket] = tpem_cast (1.0, 0.0);	/* set |j,0> */
	
	
	/* calculate |j,1> = X|j,0> */
	tpem_sparse_product_pem (matrix, jm1, jm0);
	
	sum1 = tpem_conj_mult (jm0[ket], jm1[ket]);
	moment[1] += tpem_real (sum1);
	
	sum2 = tpem_cast (0.0, 0.0);
	for (j=0;j<matrix->rank;j++)
		sum2 = tpem_add (sum2, tpem_conj_mult (jm1[j], jm1[j]));
	moment[2] += tpem_real (sum2);

	/* calculate |j,m> = 2X|j,m-1> - |j,m-2>
	 *
	 * begin (m=2) pass	jm0 = |j,0>	jm1 = |j,1>
	 * end   (m=2) pass	jm0 = |j,1>	jm1 = |j,2>
	 * ...
	 * begin (m=k) pass	jm0 = |j,k-2>	jm1 = |j,k-1>
	 * end   (m=k) pass	jm0 = |j,k-1>	jm1 = |j,k>
	 */
	for (m = 2; m < n / 2; m++) {
		/* calculate |tmp> = X|jm1> */
		tpem_sparse_product_pem (matrix, tmp, jm1);
		sum1 = sum2 = tpem_cast (0.0, 0.0);
		for (i=0;i<matrix->rank;i++) {
			keep = tpem_sub (tpem_add (tmp[i], tmp[i]), jm0[i]);
			/* for moment[2 * m    ] */
			sum1 = tpem_add (sum1, tpem_conj_mult (keep, keep));
			/* for moment[2 * m - 1] */
			sum2 = tpem_add (sum2, tpem_conj_mult (keep, jm1[i]));
			/* set |j,m-1> and |j,m-2> for next iteration */
			jm0[i] = jm1[i];
			jm1[i] = keep;
		}
		moment[m + m    ] += tpem_real (sum1);
		moment[m + m - 1] += tpem_real (sum2);
	}
}


static tpem_subspace *tpem_subspace_for_trace (
		tpem_sparse *matrix0, tpem_sparse *matrix1,
		size_t n, double *moment0, double *moment1,
		size_t nsupp, size_t *support, double eps)
{	/* f77: mkTraceState */
	static	int	firstcall = 1;
	static	tpem_subspace	*info;
	tpem_subspace	*work;
	size_t	i, j, *p;

	if (matrix0 == NULL || matrix1 == NULL) {
		tpem_subspace_free (info);
		firstcall = 1;
		return NULL;
	}
	if (firstcall) {
		firstcall = 0;
		info = new_tpem_subspace (matrix0->rank);
	} else if (info->size != matrix0->rank) {
		tpem_finalize ();
		tpem_subspace_for_trace (matrix0, matrix1, n,
				moment0, moment1, nsupp, support, eps);
	}

	tpem_subspace_reset (info);
	for (i = 0; i < nsupp; i++) {
		j = support[i];
		work = tpem_diagonal_element_tpem (matrix0, n, moment0, j, eps);
		for (p = work->stack; p < work->top; p++)
			tpem_subspace_push (info, *p);
		work = tpem_diagonal_element_tpem (matrix1, n, moment1, j, eps);
		for (p = work->stack; p < work->top; p++)
			tpem_subspace_push (info, *p);
	}
	return info;
}


void tpem_calculate_moment_tpem (tpem_sparse *matrix, size_t n, double *moment,
		double eps)
{
	size_t	i;
	
	for (i = 0; i < n; i++)
		moment[i] = 0.0;
	moment[0] = (double) matrix->rank;
	
	for (i = 0; i < matrix->rank; i++)
		tpem_diagonal_element_tpem (matrix, n, moment, i, eps);

	for (i = 2; i < n; i += 2)
		moment[i] = 2.0 * moment[i] - moment[0];
	
	for (i = 3; i < n - 1; i += 2)
		moment[i] = 2.0 * moment[i] - moment[1];
}

void tpem_calculate_moment_pem (tpem_sparse *matrix, size_t n, double *moment)
{	
	size_t	i;
	
	for (i = 0; i < n; i++)
		moment[i] = 0.0;
	moment[0] = (double) matrix->rank;
	
	for (i = 0; i < matrix->rank; i++)
		tpem_diagonal_element_pem (matrix, n, moment, i);

	for (i = 2; i < n; i += 2)
		moment[i] = 2.0 * moment[i] - moment[0];
	
	for (i = 3; i < n - 1; i += 2)
		moment[i] = 2.0 * moment[i] - moment[1];
}

void tpem_calculate_moment_diff_tpem (tpem_sparse *matrix0, tpem_sparse *matrix1,
		size_t n, double *moment, size_t nsupp, size_t *support,
		double eps_trace, double eps_prod)
{	
	static	int	firstcall = 1;
	static	size_t	n_copy;
	static	double	*moment0, *moment1;
	tpem_subspace	*info;
	size_t	i, *p;

	if (matrix0 == NULL || matrix1 == NULL || moment == NULL) {
		if (moment0) free (moment0), moment0 = NULL;
		if (moment1) free (moment1), moment1 = NULL;
		n_copy = 0;
		firstcall = 1;
		return;
	}
	if (firstcall) {
		firstcall = 0;
		moment0 = tpem_malloc (n * sizeof *moment0);
		moment1 = tpem_malloc (n * sizeof *moment1);
		n_copy = n;
	} else if (n_copy != n) {
		tpem_finalize ();
		tpem_calculate_moment_diff_tpem (matrix0, matrix1, n, moment,
				nsupp, support, eps_trace, eps_prod);
	}

	info = tpem_subspace_for_trace (matrix0, matrix1, n, moment0, moment1,
			nsupp, support, eps_trace);

	for (i = 0; i < n; i++)
		moment0[i] = moment1[i] = 0.0;
	moment0[0] = moment1[0] = (double) matrix0->rank;

	for (p = info->stack; p < info->top; p++) {
		tpem_diagonal_element_tpem (matrix0, n, moment0, *p, eps_prod);
		tpem_diagonal_element_tpem (matrix1, n, moment1, *p, eps_prod);
	}

	for (i = 2; i < n; i += 2) {
		moment0[i] = 2.0 * moment0[i] - moment0[0];
		moment1[i] = 2.0 * moment1[i] - moment1[0];
	}

	for (i = 3; i < n - 1; i += 2) {
		moment0[i] = 2.0 * moment0[i] - moment0[1];
		moment1[i] = 2.0 * moment1[i] - moment1[1];
	}

	for (i = 0; i < n; i++)
		moment[i] = moment0[i] - moment1[i];

	
}

void tpem_calculate_moment_diff_pem (tpem_sparse *matrix0, tpem_sparse *matrix1,
		size_t n, double *moment)
{	
	static	int	firstcall = 1;
	static	size_t	n_copy;
	static	double	*moment0, *moment1;
	size_t	i, j;

	if (matrix0 == NULL || matrix1 == NULL || moment == NULL) {
		if (moment0) free (moment0), moment0 = NULL;
		if (moment1) free (moment1), moment1 = NULL;
		n_copy = 0;
		firstcall = 1;
		return;
	}
	if (firstcall) {
		firstcall = 0;
		moment0 = tpem_malloc (n * sizeof *moment0);
		moment1 = tpem_malloc (n * sizeof *moment1);
		n_copy = n;
	} else if (n_copy != n) {
		tpem_finalize ();
		tpem_calculate_moment_diff_pem (matrix0, matrix1, n, moment);
	}

	for (i = 0; i < n; i++)
		moment0[i] = moment1[i] = 0.0;
		moment0[0] = moment1[0] = (double) matrix0->rank;

	for (j=0;j<matrix0->rank;j++) {
		tpem_diagonal_element_pem (matrix0, n, moment0, j);
		tpem_diagonal_element_pem (matrix1, n, moment1, j);
	}

	for (i = 2; i < n; i += 2) {
		moment0[i] = 2.0 * moment0[i] - moment0[0];
		moment1[i] = 2.0 * moment1[i] - moment1[0];
	}

	for (i = 3; i < n - 1; i += 2) {
		moment0[i] = 2.0 * moment0[i] - moment0[1];
		moment1[i] = 2.0 * moment1[i] - moment1[1];
	}

	for (i = 0; i < n; i++)
		moment[i] = moment0[i] - moment1[i];

}



double tpem_expansion (size_t n, double *moments, double *coeffs) 
{
	size_t	i;
	double	ret = 0.0;
	
	for (i = 0; i < n; i++)
		ret += moments[i] * coeffs[i];
	return ret;
}


void tpem_finalize (void)
{
	tpem_diagonal_element_tpem (NULL, 0, NULL, 0, 0.0);
	tpem_subspace_for_trace (NULL, NULL, 0, NULL, NULL, 0, NULL, 0.0);
	tpem_calculate_moment_diff_tpem (NULL, NULL, 0, NULL, 0, NULL, 0.0, 0.0);
}

double factor_alpha(int m)
{
	if (m==0) return 1;
	return 2;
}

double my_f (double x, void * p) {
   struct my_f_params * params 
     = (struct my_f_params *)p;

   /* return  params->funk(params->m,x);*/
   double tmp;
   double tmp2= (double)1.0/M_PI;
   int m=params->m;
   tmp = params->funk(x) *  tmp2 * factor_alpha(m) * tpem_chebyshev(m,x)/sqrt(1.0-x*x);
   return tmp;
}

void tpem_calculate_coeffs (size_t FKWCutOff, double *vobs, tpem_func_ptr funk)
{
	int m;
	struct my_f_params params;
	double *pts,result,result1,result2,abserr;
	int npts,limit;
	double epsabs,epsrel;

#ifndef NO_GSL
	gsl_function f;
	gsl_integration_workspace * workspace;
	
	npts = 2;
	pts = malloc(sizeof(double)*2);
	pts[0]= -1.0;
	pts[1] = 1.0;
	epsabs=1e-9;
	epsrel=1e-9;
	limit = 1e6;
	workspace =  gsl_integration_workspace_alloc (limit+2);
	result=result1=result2=abserr=0;

	params.funk = funk;
	f.function= &my_f;
	f.params = &params;
	for (m=0;m<FKWCutOff;m++) {
		params.m = m;
		/* gsl_integration_qag (&f,pts[0],pts[1],epsabs,epsrel,limit,1,workspace,&result,&abserr);  */
		gsl_integration_qagp (&f,pts,npts,epsabs,epsrel,limit,workspace,&result,&abserr);
	 	/* gsl_integration_qags(&f,pts[0],pts[1],epsabs,epsrel,limit,workspace,&result,&abserr); */
		vobs[m] = result;
		
	}
	gsl_integration_workspace_free (workspace);s
#else
        fprintf(stderr,"GSL is required to compute tpem_calculate_coeffs\n");
        fprintf(stderr,"Use tpem_calculate_coeffs_alt instead\n");
        fprintf(stderr,"At this point %s %d\n",__FILE__,__LINE__);
        exit(1);
#endif
        free(pts);
}

void tpem_calculate_coeffs_alt (size_t CutOff, double *vobs, tpem_func_ptr funk)
{
  int m,k;
 
  double *f;
  double x, sum;
  size_t FKWCutOff = 2*CutOff;
 
  f = malloc(FKWCutOff * sizeof(*f));
 
  for(m = 0; m<FKWCutOff; m++ )
    {
      x = cos(M_PI*(m+0.5)/FKWCutOff);
      f[m] = (*funk)(x);
    }
  
  for(m = 0; m<CutOff; m++)
    {
      sum = 0.0;
      for(k = 0; k<FKWCutOff; k++)
        sum += f[k]*cos(M_PI * m*(k+0.5)/FKWCutOff);
      vobs[m] = 2.0 * sum/FKWCutOff;
    }
  vobs[0]/=2.0;
  free(f);
}

void tpem_calculate_coeffs_dos (size_t FKWCutOff, double *vobs, size_t num_x, size_t interval, size_t left, size_t right)
{
  int m,k,num_integral;
  double x, sum;
  if (2*FKWCutOff > num_x*interval)
	num_integral = 2*FKWCutOff;
  else
	num_integral = num_x*interval;

  for(m = 0; m<FKWCutOff; m++)
    {

      sum = 0.0;
      for(k = 0; k< num_integral; k++){
      	x = cos(M_PI*(k+0.5)/num_integral);
		if (x >= -1.0+2.0*left/num_x && x < -1.0+2.0*right/num_x){
        	sum += cos(M_PI * m*(k+0.5)/num_integral);
	  	}
	  }
      vobs[m] = 2.0 * sum/num_integral;
    }
  vobs[0]/=2.0;
}

double tpem_chebyshev(int m,double x) {
	double tmp;
	int p;
	if (m==0) return 1;
	if (m==1) return x;
	
	if ((m%2)==0) {
		p=m/2;
		tmp=tpem_chebyshev(p,x);
		return (2*tmp*tmp-1);
	}
	else {
		p=(m-1)/2;
		return (2*tpem_chebyshev(p,x)*tpem_chebyshev(p+1,x)-x);
	}
}	

void my_handler (const char * reason, const char * file, int line, int gsl_errno)
{
	fprintf(stderr,"WARNING (Ignore), %s, file=%s, line=%d, code=%d\n",
			reason,file,line,gsl_errno);
}

void tpem_init (void){
#ifndef NO_GSL
        gsl_set_error_handler (&my_handler);
#endif
}

