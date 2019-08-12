/* bnorm.f -- translated by f2c (version 20190311).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

doublereal bnorm_(integer *n, doublereal *a, integer *nra, integer *ml, 
	integer *mu, doublereal *w)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal ret_val, d__1, d__2;

    /* Local variables */
    static integer i__, j, i1;
    static doublereal an;
    static integer jhi, jlo;
    static doublereal sum;

/* lll. optimize */
/* ----------------------------------------------------------------------- */
/* this function computes the norm of a banded n by n matrix, */
/* stored in the array a, that is consistent with the weighted max-norm */
/* on vectors, with weights stored in the array w. */
/* ml and mu are the lower and upper half-bandwidths of the matrix. */
/* nra is the first dimension of the a array, nra .ge. ml+mu+1. */
/* in terms of the matrix elements a(i,j), the norm is given by.. */
/*   bnorm = max(i=1,...,n) ( w(i) * sum(j=1,...,n) abs(a(i,j))/w(j) ) */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --w;
    a_dim1 = *nra;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    an = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = 0.;
	i1 = i__ + *mu + 1;
/* Computing MAX */
	i__2 = i__ - *ml;
	jlo = max(i__2,1);
/* Computing MIN */
	i__2 = i__ + *mu;
	jhi = min(i__2,*n);
	i__2 = jhi;
	for (j = jlo; j <= i__2; ++j) {
/* L10: */
	    sum += (d__1 = a[i1 - j + j * a_dim1], abs(d__1)) / w[j];
	}
/* Computing MAX */
	d__1 = an, d__2 = sum * w[i__];
	an = max(d__1,d__2);
/* L20: */
    }
    ret_val = an;
    return ret_val;
/* ----------------------- end of function bnorm ------------------------- */
} /* bnorm_ */

