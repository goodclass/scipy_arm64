/* fpchec.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int fpchec_(doublereal *x, integer *m, doublereal *t, 
	integer *n, integer *k, integer *ier)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, l, k1, k2;
    static doublereal tj, tl;
    static integer nk1, nk2, nk3;

/*  subroutine fpchec verifies the number and the position of the knots */
/*  t(j),j=1,2,...,n of a spline of degree k, in relation to the number */
/*  and the position of the data points x(i),i=1,2,...,m. if all of the */
/*  following conditions are fulfilled, the error parameter ier is set */
/*  to zero. if one of the conditions is violated ier is set to ten. */
/*      1) k+1 <= n-k-1 <= m */
/*      2) t(1) <= t(2) <= ... <= t(k+1) */
/*         t(n-k) <= t(n-k+1) <= ... <= t(n) */
/*      3) t(k+1) < t(k+2) < ... < t(n-k) */
/*      4) t(k+1) <= x(i) <= t(n-k) */
/*      5) the conditions specified by schoenberg and whitney must hold */
/*         for at least one subset of data points, i.e. there must be a */
/*         subset of data points y(j) such that */
/*             t(j) < y(j) < t(j+k+1), j=1,2,...,n-k-1 */
/*  .. */
/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  .. */
    /* Parameter adjustments */
    --x;
    --t;

    /* Function Body */
    k1 = *k + 1;
    k2 = k1 + 1;
    nk1 = *n - k1;
    nk2 = nk1 + 1;
    *ier = 10;
/*  check condition no 1 */
    if (nk1 < k1 || nk1 > *m) {
	goto L80;
    }
/*  check condition no 2 */
    j = *n;
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (t[i__] > t[i__ + 1]) {
	    goto L80;
	}
	if (t[j] < t[j - 1]) {
	    goto L80;
	}
	--j;
/* L20: */
    }
/*  check condition no 3 */
    i__1 = nk2;
    for (i__ = k2; i__ <= i__1; ++i__) {
	if (t[i__] <= t[i__ - 1]) {
	    goto L80;
	}
/* L30: */
    }
/*  check condition no 4 */
    if (x[1] < t[k1] || x[*m] > t[nk2]) {
	goto L80;
    }
/*  check condition no 5 */
    if (x[1] >= t[k2] || x[*m] <= t[nk1]) {
	goto L80;
    }
    i__ = 1;
    l = k2;
    nk3 = nk1 - 1;
    if (nk3 < 2) {
	goto L70;
    }
    i__1 = nk3;
    for (j = 2; j <= i__1; ++j) {
	tj = t[j];
	++l;
	tl = t[l];
L40:
	++i__;
	if (i__ >= *m) {
	    goto L80;
	}
	if (x[i__] <= tj) {
	    goto L40;
	}
	if (x[i__] >= tl) {
	    goto L80;
	}
/* L60: */
    }
L70:
    *ier = 0;
L80:
    return 0;
} /* fpchec_ */

