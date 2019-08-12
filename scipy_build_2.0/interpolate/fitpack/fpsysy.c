/* fpsysy.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int fpsysy_(doublereal *a, integer *n, doublereal *g)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, i1;
    static doublereal fac;

/* subroutine fpsysy solves a linear n x n symmetric system */
/*    (a) * (b) = (g) */
/* on input, vector g contains the right hand side ; on output it will */
/* contain the solution (b). */
/*  .. */
/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  .. */
    /* Parameter adjustments */
    --g;
    a -= 7;

    /* Function Body */
    g[1] /= a[7];
    if (*n == 1) {
	return 0;
    }
/*  decomposition of the symmetric matrix (a) = (l) * (d) *(l)' */
/*  with (l) a unit lower triangular matrix and (d) a diagonal */
/*  matrix */
    i__1 = *n;
    for (k = 2; k <= i__1; ++k) {
	a[k + 6] /= a[7];
/* L10: */
    }
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	i1 = i__ - 1;
	i__2 = *n;
	for (k = i__; k <= i__2; ++k) {
	    fac = a[k + i__ * 6];
	    i__3 = i1;
	    for (j = 1; j <= i__3; ++j) {
		fac -= a[j + j * 6] * a[k + j * 6] * a[i__ + j * 6];
/* L20: */
	    }
	    a[k + i__ * 6] = fac;
	    if (k > i__) {
		a[k + i__ * 6] = fac / a[i__ + i__ * 6];
	    }
/* L30: */
	}
/* L40: */
    }
/*  solve the system (l)*(d)*(l)'*(b) = (g). */
/*  first step : solve (l)*(d)*(c) = (g). */
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	i1 = i__ - 1;
	fac = g[i__];
	i__2 = i1;
	for (j = 1; j <= i__2; ++j) {
	    fac -= g[j] * a[j + j * 6] * a[i__ + j * 6];
/* L50: */
	}
	g[i__] = fac / a[i__ + i__ * 6];
/* L60: */
    }
/*  second step : solve (l)'*(b) = (c) */
    i__ = *n;
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
	i1 = i__;
	--i__;
	fac = g[i__];
	i__2 = *n;
	for (k = i1; k <= i__2; ++k) {
	    fac -= g[k] * a[k + i__ * 6];
/* L70: */
	}
	g[i__] = fac;
/* L80: */
    }
    return 0;
} /* fpsysy_ */

