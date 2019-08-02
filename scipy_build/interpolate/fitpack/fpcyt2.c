/* fpcyt2.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int fpcyt2_(doublereal *a, integer *n, doublereal *b, 
	doublereal *c__, integer *nn)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    static integer i__, j, j1, n1;
    static doublereal cc, sum;

/* subroutine fpcyt2 solves a linear n x n system */
/*         a * c = b */
/* where matrix a is a cyclic tridiagonal matrix, decomposed */
/* using subroutine fpsyt1. */
/*  .. */
/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  .. */
    /* Parameter adjustments */
    --c__;
    --b;
    a_dim1 = *nn;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    c__[1] = b[1] * a[(a_dim1 << 2) + 1];
    sum = c__[1] * a[a_dim1 * 5 + 1];
    n1 = *n - 1;
    i__1 = n1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	c__[i__] = (b[i__] - a[i__ + a_dim1] * c__[i__ - 1]) * a[i__ + (
		a_dim1 << 2)];
	sum += c__[i__] * a[i__ + a_dim1 * 5];
/* L10: */
    }
    cc = (b[*n] - sum) * a[*n + (a_dim1 << 2)];
    c__[*n] = cc;
    c__[n1] -= cc * a[n1 + a_dim1 * 6];
    j = n1;
    i__1 = *n;
    for (i__ = 3; i__ <= i__1; ++i__) {
	j1 = j - 1;
	c__[j1] = c__[j1] - c__[j] * a[j1 + a_dim1 * 3] * a[j1 + (a_dim1 << 2)
		] - cc * a[j1 + a_dim1 * 6];
	j = j1;
/* L20: */
    }
    return 0;
} /* fpcyt2_ */

