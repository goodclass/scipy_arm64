/* fpader.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int fpader_(doublereal *t, integer *n, doublereal *c__, 
	integer *k1, doublereal *x, integer *l, doublereal *d__)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static doublereal h__[20];
    static integer i__, j, j1, j2;
    static doublereal ak;
    static integer ik, jj, ki, kj, li, lj, lk;
    static doublereal fac, one;

/*  subroutine fpader calculates the derivatives */
/*             (j-1) */
/*     d(j) = s     (x) , j=1,2,...,k1 */
/*  of a spline of order k1 at the point t(l)<=x<t(l+1), using the */
/*  stable recurrence scheme of de boor */
/*  .. */
/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  ..local array.. */
/*  .. */
    /* Parameter adjustments */
    --c__;
    --t;
    --d__;

    /* Function Body */
    one = 1.;
    lk = *l - *k1;
    i__1 = *k1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ik = i__ + lk;
	h__[i__ - 1] = c__[ik];
/* L100: */
    }
    kj = *k1;
    fac = one;
    i__1 = *k1;
    for (j = 1; j <= i__1; ++j) {
	ki = kj;
	j1 = j + 1;
	if (j == 1) {
	    goto L300;
	}
	i__ = *k1;
	i__2 = *k1;
	for (jj = j; jj <= i__2; ++jj) {
	    li = i__ + lk;
	    lj = li + kj;
	    h__[i__ - 1] = (h__[i__ - 1] - h__[i__ - 2]) / (t[lj] - t[li]);
	    --i__;
/* L200: */
	}
L300:
	i__2 = *k1;
	for (i__ = j; i__ <= i__2; ++i__) {
	    d__[i__] = h__[i__ - 1];
/* L400: */
	}
	if (j == *k1) {
	    goto L600;
	}
	i__2 = *k1;
	for (jj = j1; jj <= i__2; ++jj) {
	    --ki;
	    i__ = *k1;
	    i__3 = *k1;
	    for (j2 = jj; j2 <= i__3; ++j2) {
		li = i__ + lk;
		lj = li + ki;
		d__[i__] = ((*x - t[li]) * d__[i__] + (t[lj] - *x) * d__[i__ 
			- 1]) / (t[lj] - t[li]);
		--i__;
/* L500: */
	    }
	}
L600:
	d__[j] = d__[*k1] * fac;
	ak = (doublereal) (*k1 - j);
	fac *= ak;
	--kj;
/* L700: */
    }
    return 0;
} /* fpader_ */

