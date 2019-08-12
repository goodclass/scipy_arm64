/* fpadpo.f -- translated by f2c (version 20190311).
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

/* Table of constant values */

static integer c__0 = 0;

/* Subroutine */ int fpadpo_(integer *idim, doublereal *t, integer *n, 
	doublereal *c__, integer *nc, integer *k, doublereal *cp, integer *np,
	 doublereal *cc, doublereal *t1, doublereal *t2)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, l, k1, l1, n1, n2, ii, jj, nk1, nk2;
    extern /* Subroutine */ int fpinst_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *);

/*  given a idim-dimensional spline curve of degree k, in its b-spline */
/*  representation ( knots t(j),j=1,...,n , b-spline coefficients c(j), */
/*  j=1,...,nc) and given also a polynomial curve in its b-spline */
/*  representation ( coefficients cp(j), j=1,...,np), subroutine fpadpo */
/*  calculates the b-spline representation (coefficients c(j),j=1,...,nc) */
/*  of the sum of the two curves. */

/*  other subroutine required : fpinst */

/*  .. */
/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  .. */
    /* Parameter adjustments */
    --t2;
    --t1;
    --t;
    --cc;
    --c__;
    --cp;

    /* Function Body */
    k1 = *k + 1;
    nk1 = *n - k1;
/*  initialization */
    j = 1;
    l = 1;
    i__1 = *idim;
    for (jj = 1; jj <= i__1; ++jj) {
	l1 = j;
	i__2 = k1;
	for (ii = 1; ii <= i__2; ++ii) {
	    cc[l1] = cp[l];
	    ++l1;
	    ++l;
/* L10: */
	}
	j += *n;
	l += k1;
/* L20: */
    }
    if (nk1 == k1) {
	goto L70;
    }
    n1 = k1 << 1;
    j = *n;
    l = n1;
    i__1 = k1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t1[i__] = t[i__];
	t1[l] = t[j];
	--l;
	--j;
/* L30: */
    }
/*  find the b-spline representation of the given polynomial curve */
/*  according to the given set of knots. */
    nk2 = nk1 - 1;
    i__1 = nk2;
    for (l = k1; l <= i__1; ++l) {
	l1 = l + 1;
	j = 1;
	i__2 = *idim;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    fpinst_(&c__0, &t1[1], &n1, &cc[j], k, &t[l1], &l, &t2[1], &n2, &
		    cc[j], n);
	    j += *n;
/* L40: */
	}
	i__2 = n2;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    t1[i__] = t2[i__];
/* L50: */
	}
	n1 = n2;
/* L60: */
    }
/*  find the b-spline representation of the resulting curve. */
L70:
    j = 1;
    i__1 = *idim;
    for (jj = 1; jj <= i__1; ++jj) {
	l = j;
	i__2 = nk1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    c__[l] = cc[l] + c__[l];
	    ++l;
/* L80: */
	}
	j += *n;
/* L90: */
    }
    return 0;
} /* fpadpo_ */

