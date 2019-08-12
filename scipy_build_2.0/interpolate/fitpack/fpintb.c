/* fpintb.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int fpintb_(doublereal *t, integer *n, doublereal *bint, 
	integer *nk1, doublereal *x, doublereal *y)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal a, b, f, h__[6];
    static integer i__, j, k, l;
    static doublereal h1[6];
    static integer j1, k1, l0, ia, ib;
    static doublereal ak;
    static integer li, lj, lk, it;
    static doublereal arg, one;
    static integer min__;
    static doublereal aint[6];

/*  subroutine fpintb calculates integrals of the normalized b-splines */
/*  nj,k+1(x) of degree k, defined on the set of knots t(j),j=1,2,...n. */
/*  it makes use of the formulae of gaffney for the calculation of */
/*  indefinite integrals of b-splines. */

/*  calling sequence: */
/*     call fpintb(t,n,bint,nk1,x,y) */

/*  input parameters: */
/*    t    : real array,length n, containing the position of the knots. */
/*    n    : integer value, giving the number of knots. */
/*    nk1  : integer value, giving the number of b-splines of degree k, */
/*           defined on the set of knots ,i.e. nk1 = n-k-1. */
/*    x,y  : real values, containing the end points of the integration */
/*           interval. */
/*  output parameter: */
/*    bint : array,length nk1, containing the integrals of the b-splines. */
/*  .. */
/*  ..scalars arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  ..local arrays.. */
/*  initialization. */
    /* Parameter adjustments */
    --t;
    --bint;

    /* Function Body */
    one = 1.;
    k1 = *n - *nk1;
    ak = (doublereal) k1;
    k = k1 - 1;
    i__1 = *nk1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	bint[i__] = 0.;
/* L10: */
    }
/*  the integration limits are arranged in increasing order. */
    a = *x;
    b = *y;
    min__ = 0;
    if (a < b) {
	goto L30;
    }
    if (a == b) {
	goto L160;
    }
    goto L20;
L20:
    a = *y;
    b = *x;
    min__ = 1;
L30:
    if (a < t[k1]) {
	a = t[k1];
    }
    if (b > t[*nk1 + 1]) {
	b = t[*nk1 + 1];
    }
/*  using the expression of gaffney for the indefinite integral of a */
/*  b-spline we find that */
/*  bint(j) = (t(j+k+1)-t(j))*(res(j,b)-res(j,a))/(k+1) */
/*    where for t(l) <= x < t(l+1) */
/*    res(j,x) = 0, j=1,2,...,l-k-1 */
/*             = 1, j=l+1,l+2,...,nk1 */
/*             = aint(j+k-l+1), j=l-k,l-k+1,...,l */
/*               = sumi((x-t(j+i))*nj+i,k+1-i(x)/(t(j+k+1)-t(j+i))) */
/*                 i=0,1,...,k */
    l = k1;
    l0 = l + 1;
/*  set arg = a. */
    arg = a;
    for (it = 1; it <= 2; ++it) {
/*  search for the knot interval t(l) <= arg < t(l+1). */
L40:
	if (arg < t[l0] || l == *nk1) {
	    goto L50;
	}
	l = l0;
	l0 = l + 1;
	goto L40;
/*  calculation of aint(j), j=1,2,...,k+1. */
/*  initialization. */
L50:
	i__1 = k1;
	for (j = 1; j <= i__1; ++j) {
	    aint[j - 1] = 0.;
/* L55: */
	}
	aint[0] = (arg - t[l]) / (t[l + 1] - t[l]);
	h1[0] = one;
	i__1 = k;
	for (j = 1; j <= i__1; ++j) {
/*  evaluation of the non-zero b-splines of degree j at arg,i.e. */
/*    h(i+1) = nl-j+i,j(arg), i=0,1,...,j. */
	    h__[0] = 0.;
	    i__2 = j;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		li = l + i__;
		lj = li - j;
		f = h1[i__ - 1] / (t[li] - t[lj]);
		h__[i__ - 1] += f * (t[li] - arg);
		h__[i__] = f * (arg - t[lj]);
/* L60: */
	    }
/*  updating of the integrals aint. */
	    j1 = j + 1;
	    i__2 = j1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		li = l + i__;
		lj = li - j1;
		aint[i__ - 1] += h__[i__ - 1] * (arg - t[lj]) / (t[li] - t[lj]
			);
		h1[i__ - 1] = h__[i__ - 1];
/* L70: */
	    }
	}
	if (it == 2) {
	    goto L100;
	}
/*  updating of the integrals bint */
	lk = l - k;
	ia = lk;
	i__2 = k1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    bint[lk] = -aint[i__ - 1];
	    ++lk;
/* L80: */
	}
/*  set arg = b. */
	arg = b;
/* L90: */
    }
/*  updating of the integrals bint. */
L100:
    lk = l - k;
    ib = lk - 1;
    i__2 = k1;
    for (i__ = 1; i__ <= i__2; ++i__) {
	bint[lk] += aint[i__ - 1];
	++lk;
/* L110: */
    }
    if (ib < ia) {
	goto L130;
    }
    i__2 = ib;
    for (i__ = ia; i__ <= i__2; ++i__) {
	bint[i__] += one;
/* L120: */
    }
/*  the scaling factors are taken into account. */
L130:
    f = one / ak;
    i__2 = *nk1;
    for (i__ = 1; i__ <= i__2; ++i__) {
	j = i__ + k1;
	bint[i__] = bint[i__] * (t[j] - t[i__]) * f;
/* L140: */
    }
/*  the order of the integration limits is taken into account. */
    if (min__ == 0) {
	goto L160;
    }
    i__2 = *nk1;
    for (i__ = 1; i__ <= i__2; ++i__) {
	bint[i__] = -bint[i__];
/* L150: */
    }
L160:
    return 0;
} /* fpintb_ */

