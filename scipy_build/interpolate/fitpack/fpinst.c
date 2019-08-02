/* fpinst.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int fpinst_(integer *iopt, doublereal *t, integer *n, 
	doublereal *c__, integer *k, doublereal *x, integer *l, doublereal *
	tt, integer *nn, doublereal *cc, integer *nest)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, m, i1, k1, mk, nk, nl, ll, nk1;
    static doublereal fac, one, per;

/*  given the b-spline representation (knots t(j),j=1,2,...,n, b-spline */
/*  coefficients c(j),j=1,2,...,n-k-1) of a spline of degree k, fpinst */
/*  calculates the b-spline representation (knots tt(j),j=1,2,...,nn, */
/*  b-spline coefficients cc(j),j=1,2,...,nn-k-1) of the same spline if */
/*  an additional knot is inserted at the point x situated in the inter- */
/*  val t(l)<=x<t(l+1). iopt denotes whether (iopt.ne.0) or not (iopt=0) */
/*  the given spline is periodic. in case of a periodic spline at least */
/*  one of the following conditions must be fulfilled: l>2*k or l<n-2*k. */

/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  .. */
    /* Parameter adjustments */
    --cc;
    --tt;
    --c__;
    --t;

    /* Function Body */
    one = 1.f;
    k1 = *k + 1;
    nk1 = *n - k1;
/*  the new knots */
    ll = *l + 1;
    i__ = *n;
    i__1 = *n;
    for (j = ll; j <= i__1; ++j) {
	tt[i__ + 1] = t[i__];
	--i__;
/* L10: */
    }
    tt[ll] = *x;
    i__1 = *l;
    for (j = 1; j <= i__1; ++j) {
	tt[j] = t[j];
/* L20: */
    }
/*  the new b-spline coefficients */
    i__ = nk1;
    i__1 = nk1;
    for (j = *l; j <= i__1; ++j) {
	cc[i__ + 1] = c__[i__];
	--i__;
/* L30: */
    }
    i__ = *l;
    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
	m = i__ + k1;
	fac = (*x - tt[i__]) / (tt[m] - tt[i__]);
	i1 = i__ - 1;
	cc[i__] = fac * c__[i__] + (one - fac) * c__[i1];
	i__ = i1;
/* L40: */
    }
    i__1 = i__;
    for (j = 1; j <= i__1; ++j) {
	cc[j] = c__[j];
/* L50: */
    }
    *nn = *n + 1;
    if (*iopt == 0) {
	return 0;
    }
/*   incorporate the boundary conditions for a periodic spline. */
    nk = *nn - *k;
    nl = nk - k1;
    per = tt[nk] - tt[k1];
    i__ = k1;
    j = nk;
    if (ll <= nl) {
	goto L70;
    }
    i__1 = *k;
    for (m = 1; m <= i__1; ++m) {
	mk = m + nl;
	cc[m] = cc[mk];
	--i__;
	--j;
	tt[i__] = tt[j] - per;
/* L60: */
    }
    return 0;
L70:
    if (ll > k1 + *k) {
	return 0;
    }
    i__1 = *k;
    for (m = 1; m <= i__1; ++m) {
	mk = m + nl;
	cc[mk] = cc[m];
	++i__;
	++j;
	tt[j] = tt[i__] + per;
/* L80: */
    }
    return 0;
} /* fpinst_ */

