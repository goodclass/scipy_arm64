/* splder.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int splder_(doublereal *t, integer *n, doublereal *c__, 
	integer *k, integer *nu, doublereal *x, doublereal *y, integer *m, 
	integer *e, doublereal *wrk, integer *ier)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal h__[6];
    static integer i__, j, l, k1, k2, l1, l2, k3;
    static doublereal ak;
    static integer kk;
    static doublereal tb;
    static integer ll;
    static doublereal te;
    static integer nn;
    static doublereal sp;
    static integer nk1, nk2;
    static doublereal fac, arg;
    extern /* Subroutine */ int fpbspl_(doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *);

/*  subroutine splder evaluates in a number of points x(i),i=1,2,...,m */
/*  the derivative of order nu of a spline s(x) of degree k,given in */
/*  its b-spline representation. */

/*  calling sequence: */
/*     call splder(t,n,c,k,nu,x,y,m,e,wrk,ier) */

/*  input parameters: */
/*    t    : array,length n, which contains the position of the knots. */
/*    n    : integer, giving the total number of knots of s(x). */
/*    c    : array,length n, which contains the b-spline coefficients. */
/*    k    : integer, giving the degree of s(x). */
/*    nu   : integer, specifying the order of the derivative. 0<=nu<=k */
/*    x    : array,length m, which contains the points where the deriv- */
/*           ative of s(x) must be evaluated. */
/*    m    : integer, giving the number of points where the derivative */
/*           of s(x) must be evaluated */
/*    e    : integer, if 0 the spline is extrapolated from the end */
/*           spans for points not in the support, if 1 the spline */
/*           evaluates to zero for those points, and if 2 ier is set to */
/*           1 and the subroutine returns. */
/*    wrk  : real array of dimension n. used as working space. */

/*  output parameters: */
/*    y    : array,length m, giving the value of the derivative of s(x) */
/*           at the different points. */
/*    ier  : error flag */
/*      ier = 0 : normal return */
/*      ier = 1 : argument out of bounds and e == 2 */
/*      ier =10 : invalid input data (see restrictions) */

/*  restrictions: */
/*    0 <= nu <= k */
/*    m >= 1 */
/*    t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1. */

/*  other subroutines required: fpbspl */

/*  references : */
/*    de boor c : on calculating with b-splines, j. approximation theory */
/*                6 (1972) 50-62. */
/*    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths */
/*                applics 10 (1972) 134-149. */
/*   dierckx p. : curve and surface fitting with splines, monographs on */
/*                numerical analysis, oxford university press, 1993. */

/*  author : */
/*    p.dierckx */
/*    dept. computer science, k.u.leuven */
/*    celestijnenlaan 200a, b-3001 heverlee, belgium. */
/*    e-mail : Paul.Dierckx@cs.kuleuven.ac.be */

/*  latest update : march 1987 */

/* ++ pearu: 13 aug 20003 */
/* ++   - disabled cliping x values to interval [min(t),max(t)] */
/* ++   - removed the restriction of the orderness of x values */
/* ++   - fixed initialization of sp to double precision value */

/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/* ++.. */
/* ..++ */
/*  ..local arrays .. */
/*  before starting computations a data check is made. if the input data */
/*  are invalid control is immediately repassed to the calling program. */
    /* Parameter adjustments */
    --wrk;
    --c__;
    --t;
    --y;
    --x;

    /* Function Body */
    *ier = 10;
    if (*nu < 0 || *nu > *k) {
	goto L200;
    }
/* --      if(m-1) 200,30,10 */
/* ++.. */
    if (*m < 1) {
	goto L200;
    }
/* ..++ */
/* --  10  do 20 i=2,m */
/* --        if(x(i).lt.x(i-1)) go to 200 */
/* --  20  continue */
    *ier = 0;
/*  fetch tb and te, the boundaries of the approximation interval. */
    k1 = *k + 1;
    k3 = k1 + 1;
    nk1 = *n - k1;
    tb = t[k1];
    te = t[nk1 + 1];
/*  the derivative of order nu of a spline of degree k is a spline of */
/*  degree k-nu,the b-spline coefficients wrk(i) of which can be found */
/*  using the recurrence scheme of de boor. */
    l = 1;
    kk = *k;
    nn = *n;
    i__1 = nk1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wrk[i__] = c__[i__];
/* L40: */
    }
    if (*nu == 0) {
	goto L100;
    }
    nk2 = nk1;
    i__1 = *nu;
    for (j = 1; j <= i__1; ++j) {
	ak = (doublereal) kk;
	--nk2;
	l1 = l;
	i__2 = nk2;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ++l1;
	    l2 = l1 + kk;
	    fac = t[l2] - t[l1];
	    if (fac <= 0.f) {
		goto L50;
	    }
	    wrk[i__] = ak * (wrk[i__ + 1] - wrk[i__]) / fac;
L50:
	    ;
	}
	++l;
	--kk;
/* L60: */
    }
    if (kk != 0) {
	goto L100;
    }
/*  if nu=k the derivative is a piecewise constant function */
    j = 1;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	arg = x[i__];
/* ++.. */
/*  check if arg is in the support */
	if (arg < tb || arg > te) {
	    if (*e == 0) {
		goto L65;
	    } else if (*e == 1) {
		y[i__] = 0.;
		goto L90;
	    } else if (*e == 2) {
		*ier = 1;
		goto L200;
	    }
	}
/*  search for knot interval t(l) <= arg < t(l+1) */
L65:
	if (arg >= t[l] || l + 1 == k3) {
	    goto L70;
	}
	l1 = l;
	--l;
	--j;
	goto L65;
/* ..++ */
L70:
	if (arg < t[l + 1] || l == nk1) {
	    goto L80;
	}
	++l;
	++j;
	goto L70;
L80:
	y[i__] = wrk[j];
L90:
	;
    }
    goto L200;
L100:
    l = k1;
    l1 = l + 1;
    k2 = k1 - *nu;
/*  main loop for the different points. */
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*  fetch a new x-value arg. */
	arg = x[i__];
/*  check if arg is in the support */
	if (arg < tb || arg > te) {
	    if (*e == 0) {
		goto L135;
	    } else if (*e == 1) {
		y[i__] = 0.;
		goto L180;
	    } else if (*e == 2) {
		*ier = 1;
		goto L200;
	    }
	}
/*  search for knot interval t(l) <= arg < t(l+1) */
L135:
	if (arg >= t[l] || l1 == k3) {
	    goto L140;
	}
	l1 = l;
	--l;
	goto L135;
/* ..++ */
L140:
	if (arg < t[l1] || l == nk1) {
	    goto L150;
	}
	l = l1;
	l1 = l + 1;
	goto L140;
/*  evaluate the non-zero b-splines of degree k-nu at arg. */
L150:
	fpbspl_(&t[1], n, &kk, &arg, &l, h__);
/*  find the value of the derivative at x=arg. */
	sp = 0.;
	ll = l - k1;
	i__2 = k2;
	for (j = 1; j <= i__2; ++j) {
	    ++ll;
	    sp += wrk[ll] * h__[j - 1];
/* L160: */
	}
	y[i__] = sp;
L180:
	;
    }
L200:
    return 0;
} /* splder_ */

