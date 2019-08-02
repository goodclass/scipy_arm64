/* splev.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int splev_(doublereal *t, integer *n, doublereal *c__, 
	integer *k, doublereal *x, doublereal *y, integer *m, integer *e, 
	integer *ier)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal h__[20];
    static integer i__, j, l, k1, l1, k2;
    static doublereal tb;
    static integer ll;
    static doublereal te, sp;
    static integer nk1;
    static doublereal arg;
    extern /* Subroutine */ int fpbspl_(doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *);

/*  subroutine splev evaluates in a number of points x(i),i=1,2,...,m */
/*  a spline s(x) of degree k, given in its b-spline representation. */

/*  calling sequence: */
/*     call splev(t,n,c,k,x,y,m,e,ier) */

/*  input parameters: */
/*    t    : array,length n, which contains the position of the knots. */
/*    n    : integer, giving the total number of knots of s(x). */
/*    c    : array,length n, which contains the b-spline coefficients. */
/*    k    : integer, giving the degree of s(x). */
/*    x    : array,length m, which contains the points where s(x) must */
/*           be evaluated. */
/*    m    : integer, giving the number of points where s(x) must be */
/*           evaluated. */
/*    e    : integer, if 0 the spline is extrapolated from the end */
/*           spans for points not in the support, if 1 the spline */
/*           evaluates to zero for those points, if 2 ier is set to */
/*           1 and the subroutine returns, and if 3 the spline evaluates */
/*           to the value of the nearest boundary point. */

/*  output parameter: */
/*    y    : array,length m, giving the value of s(x) at the different */
/*           points. */
/*    ier  : error flag */
/*      ier = 0 : normal return */
/*      ier = 1 : argument out of bounds and e == 2 */
/*      ier =10 : invalid input data (see restrictions) */

/*  restrictions: */
/*    m >= 1 */
/* --    t(k+1) <= x(i) <= x(i+1) <= t(n-k) , i=1,2,...,m-1. */

/*  other subroutines required: fpbspl. */

/*  references : */
/*    de boor c  : on calculating with b-splines, j. approximation theory */
/*                 6 (1972) 50-62. */
/*    cox m.g.   : the numerical evaluation of b-splines, j. inst. maths */
/*                 applics 10 (1972) 134-149. */
/*    dierckx p. : curve and surface fitting with splines, monographs on */
/*                 numerical analysis, oxford university press, 1993. */

/*  author : */
/*    p.dierckx */
/*    dept. computer science, k.u.leuven */
/*    celestijnenlaan 200a, b-3001 heverlee, belgium. */
/*    e-mail : Paul.Dierckx@cs.kuleuven.ac.be */

/*  latest update : march 1987 */

/* ++ pearu: 11 aug 2003 */
/* ++   - disabled cliping x values to interval [min(t),max(t)] */
/* ++   - removed the restriction of the orderness of x values */
/* ++   - fixed initialization of sp to double precision value */

/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/* ++.. */
/* ..++ */
/*  ..local array.. */
/*  .. */
/*  before starting computations a data check is made. if the input data */
/*  are invalid control is immediately repassed to the calling program. */
    /* Parameter adjustments */
    --c__;
    --t;
    --y;
    --x;

    /* Function Body */
    *ier = 10;
/* --      if(m-1) 100,30,10 */
/* ++.. */
    if (*m < 1) {
	goto L100;
    }
/* ..++ */
/* --  10  do 20 i=2,m */
/* --        if(x(i).lt.x(i-1)) go to 100 */
/* --  20  continue */
    *ier = 0;
/*  fetch tb and te, the boundaries of the approximation interval. */
    k1 = *k + 1;
/* ++.. */
    k2 = k1 + 1;
/* ..++ */
    nk1 = *n - k1;
    tb = t[k1];
    te = t[nk1 + 1];
    l = k1;
    l1 = l + 1;
/*  main loop for the different points. */
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*  fetch a new x-value arg. */
	arg = x[i__];
/*  check if arg is in the support */
	if (arg < tb || arg > te) {
	    if (*e == 0) {
		goto L35;
	    } else if (*e == 1) {
		y[i__] = 0.;
		goto L80;
	    } else if (*e == 2) {
		*ier = 1;
		goto L100;
	    } else if (*e == 3) {
		if (arg < tb) {
		    arg = tb;
		} else {
		    arg = te;
		}
	    }
	}
/*  search for knot interval t(l) <= arg < t(l+1) */
/* ++.. */
L35:
	if (arg >= t[l] || l1 == k2) {
	    goto L40;
	}
	l1 = l;
	--l;
	goto L35;
/* ..++ */
L40:
	if (arg < t[l1] || l == nk1) {
	    goto L50;
	}
	l = l1;
	l1 = l + 1;
	goto L40;
/*  evaluate the non-zero b-splines at arg. */
L50:
	fpbspl_(&t[1], n, k, &arg, &l, h__);
/*  find the value of s(x) at x=arg. */
	sp = 0.;
	ll = l - k1;
	i__2 = k1;
	for (j = 1; j <= i__2; ++j) {
	    ++ll;
	    sp += c__[ll] * h__[j - 1];
/* L60: */
	}
	y[i__] = sp;
L80:
	;
    }
L100:
    return 0;
} /* splev_ */

