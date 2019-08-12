/* spalde.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int spalde_(doublereal *t, integer *n, doublereal *c__, 
	integer *k1, doublereal *x, doublereal *d__, integer *ier)
{
    static integer l, nk1;
    extern /* Subroutine */ int fpader_(doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *);

/*  subroutine spalde evaluates at a point x all the derivatives */
/*              (j-1) */
/*      d(j) = s     (x) , j=1,2,...,k1 */
/*  of a spline s(x) of order k1 (degree k=k1-1), given in its b-spline */
/*  representation. */

/*  calling sequence: */
/*     call spalde(t,n,c,k1,x,d,ier) */

/*  input parameters: */
/*    t    : array,length n, which contains the position of the knots. */
/*    n    : integer, giving the total number of knots of s(x). */
/*    c    : array,length n, which contains the b-spline coefficients. */
/*    k1   : integer, giving the order of s(x) (order=degree+1) */
/*    x    : real, which contains the point where the derivatives must */
/*           be evaluated. */

/*  output parameters: */
/*    d    : array,length k1, containing the derivative values of s(x). */
/*    ier  : error flag */
/*      ier = 0 : normal return */
/*      ier =10 : invalid input data (see restrictions) */

/*  restrictions: */
/*    t(k1) <= x <= t(n-k1+1) */

/*  further comments: */
/*    if x coincides with a knot, right derivatives are computed */
/*    ( left derivatives if x = t(n-k1+1) ). */

/*  other subroutines required: fpader. */

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

/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  .. */
/*  before starting computations a data check is made. if the input data */
/*  are invalid control is immediately repassed to the calling program. */
    /* Parameter adjustments */
    --c__;
    --t;
    --d__;

    /* Function Body */
    *ier = 10;
    nk1 = *n - *k1;
    if (*x < t[*k1] || *x > t[nk1 + 1]) {
	goto L300;
    }
/*  search for knot interval t(l) <= x < t(l+1) */
    l = *k1;
L100:
    if (*x < t[l + 1] || l == nk1) {
	goto L200;
    }
    ++l;
    goto L100;
L200:
    if (t[l] >= t[l + 1]) {
	goto L300;
    }
    *ier = 0;
/*  calculate the derivatives. */
    fpader_(&t[1], n, &c__[1], k1, x, &l, &d__[1]);
L300:
    return 0;
} /* spalde_ */

