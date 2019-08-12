/* splint.f -- translated by f2c (version 20190311).
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

doublereal splint_(doublereal *t, integer *n, doublereal *c__, integer *k, 
	doublereal *a, doublereal *b, doublereal *wrk)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static integer i__, nk1;
    extern /* Subroutine */ int fpintb_(doublereal *, integer *, doublereal *,
	     integer *, doublereal *, doublereal *);

/*  function splint calculates the integral of a spline function s(x) */
/*  of degree k, which is given in its normalized b-spline representation */

/*  calling sequence: */
/*     aint = splint(t,n,c,k,a,b,wrk) */

/*  input parameters: */
/*    t    : array,length n,which contains the position of the knots */
/*           of s(x). */
/*    n    : integer, giving the total number of knots of s(x). */
/*    c    : array,length n, containing the b-spline coefficients. */
/*    k    : integer, giving the degree of s(x). */
/*    a,b  : real values, containing the end points of the integration */
/*           interval. s(x) is considered to be identically zero outside */
/*           the interval (t(k+1),t(n-k)). */

/*  output parameter: */
/*    aint : real, containing the integral of s(x) between a and b. */
/*    wrk  : real array, length n.  used as working space */
/*           on output, wrk will contain the integrals of the normalized */
/*           b-splines defined on the set of knots. */

/*  other subroutines required: fpintb. */

/*  references : */
/*    gaffney p.w. : the calculation of indefinite integrals of b-splines */
/*                   j. inst. maths applics 17 (1976) 37-41. */
/*    dierckx p. : curve and surface fitting with splines, monographs on */
/*                 numerical analysis, oxford university press, 1993. */

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
    /* Parameter adjustments */
    --wrk;
    --c__;
    --t;

    /* Function Body */
    nk1 = *n - *k - 1;
/*  calculate the integrals wrk(i) of the normalized b-splines */
/*  ni,k+1(x), i=1,2,...nk1. */
    fpintb_(&t[1], n, &wrk[1], &nk1, a, b);
/*  calculate the integral of s(x). */
    ret_val = 0.;
    i__1 = nk1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ret_val += c__[i__] * wrk[i__];
/* L10: */
    }
    return ret_val;
} /* splint_ */

