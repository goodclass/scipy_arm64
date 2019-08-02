/* cualde.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int cualde_(integer *idim, doublereal *t, integer *n, 
	doublereal *c__, integer *nc, integer *k1, doublereal *u, doublereal *
	d__, integer *nd, integer *ier)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal h__[6];
    static integer i__, j, l, m, kk, nk1;
    extern /* Subroutine */ int fpader_(doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *);

/*  subroutine cualde evaluates at the point u all the derivatives */
/*                     (l) */
/*     d(idim*l+j) = sj   (u) ,l=0,1,...,k, j=1,2,...,idim */
/*  of a spline curve s(u) of order k1 (degree k=k1-1) and dimension idim */
/*  given in its b-spline representation. */

/*  calling sequence: */
/*     call cualde(idim,t,n,c,nc,k1,u,d,nd,ier) */

/*  input parameters: */
/*    idim : integer, giving the dimension of the spline curve. */
/*    t    : array,length n, which contains the position of the knots. */
/*    n    : integer, giving the total number of knots of s(u). */
/*    c    : array,length nc, which contains the b-spline coefficients. */
/*    nc   : integer, giving the total number of coefficients of s(u). */
/*    k1   : integer, giving the order of s(u) (order=degree+1). */
/*    u    : real, which contains the point where the derivatives must */
/*           be evaluated. */
/*    nd   : integer, giving the dimension of the array d. nd >= k1*idim */

/*  output parameters: */
/*    d    : array,length nd,giving the different curve derivatives. */
/*           d(idim*l+j) will contain the j-th coordinate of the l-th */
/*           derivative of the curve at the point u. */
/*    ier  : error flag */
/*      ier = 0 : normal return */
/*      ier =10 : invalid input data (see restrictions) */

/*  restrictions: */
/*    nd >= k1*idim */
/*    t(k1) <= u <= t(n-k1+1) */

/*  further comments: */
/*    if u coincides with a knot, right derivatives are computed */
/*    ( left derivatives if u = t(n-k1+1) ). */

/*  other subroutines required: fpader. */

/*  references : */
/*    de boor c : on calculating with b-splines, j. approximation theory */
/*                6 (1972) 50-62. */
/*    cox m.g.  : the numerical evaluation of b-splines, j. inst. maths */
/*                applics 10 (1972) 134-149. */
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
/*  ..local array.. */
/*  .. */
/*  before starting computations a data check is made. if the input data */
/*  are invalid control is immediately repassed to the calling program. */
    /* Parameter adjustments */
    --t;
    --c__;
    --d__;

    /* Function Body */
    *ier = 10;
    if (*nd < *k1 * *idim) {
	goto L500;
    }
    nk1 = *n - *k1;
    if (*u < t[*k1] || *u > t[nk1 + 1]) {
	goto L500;
    }
/*  search for knot interval t(l) <= u < t(l+1) */
    l = *k1;
L100:
    if (*u < t[l + 1] || l == nk1) {
	goto L200;
    }
    ++l;
    goto L100;
L200:
    if (t[l] >= t[l + 1]) {
	goto L500;
    }
    *ier = 0;
/*  calculate the derivatives. */
    j = 1;
    i__1 = *idim;
    for (i__ = 1; i__ <= i__1; ++i__) {
	fpader_(&t[1], n, &c__[j], k1, u, &l, h__);
	m = i__;
	i__2 = *k1;
	for (kk = 1; kk <= i__2; ++kk) {
	    d__[m] = h__[kk - 1];
	    m += *idim;
/* L300: */
	}
	j += *n;
/* L400: */
    }
L500:
    return 0;
} /* cualde_ */

