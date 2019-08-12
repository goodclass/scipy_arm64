/* profil.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int profil_(integer *iopt, doublereal *tx, integer *nx, 
	doublereal *ty, integer *ny, doublereal *c__, integer *kx, integer *
	ky, doublereal *u, integer *nu, doublereal *cu, integer *ier)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal h__[6];
    static integer i__, j, l, m, l1, m0, kx1, ky1;
    static doublereal sum;
    static integer nkx1, nky1;
    extern /* Subroutine */ int fpbspl_(doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *);

/*  if iopt=0 subroutine profil calculates the b-spline coefficients of */
/*  the univariate spline f(y) = s(u,y) with s(x,y) a bivariate spline of */
/*  degrees kx and ky, given in the b-spline representation. */
/*  if iopt = 1 it calculates the b-spline coefficients of the univariate */
/*  spline g(x) = s(x,u) */

/*  calling sequence: */
/*     call profil(iopt,tx,nx,ty,ny,c,kx,ky,u,nu,cu,ier) */

/*  input parameters: */
/*   iopt  : integer flag, specifying whether the profile f(y) (iopt=0) */
/*           or the profile g(x) (iopt=1) must be determined. */
/*   tx    : real array, length nx, which contains the position of the */
/*           knots in the x-direction. */
/*   nx    : integer, giving the total number of knots in the x-direction */
/*   ty    : real array, length ny, which contains the position of the */
/*           knots in the y-direction. */
/*   ny    : integer, giving the total number of knots in the y-direction */
/*   c     : real array, length (nx-kx-1)*(ny-ky-1), which contains the */
/*           b-spline coefficients. */
/*   kx,ky : integer values, giving the degrees of the spline. */
/*   u     : real value, specifying the requested profile. */
/*           tx(kx+1)<=u<=tx(nx-kx), if iopt=0. */
/*           ty(ky+1)<=u<=ty(ny-ky), if iopt=1. */
/*   nu    : on entry nu must specify the dimension of the array cu. */
/*           nu >= ny if iopt=0, nu >= nx if iopt=1. */

/*  output parameters: */
/*   cu    : real array of dimension (nu). */
/*           on successful exit this array contains the b-spline */
/*   ier   : integer error flag */
/*    ier=0 : normal return */
/*    ier=10: invalid input data (see restrictions) */

/*  restrictions: */
/*   if iopt=0 : tx(kx+1) <= u <= tx(nx-kx), nu >=ny. */
/*   if iopt=1 : ty(ky+1) <= u <= ty(ny-ky), nu >=nx. */

/*  other subroutines required: */
/*    fpbspl */

/*  author : */
/*    p.dierckx */
/*    dept. computer science, k.u.leuven */
/*    celestijnenlaan 200a, b-3001 heverlee, belgium. */
/*    e-mail : Paul.Dierckx@cs.kuleuven.ac.be */

/*  latest update : march 1987 */

/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  ..local array */
/*  .. */
/*  before starting computations a data check is made. if the input data */
/*  are invalid control is immediately repassed to the calling program. */
    /* Parameter adjustments */
    --tx;
    --ty;
    --c__;
    --cu;

    /* Function Body */
    kx1 = *kx + 1;
    ky1 = *ky + 1;
    nkx1 = *nx - kx1;
    nky1 = *ny - ky1;
    *ier = 10;
    if (*iopt != 0) {
	goto L200;
    }
    if (*nu < *ny) {
	goto L300;
    }
    if (*u < tx[kx1] || *u > tx[nkx1 + 1]) {
	goto L300;
    }
/*  the b-splinecoefficients of f(y) = s(u,y). */
    *ier = 0;
    l = kx1;
    l1 = l + 1;
L110:
    if (*u < tx[l1] || l == nkx1) {
	goto L120;
    }
    l = l1;
    l1 = l + 1;
    goto L110;
L120:
    fpbspl_(&tx[1], nx, kx, u, &l, h__);
    m0 = (l - kx1) * nky1 + 1;
    i__1 = nky1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	m = m0;
	sum = 0.f;
	i__2 = kx1;
	for (j = 1; j <= i__2; ++j) {
	    sum += h__[j - 1] * c__[m];
	    m += nky1;
/* L130: */
	}
	cu[i__] = sum;
	++m0;
/* L140: */
    }
    goto L300;
L200:
    if (*nu < *nx) {
	goto L300;
    }
    if (*u < ty[ky1] || *u > ty[nky1 + 1]) {
	goto L300;
    }
/*  the b-splinecoefficients of g(x) = s(x,u). */
    *ier = 0;
    l = ky1;
    l1 = l + 1;
L210:
    if (*u < ty[l1] || l == nky1) {
	goto L220;
    }
    l = l1;
    l1 = l + 1;
    goto L210;
L220:
    fpbspl_(&ty[1], ny, ky, u, &l, h__);
    m0 = l - *ky;
    i__1 = nkx1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	m = m0;
	sum = 0.f;
	i__2 = ky1;
	for (j = 1; j <= i__2; ++j) {
	    sum += h__[j - 1] * c__[m];
	    ++m;
/* L230: */
	}
	cu[i__] = sum;
	m0 += nky1;
/* L240: */
    }
L300:
    return 0;
} /* profil_ */

