/* parder.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int parder_(doublereal *tx, integer *nx, doublereal *ty, 
	integer *ny, doublereal *c__, integer *kx, integer *ky, integer *nux, 
	integer *nuy, doublereal *x, integer *mx, doublereal *y, integer *my, 
	doublereal *z__, doublereal *wrk, integer *lwrk, integer *iwrk, 
	integer *kwrk, integer *ier)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, m, l1, l2, m0, m1;
    static doublereal ak;
    static integer nc, lx, ly, kx1, ky1;
    static doublereal fac;
    static integer kkx, kky, iwx, iwy, nxx, nyy, nkx1, nky1, lwest;
    extern /* Subroutine */ int fpbisp_(doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *);

/*  subroutine parder evaluates on a grid (x(i),y(j)),i=1,...,mx; j=1,... */
/*  ,my the partial derivative ( order nux,nuy) of a bivariate spline */
/*  s(x,y) of degrees kx and ky, given in the b-spline representation. */

/*  calling sequence: */
/*     call parder(tx,nx,ty,ny,c,kx,ky,nux,nuy,x,mx,y,my,z,wrk,lwrk, */
/*    * iwrk,kwrk,ier) */

/*  input parameters: */
/*   tx    : real array, length nx, which contains the position of the */
/*           knots in the x-direction. */
/*   nx    : integer, giving the total number of knots in the x-direction */
/*   ty    : real array, length ny, which contains the position of the */
/*           knots in the y-direction. */
/*   ny    : integer, giving the total number of knots in the y-direction */
/*   c     : real array, length (nx-kx-1)*(ny-ky-1), which contains the */
/*           b-spline coefficients. */
/*   kx,ky : integer values, giving the degrees of the spline. */
/*   nux   : integer values, specifying the order of the partial */
/*   nuy     derivative. 0<=nux<kx, 0<=nuy<ky. */
/*   x     : real array of dimension (mx). */
/*           before entry x(i) must be set to the x co-ordinate of the */
/*           i-th grid point along the x-axis. */
/*           tx(kx+1)<=x(i-1)<=x(i)<=tx(nx-kx), i=2,...,mx. */
/*   mx    : on entry mx must specify the number of grid points along */
/*           the x-axis. mx >=1. */
/*   y     : real array of dimension (my). */
/*           before entry y(j) must be set to the y co-ordinate of the */
/*           j-th grid point along the y-axis. */
/*           ty(ky+1)<=y(j-1)<=y(j)<=ty(ny-ky), j=2,...,my. */
/*   my    : on entry my must specify the number of grid points along */
/*           the y-axis. my >=1. */
/*   wrk   : real array of dimension lwrk. used as workspace. */
/*   lwrk  : integer, specifying the dimension of wrk. */
/*           lwrk >= mx*(kx+1-nux)+my*(ky+1-nuy)+(nx-kx-1)*(ny-ky-1) */
/*   iwrk  : integer array of dimension kwrk. used as workspace. */
/*   kwrk  : integer, specifying the dimension of iwrk. kwrk >= mx+my. */

/*  output parameters: */
/*   z     : real array of dimension (mx*my). */
/*           on successful exit z(my*(i-1)+j) contains the value of the */
/*           specified partial derivative of s(x,y) at the point */
/*           (x(i),y(j)),i=1,...,mx;j=1,...,my. */
/*   ier   : integer error flag */
/*    ier=0 : normal return */
/*    ier=10: invalid input data (see restrictions) */

/*  restrictions: */
/*   mx >=1, my >=1, 0 <= nux < kx, 0 <= nuy < ky, kwrk>=mx+my */
/*   lwrk>=mx*(kx+1-nux)+my*(ky+1-nuy)+(nx-kx-1)*(ny-ky-1), */
/*   tx(kx+1) <= x(i-1) <= x(i) <= tx(nx-kx), i=2,...,mx */
/*   ty(ky+1) <= y(j-1) <= y(j) <= ty(ny-ky), j=2,...,my */

/*  other subroutines required: */
/*    fpbisp,fpbspl */

/*  references : */
/*    de boor c : on calculating with b-splines, j. approximation theory */
/*                6 (1972) 50-62. */
/*   dierckx p. : curve and surface fitting with splines, monographs on */
/*                numerical analysis, oxford university press, 1993. */

/*  author : */
/*    p.dierckx */
/*    dept. computer science, k.u.leuven */
/*    celestijnenlaan 200a, b-3001 heverlee, belgium. */
/*    e-mail : Paul.Dierckx@cs.kuleuven.ac.be */

/*  latest update : march 1989 */

/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  .. */
/*  before starting computations a data check is made. if the input data */
/*  are invalid control is immediately repassed to the calling program. */
    /* Parameter adjustments */
    --tx;
    --ty;
    --c__;
    --x;
    --z__;
    --y;
    --wrk;
    --iwrk;

    /* Function Body */
    *ier = 10;
    kx1 = *kx + 1;
    ky1 = *ky + 1;
    nkx1 = *nx - kx1;
    nky1 = *ny - ky1;
    nc = nkx1 * nky1;
    if (*nux < 0 || *nux >= *kx) {
	goto L400;
    }
    if (*nuy < 0 || *nuy >= *ky) {
	goto L400;
    }
    lwest = nc + (kx1 - *nux) * *mx + (ky1 - *nuy) * *my;
    if (*lwrk < lwest) {
	goto L400;
    }
    if (*kwrk < *mx + *my) {
	goto L400;
    }
    if (*mx < 1) {
	goto L400;
    }
    if (*mx == 1) {
	goto L30;
    }
    goto L10;
L10:
    i__1 = *mx;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (x[i__] < x[i__ - 1]) {
	    goto L400;
	}
/* L20: */
    }
L30:
    if (*my < 1) {
	goto L400;
    }
    if (*my == 1) {
	goto L60;
    }
    goto L40;
L40:
    i__1 = *my;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (y[i__] < y[i__ - 1]) {
	    goto L400;
	}
/* L50: */
    }
L60:
    *ier = 0;
    nxx = nkx1;
    nyy = nky1;
    kkx = *kx;
    kky = *ky;
/*  the partial derivative of order (nux,nuy) of a bivariate spline of */
/*  degrees kx,ky is a bivariate spline of degrees kx-nux,ky-nuy. */
/*  we calculate the b-spline coefficients of this spline */
    i__1 = nc;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wrk[i__] = c__[i__];
/* L70: */
    }
    if (*nux == 0) {
	goto L200;
    }
    lx = 1;
    i__1 = *nux;
    for (j = 1; j <= i__1; ++j) {
	ak = (doublereal) kkx;
	--nxx;
	l1 = lx;
	m0 = 1;
	i__2 = nxx;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ++l1;
	    l2 = l1 + kkx;
	    fac = tx[l2] - tx[l1];
	    if (fac <= 0.f) {
		goto L90;
	    }
	    i__3 = nyy;
	    for (m = 1; m <= i__3; ++m) {
		m1 = m0 + nyy;
		wrk[m0] = (wrk[m1] - wrk[m0]) * ak / fac;
		++m0;
/* L80: */
	    }
L90:
	    ;
	}
	++lx;
	--kkx;
/* L100: */
    }
L200:
    if (*nuy == 0) {
	goto L300;
    }
    ly = 1;
    i__1 = *nuy;
    for (j = 1; j <= i__1; ++j) {
	ak = (doublereal) kky;
	--nyy;
	l1 = ly;
	i__2 = nyy;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ++l1;
	    l2 = l1 + kky;
	    fac = ty[l2] - ty[l1];
	    if (fac <= 0.f) {
		goto L220;
	    }
	    m0 = i__;
	    i__3 = nxx;
	    for (m = 1; m <= i__3; ++m) {
		m1 = m0 + 1;
		wrk[m0] = (wrk[m1] - wrk[m0]) * ak / fac;
		m0 += nky1;
/* L210: */
	    }
L220:
	    ;
	}
	++ly;
	--kky;
/* L230: */
    }
    m0 = nyy;
    m1 = nky1;
    i__1 = nxx;
    for (m = 2; m <= i__1; ++m) {
	i__2 = nyy;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ++m0;
	    ++m1;
	    wrk[m0] = wrk[m1];
/* L240: */
	}
	m1 += *nuy;
/* L250: */
    }
/*  we partition the working space and evaluate the partial derivative */
L300:
    iwx = nxx * nyy + 1;
    iwy = iwx + *mx * (kx1 - *nux);
    i__1 = *nx - (*nux << 1);
    i__2 = *ny - (*nuy << 1);
    fpbisp_(&tx[*nux + 1], &i__1, &ty[*nuy + 1], &i__2, &wrk[1], &kkx, &kky, &
	    x[1], mx, &y[1], my, &z__[1], &wrk[iwx], &wrk[iwy], &iwrk[1], &
	    iwrk[*mx + 1]);
L400:
    return 0;
} /* parder_ */

