/* bispev.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int bispev_(doublereal *tx, integer *nx, doublereal *ty, 
	integer *ny, doublereal *c__, integer *kx, integer *ky, doublereal *x,
	 integer *mx, doublereal *y, integer *my, doublereal *z__, doublereal 
	*wrk, integer *lwrk, integer *iwrk, integer *kwrk, integer *ier)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, iw, lwest;
    extern /* Subroutine */ int fpbisp_(doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *);

/*  subroutine bispev evaluates on a grid (x(i),y(j)),i=1,...,mx; j=1,... */
/*  ,my a bivariate spline s(x,y) of degrees kx and ky, given in the */
/*  b-spline representation. */

/*  calling sequence: */
/*     call bispev(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wrk,lwrk, */
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
/*           lwrk >= mx*(kx+1)+my*(ky+1) */
/*   iwrk  : integer array of dimension kwrk. used as workspace. */
/*   kwrk  : integer, specifying the dimension of iwrk. kwrk >= mx+my. */

/*  output parameters: */
/*   z     : real array of dimension (mx*my). */
/*           on successful exit z(my*(i-1)+j) contains the value of s(x,y) */
/*           at the point (x(i),y(j)),i=1,...,mx;j=1,...,my. */
/*   ier   : integer error flag */
/*    ier=0 : normal return */
/*    ier=10: invalid input data (see restrictions) */

/*  restrictions: */
/*   mx >=1, my >=1, lwrk>=mx*(kx+1)+my*(ky+1), kwrk>=mx+my */
/*   tx(kx+1) <= x(i-1) <= x(i) <= tx(nx-kx), i=2,...,mx */
/*   ty(ky+1) <= y(j-1) <= y(j) <= ty(ny-ky), j=2,...,my */

/*  other subroutines required: */
/*    fpbisp,fpbspl */

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
    lwest = (*kx + 1) * *mx + (*ky + 1) * *my;
    if (*lwrk < lwest) {
	goto L100;
    }
    if (*kwrk < *mx + *my) {
	goto L100;
    }
    if (*mx < 1) {
	goto L100;
    }
    if (*mx == 1) {
	goto L30;
    }
    goto L10;
L10:
    i__1 = *mx;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (x[i__] < x[i__ - 1]) {
	    goto L100;
	}
/* L20: */
    }
L30:
    if (*my < 1) {
	goto L100;
    }
    if (*my == 1) {
	goto L60;
    }
    goto L40;
L40:
    i__1 = *my;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (y[i__] < y[i__ - 1]) {
	    goto L100;
	}
/* L50: */
    }
L60:
    *ier = 0;
    iw = *mx * (*kx + 1) + 1;
    fpbisp_(&tx[1], nx, &ty[1], ny, &c__[1], kx, ky, &x[1], mx, &y[1], my, &
	    z__[1], &wrk[1], &wrk[iw], &iwrk[1], &iwrk[*mx + 1]);
L100:
    return 0;
} /* bispev_ */

