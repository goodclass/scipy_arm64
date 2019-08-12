/* bispeu.f -- translated by f2c (version 20190311).
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

/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int bispeu_(doublereal *tx, integer *nx, doublereal *ty, 
	integer *ny, doublereal *c__, integer *kx, integer *ky, doublereal *x,
	 doublereal *y, doublereal *z__, integer *m, doublereal *wrk, integer 
	*lwrk, integer *ier)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, iwrk[2], lwest;
    extern /* Subroutine */ int fpbisp_(doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *);

/*  subroutine bispeu evaluates on a set of points (x(i),y(i)),i=1,...,m */
/*  a bivariate spline s(x,y) of degrees kx and ky, given in the */
/*  b-spline representation. */

/*  calling sequence: */
/*     call bispeu(tx,nx,ty,ny,c,kx,ky,x,y,z,m,wrk,lwrk, */
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
/*   y     : real array of dimension (my). */
/*   m     : on entry m must specify the number points. m >= 1. */
/*   wrk   : real array of dimension lwrk. used as workspace. */
/*   lwrk  : integer, specifying the dimension of wrk. */
/*           lwrk >= kx+ky+2 */

/*  output parameters: */
/*   z     : real array of dimension m. */
/*           on successful exit z(i) contains the value of s(x,y) */
/*           at the point (x(i),y(i)), i=1,...,m. */
/*   ier   : integer error flag */
/*    ier=0 : normal return */
/*    ier=10: invalid input data (see restrictions) */

/*  restrictions: */
/*   m >=1, lwrk>=mx*(kx+1)+my*(ky+1), kwrk>=mx+my */
/*   tx(kx+1) <= x(i-1) <= x(i) <= tx(nx-kx), i=2,...,mx */
/*   ty(ky+1) <= y(j-1) <= y(j) <= ty(ny-ky), j=2,...,my */

/*  other subroutines required: */
/*    fpbisp,fpbspl */

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
    --z__;
    --y;
    --x;
    --wrk;

    /* Function Body */
    *ier = 10;
    lwest = *kx + *ky + 2;
    if (*lwrk < lwest) {
	goto L100;
    }
    if (*m < 1) {
	goto L100;
    }
    *ier = 0;
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	fpbisp_(&tx[1], nx, &ty[1], ny, &c__[1], kx, ky, &x[i__], &c__1, &y[
		i__], &c__1, &z__[i__], &wrk[1], &wrk[*kx + 2], iwrk, &iwrk[1]
		);
/* L10: */
    }
L100:
    return 0;
} /* bispeu_ */

