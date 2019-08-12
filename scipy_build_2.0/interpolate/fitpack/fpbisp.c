/* fpbisp.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int fpbisp_(doublereal *tx, integer *nx, doublereal *ty, 
	integer *ny, doublereal *c__, integer *kx, integer *ky, doublereal *x,
	 integer *mx, doublereal *y, integer *my, doublereal *z__, doublereal 
	*wx, doublereal *wy, integer *lx, integer *ly)
{
    /* System generated locals */
    integer wx_dim1, wx_offset, wy_dim1, wy_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static doublereal h__[6];
    static integer i__, j, l, m, i1, j1, l1, l2;
    static doublereal tb, te, sp;
    static integer kx1, ky1;
    static doublereal arg;
    static integer nkx1, nky1;
    extern /* Subroutine */ int fpbspl_(doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *);

/*  ..scalar arguments.. */
/*  ..array arguments.. */
/*  ..local scalars.. */
/*  ..local arrays.. */
/*  ..subroutine references.. */
/*    fpbspl */
/*  .. */
    /* Parameter adjustments */
    --tx;
    --ty;
    --c__;
    --lx;
    wx_dim1 = *mx;
    wx_offset = 1 + wx_dim1;
    wx -= wx_offset;
    --x;
    --ly;
    wy_dim1 = *my;
    wy_offset = 1 + wy_dim1;
    wy -= wy_offset;
    --z__;
    --y;

    /* Function Body */
    kx1 = *kx + 1;
    nkx1 = *nx - kx1;
    tb = tx[kx1];
    te = tx[nkx1 + 1];
    l = kx1;
    l1 = l + 1;
    i__1 = *mx;
    for (i__ = 1; i__ <= i__1; ++i__) {
	arg = x[i__];
	if (arg < tb) {
	    arg = tb;
	}
	if (arg > te) {
	    arg = te;
	}
L10:
	if (arg < tx[l1] || l == nkx1) {
	    goto L20;
	}
	l = l1;
	l1 = l + 1;
	goto L10;
L20:
	fpbspl_(&tx[1], nx, kx, &arg, &l, h__);
	lx[i__] = l - kx1;
	i__2 = kx1;
	for (j = 1; j <= i__2; ++j) {
	    wx[i__ + j * wx_dim1] = h__[j - 1];
/* L30: */
	}
/* L40: */
    }
    ky1 = *ky + 1;
    nky1 = *ny - ky1;
    tb = ty[ky1];
    te = ty[nky1 + 1];
    l = ky1;
    l1 = l + 1;
    i__1 = *my;
    for (i__ = 1; i__ <= i__1; ++i__) {
	arg = y[i__];
	if (arg < tb) {
	    arg = tb;
	}
	if (arg > te) {
	    arg = te;
	}
L50:
	if (arg < ty[l1] || l == nky1) {
	    goto L60;
	}
	l = l1;
	l1 = l + 1;
	goto L50;
L60:
	fpbspl_(&ty[1], ny, ky, &arg, &l, h__);
	ly[i__] = l - ky1;
	i__2 = ky1;
	for (j = 1; j <= i__2; ++j) {
	    wy[i__ + j * wy_dim1] = h__[j - 1];
/* L70: */
	}
/* L80: */
    }
    m = 0;
    i__1 = *mx;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l = lx[i__] * nky1;
	i__2 = kx1;
	for (i1 = 1; i1 <= i__2; ++i1) {
	    h__[i1 - 1] = wx[i__ + i1 * wx_dim1];
/* L90: */
	}
	i__2 = *my;
	for (j = 1; j <= i__2; ++j) {
	    l1 = l + ly[j];
	    sp = 0.f;
	    i__3 = kx1;
	    for (i1 = 1; i1 <= i__3; ++i1) {
		l2 = l1;
		i__4 = ky1;
		for (j1 = 1; j1 <= i__4; ++j1) {
		    ++l2;
		    sp += c__[l2] * h__[i1 - 1] * wy[j + j1 * wy_dim1];
/* L100: */
		}
		l1 += nky1;
/* L110: */
	    }
	    ++m;
	    z__[m] = sp;
/* L120: */
	}
/* L130: */
    }
    return 0;
} /* fpbisp_ */

