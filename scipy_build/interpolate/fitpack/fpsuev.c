/* fpsuev.f -- translated by f2c (version 20190311).
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

static integer c__3 = 3;

/* Subroutine */ int fpsuev_(integer *idim, doublereal *tu, integer *nu, 
	doublereal *tv, integer *nv, doublereal *c__, doublereal *u, integer *
	mu, doublereal *v, integer *mv, doublereal *f, doublereal *wu, 
	doublereal *wv, integer *lu, integer *lv)
{
    /* System generated locals */
    integer wu_dim1, wu_offset, wv_dim1, wv_offset, i__1, i__2, i__3;

    /* Local variables */
    static doublereal h__[4];
    static integer i__, j, k, l, m, i1, j1, l1, l2, l3;
    static doublereal tb, te, sp;
    static integer nu4, nv4;
    static doublereal arg;
    static integer nuv;
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
    --tu;
    --c__;
    --tv;
    --lu;
    wu_dim1 = *mu;
    wu_offset = 1 + wu_dim1;
    wu -= wu_offset;
    --u;
    --lv;
    wv_dim1 = *mv;
    wv_offset = 1 + wv_dim1;
    wv -= wv_offset;
    --f;
    --v;

    /* Function Body */
    nu4 = *nu - 4;
    tb = tu[4];
    te = tu[nu4 + 1];
    l = 4;
    l1 = l + 1;
    i__1 = *mu;
    for (i__ = 1; i__ <= i__1; ++i__) {
	arg = u[i__];
	if (arg < tb) {
	    arg = tb;
	}
	if (arg > te) {
	    arg = te;
	}
L10:
	if (arg < tu[l1] || l == nu4) {
	    goto L20;
	}
	l = l1;
	l1 = l + 1;
	goto L10;
L20:
	fpbspl_(&tu[1], nu, &c__3, &arg, &l, h__);
	lu[i__] = l - 4;
	for (j = 1; j <= 4; ++j) {
	    wu[i__ + j * wu_dim1] = h__[j - 1];
/* L30: */
	}
/* L40: */
    }
    nv4 = *nv - 4;
    tb = tv[4];
    te = tv[nv4 + 1];
    l = 4;
    l1 = l + 1;
    i__1 = *mv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	arg = v[i__];
	if (arg < tb) {
	    arg = tb;
	}
	if (arg > te) {
	    arg = te;
	}
L50:
	if (arg < tv[l1] || l == nv4) {
	    goto L60;
	}
	l = l1;
	l1 = l + 1;
	goto L50;
L60:
	fpbspl_(&tv[1], nv, &c__3, &arg, &l, h__);
	lv[i__] = l - 4;
	for (j = 1; j <= 4; ++j) {
	    wv[i__ + j * wv_dim1] = h__[j - 1];
/* L70: */
	}
/* L80: */
    }
    m = 0;
    nuv = nu4 * nv4;
    i__1 = *idim;
    for (k = 1; k <= i__1; ++k) {
	l3 = (k - 1) * nuv;
	i__2 = *mu;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    l = lu[i__] * nv4 + l3;
	    for (i1 = 1; i1 <= 4; ++i1) {
		h__[i1 - 1] = wu[i__ + i1 * wu_dim1];
/* L90: */
	    }
	    i__3 = *mv;
	    for (j = 1; j <= i__3; ++j) {
		l1 = l + lv[j];
		sp = 0.f;
		for (i1 = 1; i1 <= 4; ++i1) {
		    l2 = l1;
		    for (j1 = 1; j1 <= 4; ++j1) {
			++l2;
			sp += c__[l2] * h__[i1 - 1] * wv[j + j1 * wv_dim1];
/* L100: */
		    }
		    l1 += nv4;
/* L110: */
		}
		++m;
		f[m] = sp;
/* L120: */
	    }
/* L130: */
	}
/* L140: */
    }
    return 0;
} /* fpsuev_ */

