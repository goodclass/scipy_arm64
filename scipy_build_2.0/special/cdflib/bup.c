/* bup.f -- translated by f2c (version 20190311).
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
static integer c__0 = 0;

doublereal bup_(doublereal *a, doublereal *b, doublereal *x, doublereal *y, 
	integer *n, doublereal *eps)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static doublereal d__;
    static integer i__, k;
    static doublereal l, r__, t, w;
    static integer mu;
    static doublereal ap1;
    static integer kp1, nm1;
    static doublereal apb;
    extern doublereal brcmp1_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), exparg_(integer *);

/* ----------------------------------------------------------------------- */
/*     EVALUATION OF IX(A,B) - IX(A+N,B) WHERE N IS A POSITIVE INTEGER. */
/*     EPS IS THE TOLERANCE USED. */
/* ----------------------------------------------------------------------- */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*          OBTAIN THE SCALING FACTOR EXP(-MU) AND */
/*             EXP(MU)*(X**A*Y**B/BETA(A,B))/A */

    apb = *a + *b;
    ap1 = *a + 1.;
    mu = 0;
    d__ = 1.;
    if (*n == 1 || *a < 1.) {
	goto L10;
    }
    if (apb < ap1 * 1.1) {
	goto L10;
    }
    mu = (d__1 = exparg_(&c__1), (integer) abs(d__1));
    k = (integer) exparg_(&c__0);
    if (k < mu) {
	mu = k;
    }
    t = (doublereal) mu;
    d__ = exp(-t);

L10:
    ret_val = brcmp1_(&mu, a, b, x, y) / *a;
    if (*n == 1 || ret_val == 0.) {
	return ret_val;
    }
    nm1 = *n - 1;
    w = d__;

/*          LET K BE THE INDEX OF THE MAXIMUM TERM */

    k = 0;
    if (*b <= 1.) {
	goto L50;
    }
    if (*y > 1e-4) {
	goto L20;
    }
    k = nm1;
    goto L30;
L20:
    r__ = (*b - 1.) * *x / *y - *a;
    if (r__ < 1.) {
	goto L50;
    }
    k = nm1;
    t = (doublereal) nm1;
    if (r__ < t) {
	k = (integer) r__;
    }

/*          ADD THE INCREASING TERMS OF THE SERIES */

L30:
    i__1 = k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l = (doublereal) (i__ - 1);
	d__ = (apb + l) / (ap1 + l) * *x * d__;
	w += d__;
/* L40: */
    }
    if (k == nm1) {
	goto L70;
    }

/*          ADD THE REMAINING TERMS OF THE SERIES */

L50:
    kp1 = k + 1;
    i__1 = nm1;
    for (i__ = kp1; i__ <= i__1; ++i__) {
	l = (doublereal) (i__ - 1);
	d__ = (apb + l) / (ap1 + l) * *x * d__;
	w += d__;
	if (d__ <= *eps * w) {
	    goto L70;
	}
/* L60: */
    }

/*               TERMINATE THE PROCEDURE */

L70:
    ret_val *= w;
    return ret_val;
} /* bup_ */

