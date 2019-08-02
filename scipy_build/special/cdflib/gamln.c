/* gamln.f -- translated by f2c (version 20190311).
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

doublereal gamln_(doublereal *a)
{
    /* Initialized data */

    static doublereal d__ = .418938533204673;
    static doublereal c0 = .0833333333333333;
    static doublereal c1 = -.00277777777760991;
    static doublereal c2 = 7.9365066682539e-4;
    static doublereal c3 = -5.9520293135187e-4;
    static doublereal c4 = 8.37308034031215e-4;
    static doublereal c5 = -.00165322962780713;

    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static integer i__, n;
    static doublereal t, w;
    extern doublereal gamln1_(doublereal *);

/* ----------------------------------------------------------------------- */
/*            EVALUATION OF LN(GAMMA(A)) FOR POSITIVE A */
/* ----------------------------------------------------------------------- */
/*     WRITTEN BY ALFRED H. MORRIS */
/*          NAVAL SURFACE WARFARE CENTER */
/*          DAHLGREN, VIRGINIA */
/* -------------------------- */
/*     D = 0.5*(LN(2*PI) - 1) */
/* -------------------------- */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Data statements .. */
/* -------------------------- */
/*     .. */
/*     .. Executable Statements .. */
/* ----------------------------------------------------------------------- */
    if (*a > .8) {
	goto L10;
    }
    ret_val = gamln1_(a) - log(*a);
    return ret_val;
L10:
    if (*a > 2.25) {
	goto L20;
    }
    t = *a - .5 - .5;
    ret_val = gamln1_(&t);
    return ret_val;

L20:
    if (*a >= 10.) {
	goto L40;
    }
    n = (integer) (*a - 1.25);
    t = *a;
    w = 1.;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	t += -1.;
	w = t * w;
/* L30: */
    }
    d__1 = t - 1.;
    ret_val = gamln1_(&d__1) + log(w);
    return ret_val;

L40:
/* Computing 2nd power */
    d__1 = 1. / *a;
    t = d__1 * d__1;
    w = (((((c5 * t + c4) * t + c3) * t + c2) * t + c1) * t + c0) / *a;
    ret_val = d__ + w + (*a - .5) * (log(*a) - 1.);
    return ret_val;
} /* gamln_ */

