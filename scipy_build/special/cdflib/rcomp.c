/* rcomp.f -- translated by f2c (version 20190311).
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

doublereal rcomp_(doublereal *a, doublereal *x)
{
    /* Initialized data */

    static doublereal rt2pin = .398942280401433;

    /* System generated locals */
    doublereal ret_val, d__1;

    /* Builtin functions */
    double log(doublereal), exp(doublereal), sqrt(doublereal);

    /* Local variables */
    static doublereal t, u, t1;
    extern doublereal gam1_(doublereal *), rlog_(doublereal *), gamma_(
	    doublereal *);

/*     ------------------- */
/*     EVALUATION OF EXP(-X)*X**A/GAMMA(A) */
/*     ------------------- */
/*     RT2PIN = 1/SQRT(2*PI) */
/*     ------------------- */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Data statements .. */
/*     .. */
/*     .. Executable Statements .. */
/*     ------------------- */
    ret_val = 0.;
    if (*a >= 20.) {
	goto L20;
    }
    t = *a * log(*x) - *x;
    if (*a >= 1.) {
	goto L10;
    }
    ret_val = *a * exp(t) * (gam1_(a) + 1.);
    return ret_val;
L10:
    ret_val = exp(t) / gamma_(a);
    return ret_val;

L20:
    u = *x / *a;
    if (u == 0.) {
	return ret_val;
    }
/* Computing 2nd power */
    d__1 = 1. / *a;
    t = d__1 * d__1;
    t1 = (((t * .75 - 1.) * t + 3.5) * t - 105.) / (*a * 1260.);
    t1 -= *a * rlog_(&u);
    ret_val = rt2pin * sqrt(*a) * exp(t1);
    return ret_val;
} /* rcomp_ */

