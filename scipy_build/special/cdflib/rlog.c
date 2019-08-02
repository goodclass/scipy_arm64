/* rlog.f -- translated by f2c (version 20190311).
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

doublereal rlog_(doublereal *x)
{
    /* Initialized data */

    static doublereal a = .0566749439387324;
    static doublereal b = .0456512608815524;
    static doublereal p0 = .333333333333333;
    static doublereal p1 = -.224696413112536;
    static doublereal p2 = .00620886815375787;
    static doublereal q1 = -1.27408923933623;
    static doublereal q2 = .354508718369557;

    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static doublereal r__, t, u, w, w1;

/*     ------------------- */
/*     COMPUTATION OF  X - 1 - LN(X) */
/*     ------------------- */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Data statements .. */
/*     ------------------- */
/*     .. */
/*     .. Executable Statements .. */
/*     ------------------- */
    if (*x < .61 || *x > 1.57) {
	goto L40;
    }
    if (*x < .82) {
	goto L10;
    }
    if (*x > 1.18) {
	goto L20;
    }

/*              ARGUMENT REDUCTION */

    u = *x - .5 - .5;
    w1 = 0.;
    goto L30;

L10:
    u = *x - .7;
    u /= .7;
    w1 = a - u * .3;
    goto L30;

L20:
    u = *x * .75 - 1.;
    w1 = b + u / 3.;

/*               SERIES EXPANSION */

L30:
    r__ = u / (u + 2.);
    t = r__ * r__;
    w = ((p2 * t + p1) * t + p0) / ((q2 * t + q1) * t + 1.);
    ret_val = t * 2. * (1. / (1. - r__) - r__ * w) + w1;
    return ret_val;


L40:
    r__ = *x - .5 - .5;
    ret_val = r__ - log(*x);
    return ret_val;
} /* rlog_ */

