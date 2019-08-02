/* alnrel.f -- translated by f2c (version 20190311).
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

doublereal alnrel_(doublereal *a)
{
    /* Initialized data */

    static doublereal p1 = -1.29418923021993;
    static doublereal p2 = .405303492862024;
    static doublereal p3 = -.0178874546012214;
    static doublereal q1 = -1.62752256355323;
    static doublereal q2 = .747811014037616;
    static doublereal q3 = -.0845104217945565;

    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static doublereal t, w, x, t2;

/* ----------------------------------------------------------------------- */
/*            EVALUATION OF THE FUNCTION LN(1 + A) */
/* ----------------------------------------------------------------------- */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Data statements .. */
/*     .. */
/*     .. Executable Statements .. */
/* -------------------------- */
    if (abs(*a) > .375) {
	goto L10;
    }
    t = *a / (*a + 2.);
    t2 = t * t;
    w = (((p3 * t2 + p2) * t2 + p1) * t2 + 1.) / (((q3 * t2 + q2) * t2 + q1) *
	     t2 + 1.);
    ret_val = t * 2. * w;
    return ret_val;

L10:
    x = *a + 1.;
    ret_val = log(x);
    return ret_val;
} /* alnrel_ */

