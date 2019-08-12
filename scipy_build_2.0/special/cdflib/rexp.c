/* rexp.f -- translated by f2c (version 20190311).
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

doublereal rexp_(doublereal *x)
{
    /* Initialized data */

    static doublereal p1 = 9.14041914819518e-10;
    static doublereal p2 = .0238082361044469;
    static doublereal q1 = -.499999999085958;
    static doublereal q2 = .107141568980644;
    static doublereal q3 = -.0119041179760821;
    static doublereal q4 = 5.95130811860248e-4;

    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static doublereal w;

/* ----------------------------------------------------------------------- */
/*            EVALUATION OF THE FUNCTION EXP(X) - 1 */
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
/* ----------------------- */
    if (abs(*x) > .15) {
	goto L10;
    }
    ret_val = *x * (((p2 * *x + p1) * *x + 1.) / ((((q4 * *x + q3) * *x + q2) 
	    * *x + q1) * *x + 1.));
    return ret_val;

L10:
    w = exp(*x);
    if (*x > 0.) {
	goto L20;
    }
    ret_val = w - .5 - .5;
    return ret_val;
L20:
    ret_val = w * (.5 - 1. / w + .5);
    return ret_val;
} /* rexp_ */

