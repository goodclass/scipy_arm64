/* fpser.f -- translated by f2c (version 20190311).
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

doublereal fpser_(doublereal *a, doublereal *b, doublereal *x, doublereal *
	eps)
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double log(doublereal), exp(doublereal);

    /* Local variables */
    static doublereal c__, s, t, an, tol;
    extern doublereal exparg_(integer *);

/* ----------------------------------------------------------------------- */

/*                 EVALUATION OF I (A,B) */
/*                                X */

/*          FOR B .LT. MIN(EPS,EPS*A) AND X .LE. 0.5. */

/* ----------------------------------------------------------------------- */

/*                  SET  FPSER = X**A */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */
    ret_val = 1.;
    if (*a <= *eps * .001) {
	goto L10;
    }
    ret_val = 0.;
    t = *a * log(*x);
    if (t < exparg_(&c__1)) {
	return ret_val;
    }
    ret_val = exp(t);

/*                NOTE THAT 1/B(A,B) = B */

L10:
    ret_val = *b / *a * ret_val;
    tol = *eps / *a;
    an = *a + 1.;
    t = *x;
    s = t / an;
L20:
    an += 1.;
    t = *x * t;
    c__ = t / an;
    s += c__;
    if (abs(c__) > tol) {
	goto L20;
    }

    ret_val *= *a * s + 1.;
    return ret_val;
} /* fpser_ */

