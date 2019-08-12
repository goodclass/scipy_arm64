/* apser.f -- translated by f2c (version 20190311).
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

doublereal apser_(doublereal *a, doublereal *b, doublereal *x, doublereal *
	eps)
{
    /* Initialized data */

    static doublereal g = .577215664901533;

    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static doublereal c__, j, s, t, aj, bx;
    extern doublereal psi_(doublereal *);
    static doublereal tol;

/* ----------------------------------------------------------------------- */
/*     APSER YIELDS THE INCOMPLETE BETA RATIO I(SUB(1-X))(B,A) FOR */
/*     A .LE. MIN(EPS,EPS*B), B*X .LE. 1, AND X .LE. 0.5. USED WHEN */
/*     A IS VERY SMALL. USE ONLY IF ABOVE INEQUALITIES ARE SATISFIED. */
/* ----------------------------------------------------------------------- */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Data statements .. */
/* -------------------- */
/*     .. */
/*     .. Executable Statements .. */
/* -------------------- */
    bx = *b * *x;
    t = *x - bx;
    if (*b * *eps > .02) {
	goto L10;
    }
    c__ = log(*x) + psi_(b) + g + t;
    goto L20;
L10:
    c__ = log(bx) + g + t;

L20:
    tol = *eps * 5. * abs(c__);
    j = 1.;
    s = 0.;
L30:
    j += 1.;
    t *= *x - bx / j;
    aj = t / j;
    s += aj;
    if (abs(aj) > tol) {
	goto L30;
    }

    ret_val = -(*a) * (c__ + s);
    return ret_val;
} /* apser_ */

