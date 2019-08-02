/* gamma_fort.f -- translated by f2c (version 20190311).
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
static integer c__0 = 0;

doublereal gamma_(doublereal *a)
{
    /* Initialized data */

    static doublereal pi = 3.1415926535898;
    static doublereal d__ = .41893853320467274178;
    static doublereal p[7] = { 5.39637273585445e-4,.0026193926004269,
	    .020449366759492,.0730981088720487,.279648642639792,
	    .553413866010467,1. };
    static doublereal q[7] = { -8.32979206704073e-4,.00470059485860584,
	    .022521113103534,-.17045896931336,-.056790276197494,
	    1.13062953091122,1. };
    static doublereal r1 = 8.20756370353826e-4;
    static doublereal r2 = -5.95156336428591e-4;
    static doublereal r3 = 7.93650663183693e-4;
    static doublereal r4 = -.00277777777770481;
    static doublereal r5 = .0833333333333333;

    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Builtin functions */
    double sin(doublereal), log(doublereal), exp(doublereal);

    /* Local variables */
    static doublereal g;
    static integer i__, j, m, n;
    static doublereal s, t, w, x, z__, bot, lnx, top;
    extern doublereal exparg_(integer *), spmpar_(integer *);

/* ----------------------------------------------------------------------- */

/*         EVALUATION OF THE GAMMA FUNCTION FOR REAL ARGUMENTS */

/*                           ----------- */

/*     GAMMA(A) IS ASSIGNED THE VALUE 0 WHEN THE GAMMA FUNCTION CANNOT */
/*     BE COMPUTED. */

/* ----------------------------------------------------------------------- */
/*     WRITTEN BY ALFRED H. MORRIS, JR. */
/*          NAVAL SURFACE WEAPONS CENTER */
/*          DAHLGREN, VIRGINIA */
/* ----------------------------------------------------------------------- */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Data statements .. */
/* -------------------------- */
/*     D = 0.5*(LN(2*PI) - 1) */
/* -------------------------- */
/* -------------------------- */
/* -------------------------- */
/*     .. */
/*     .. Executable Statements .. */
/* -------------------------- */
    ret_val = 0.;
    x = *a;
    if (abs(*a) >= 15.) {
	goto L110;
    }
/* ----------------------------------------------------------------------- */
/*            EVALUATION OF GAMMA(A) FOR ABS(A) .LT. 15 */
/* ----------------------------------------------------------------------- */
    t = 1.;
    m = (integer) (*a) - 1;

/*     LET T BE THE PRODUCT OF A-J WHEN A .GE. 2 */

    if (m < 0) {
	goto L40;
    }
    if (m == 0) {
	goto L30;
    }
    goto L10;
L10:
    i__1 = m;
    for (j = 1; j <= i__1; ++j) {
	x += -1.;
	t = x * t;
/* L20: */
    }
L30:
    x += -1.;
    goto L80;

/*     LET T BE THE PRODUCT OF A+J WHEN A .LT. 1 */

L40:
    t = *a;
    if (*a > 0.) {
	goto L70;
    }
    m = -m - 1;
    if (m == 0) {
	goto L60;
    }
    i__1 = m;
    for (j = 1; j <= i__1; ++j) {
	x += 1.;
	t = x * t;
/* L50: */
    }
L60:
    x = x + .5 + .5;
    t = x * t;
    if (t == 0.) {
	return ret_val;
    }

L70:

/*     THE FOLLOWING CODE CHECKS IF 1/T CAN OVERFLOW. THIS */
/*     CODE MAY BE OMITTED IF DESIRED. */

    if (abs(t) >= 1e-30) {
	goto L80;
    }
    if (abs(t) * spmpar_(&c__3) <= 1.0001) {
	return ret_val;
    }
    ret_val = 1. / t;
    return ret_val;

/*     COMPUTE GAMMA(1 + X) FOR  0 .LE. X .LT. 1 */

L80:
    top = p[0];
    bot = q[0];
    for (i__ = 2; i__ <= 7; ++i__) {
	top = p[i__ - 1] + x * top;
	bot = q[i__ - 1] + x * bot;
/* L90: */
    }
    ret_val = top / bot;

/*     TERMINATION */

    if (*a < 1.) {
	goto L100;
    }
    ret_val *= t;
    return ret_val;
L100:
    ret_val /= t;
    return ret_val;
/* ----------------------------------------------------------------------- */
/*            EVALUATION OF GAMMA(A) FOR ABS(A) .GE. 15 */
/* ----------------------------------------------------------------------- */
L110:
    if (abs(*a) >= 1e3) {
	return ret_val;
    }
    if (*a > 0.) {
	goto L120;
    }
    x = -(*a);
    n = (integer) x;
    t = x - n;
    if (t > .9) {
	t = 1. - t;
    }
    s = sin(pi * t) / pi;
    if (n % 2 == 0) {
	s = -s;
    }
    if (s == 0.) {
	return ret_val;
    }

/*     COMPUTE THE MODIFIED ASYMPTOTIC SUM */

L120:
    t = 1. / (x * x);
    g = ((((r1 * t + r2) * t + r3) * t + r4) * t + r5) / x;

/*     ONE MAY REPLACE THE NEXT STATEMENT WITH  LNX = ALOG(X) */
/*     BUT LESS ACCURACY WILL NORMALLY BE OBTAINED. */

    lnx = log(x);

/*     FINAL ASSEMBLY */

    z__ = x;
    g = d__ + g + (z__ - .5) * (lnx - 1.);
    w = g;
    t = g - w;
    if (w > exparg_(&c__0) * .99999) {
	return ret_val;
    }
    ret_val = exp(w) * (t + 1.);
    if (*a < 0.) {
	ret_val = 1. / (ret_val * s) / x;
    }
    return ret_val;
} /* gamma_ */

