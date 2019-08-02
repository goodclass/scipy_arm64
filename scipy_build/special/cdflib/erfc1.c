/* erfc1.f -- translated by f2c (version 20190311).
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

doublereal erfc1_(integer *ind, doublereal *x)
{
    /* Initialized data */

    static doublereal c__ = .564189583547756;
    static doublereal a[5] = { 7.7105849500132e-5,-.00133733772997339,
	    .0323076579225834,.0479137145607681,.128379167095513 };
    static doublereal b[3] = { .00301048631703895,.0538971687740286,
	    .375795757275549 };
    static doublereal p[8] = { -1.36864857382717e-7,.564195517478974,
	    7.21175825088309,43.1622272220567,152.98928504694,
	    339.320816734344,451.918953711873,300.459261020162 };
    static doublereal q[8] = { 1.,12.7827273196294,77.0001529352295,
	    277.585444743988,638.980264465631,931.35409485061,
	    790.950925327898,300.459260956983 };
    static doublereal r__[5] = { 2.10144126479064,26.2370141675169,
	    21.3688200555087,4.6580782871847,.282094791773523 };
    static doublereal s[4] = { 94.153775055546,187.11481179959,
	    99.0191814623914,18.0124575948747 };

    /* System generated locals */
    doublereal ret_val, d__1;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static doublereal e, t, w, ax, bot, top;
    extern doublereal exparg_(integer *);

/* ----------------------------------------------------------------------- */
/*         EVALUATION OF THE COMPLEMENTARY ERROR FUNCTION */

/*          ERFC1(IND,X) = ERFC(X)            IF IND = 0 */
/*          ERFC1(IND,X) = EXP(X*X)*ERFC(X)   OTHERWISE */
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
/* ------------------------- */
/* ------------------------- */
/* ------------------------- */
/* ------------------------- */
/*     .. */
/*     .. Executable Statements .. */
/* ------------------------- */

/*                     ABS(X) .LE. 0.5 */

    ax = abs(*x);
    if (ax > .5) {
	goto L10;
    }
    t = *x * *x;
    top = (((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4] + 1.;
    bot = ((b[0] * t + b[1]) * t + b[2]) * t + 1.;
    ret_val = .5 - *x * (top / bot) + .5;
    if (*ind != 0) {
	ret_val = exp(t) * ret_val;
    }
    return ret_val;

/*                  0.5 .LT. ABS(X) .LE. 4 */

L10:
    if (ax > 4.) {
	goto L20;
    }
    top = ((((((p[0] * ax + p[1]) * ax + p[2]) * ax + p[3]) * ax + p[4]) * ax 
	    + p[5]) * ax + p[6]) * ax + p[7];
    bot = ((((((q[0] * ax + q[1]) * ax + q[2]) * ax + q[3]) * ax + q[4]) * ax 
	    + q[5]) * ax + q[6]) * ax + q[7];
    ret_val = top / bot;
    goto L40;

/*                      ABS(X) .GT. 4 */

L20:
    if (*x <= -5.6) {
	goto L60;
    }
    if (*ind != 0) {
	goto L30;
    }
    if (*x > 100.) {
	goto L70;
    }
    if (*x * *x > -exparg_(&c__1)) {
	goto L70;
    }

L30:
/* Computing 2nd power */
    d__1 = 1. / *x;
    t = d__1 * d__1;
    top = (((r__[0] * t + r__[1]) * t + r__[2]) * t + r__[3]) * t + r__[4];
    bot = (((s[0] * t + s[1]) * t + s[2]) * t + s[3]) * t + 1.;
    ret_val = (c__ - t * top / bot) / ax;

/*                      FINAL ASSEMBLY */

L40:
    if (*ind == 0) {
	goto L50;
    }
    if (*x < 0.) {
	ret_val = exp(*x * *x) * 2. - ret_val;
    }
    return ret_val;
L50:
    w = *x * *x;
    t = w;
    e = w - t;
    ret_val = (.5 - e + .5) * exp(-t) * ret_val;
    if (*x < 0.) {
	ret_val = 2. - ret_val;
    }
    return ret_val;

/*             LIMIT VALUE FOR LARGE NEGATIVE X */

L60:
    ret_val = 2.;
    if (*ind != 0) {
	ret_val = exp(*x * *x) * 2.;
    }
    return ret_val;

/*             LIMIT VALUE FOR LARGE POSITIVE X */
/*                       WHEN IND = 0 */

L70:
    ret_val = 0.;
    return ret_val;
} /* erfc1_ */

