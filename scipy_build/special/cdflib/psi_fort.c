/* psi_fort.f -- translated by f2c (version 20190311).
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
static integer c__1 = 1;

doublereal psi_(doublereal *xx)
{
    /* Initialized data */

    static doublereal piov4 = .785398163397448;
    static doublereal dx0 = 1.461632144968362341262659542325721325;
    static doublereal p1[7] = { .0089538502298197,4.77762828042627,
	    142.441585084029,1186.45200713425,3633.51846806499,
	    4138.10161269013,1305.60269827897 };
    static doublereal q1[6] = { 44.8452573429826,520.752771467162,
	    2210.0079924783,3641.27349079381,1908.310765963,
	    6.91091682714533e-6 };
    static doublereal p2[4] = { -2.12940445131011,-7.01677227766759,
	    -4.48616543918019,-.648157123766197 };
    static doublereal q2[4] = { 32.2703493791143,89.2920700481861,
	    54.6117738103215,7.77788548522962 };

    /* System generated locals */
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal), log(doublereal);

    /* Local variables */
    static integer i__, m, n;
    static doublereal w, x, z__;
    static integer nq;
    static doublereal den, aug, sgn, xmx0, xmax1, upper;
    extern integer ipmpar_(integer *);
    static doublereal xsmall;
    extern doublereal spmpar_(integer *);

/* --------------------------------------------------------------------- */

/*                 EVALUATION OF THE DIGAMMA FUNCTION */

/*                           ----------- */

/*     PSI(XX) IS ASSIGNED THE VALUE 0 WHEN THE DIGAMMA FUNCTION CANNOT */
/*     BE COMPUTED. */

/*     THE MAIN COMPUTATION INVOLVES EVALUATION OF RATIONAL CHEBYSHEV */
/*     APPROXIMATIONS PUBLISHED IN MATH. COMP. 27, 123-127(1973) BY */
/*     CODY, STRECOK AND THACHER. */

/* --------------------------------------------------------------------- */
/*     PSI WAS WRITTEN AT ARGONNE NATIONAL LABORATORY FOR THE FUNPACK */
/*     PACKAGE OF SPECIAL FUNCTION SUBROUTINES. PSI WAS MODIFIED BY */
/*     A.H. MORRIS (NSWC). */
/* --------------------------------------------------------------------- */
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
/* --------------------------------------------------------------------- */

/*     PIOV4 = PI/4 */
/*     DX0 = ZERO OF PSI TO EXTENDED PRECISION */

/* --------------------------------------------------------------------- */
/* --------------------------------------------------------------------- */

/*     COEFFICIENTS FOR RATIONAL APPROXIMATION OF */
/*     PSI(X) / (X - X0),  0.5 .LE. X .LE. 3.0 */

/* --------------------------------------------------------------------- */
/* --------------------------------------------------------------------- */

/*     COEFFICIENTS FOR RATIONAL APPROXIMATION OF */
/*     PSI(X) - LN(X) + 1 / (2*X),  X .GT. 3.0 */

/* --------------------------------------------------------------------- */
/*     .. */
/*     .. Executable Statements .. */
/* --------------------------------------------------------------------- */

/*     MACHINE DEPENDENT CONSTANTS ... */

/*        XMAX1  = THE SMALLEST POSITIVE FLOATING POINT CONSTANT */
/*                 WITH ENTIRELY INTEGER REPRESENTATION.  ALSO USED */
/*                 AS NEGATIVE OF LOWER BOUND ON ACCEPTABLE NEGATIVE */
/*                 ARGUMENTS AND AS THE POSITIVE ARGUMENT BEYOND WHICH */
/*                 PSI MAY BE REPRESENTED AS ALOG(X). */

/*        XSMALL = ABSOLUTE ARGUMENT BELOW WHICH PI*COTAN(PI*X) */
/*                 MAY BE REPRESENTED BY 1/X. */

/* --------------------------------------------------------------------- */
    xmax1 = (doublereal) ipmpar_(&c__3);
/* Computing MIN */
    d__1 = xmax1, d__2 = 1. / spmpar_(&c__1);
    xmax1 = min(d__1,d__2);
    xsmall = 1e-9;
/* --------------------------------------------------------------------- */
    x = *xx;
    aug = 0.;
    if (x >= .5) {
	goto L50;
    }
/* --------------------------------------------------------------------- */
/*     X .LT. 0.5,  USE REFLECTION FORMULA */
/*     PSI(1-X) = PSI(X) + PI * COTAN(PI*X) */
/* --------------------------------------------------------------------- */
    if (abs(x) > xsmall) {
	goto L10;
    }
    if (x == 0.) {
	goto L100;
    }
/* --------------------------------------------------------------------- */
/*     0 .LT. ABS(X) .LE. XSMALL.  USE 1/X AS A SUBSTITUTE */
/*     FOR  PI*COTAN(PI*X) */
/* --------------------------------------------------------------------- */
    aug = -1. / x;
    goto L40;
/* --------------------------------------------------------------------- */
/*     REDUCTION OF ARGUMENT FOR COTAN */
/* --------------------------------------------------------------------- */
L10:
    w = -x;
    sgn = piov4;
    if (w > 0.) {
	goto L20;
    }
    w = -w;
    sgn = -sgn;
/* --------------------------------------------------------------------- */
/*     MAKE AN ERROR EXIT IF X .LE. -XMAX1 */
/* --------------------------------------------------------------------- */
L20:
    if (w >= xmax1) {
	goto L100;
    }
    nq = (integer) w;
    w -= (doublereal) nq;
    nq = (integer) (w * 4.);
    w = (w - (doublereal) nq * .25) * 4.;
/* --------------------------------------------------------------------- */
/*     W IS NOW RELATED TO THE FRACTIONAL PART OF  4.0 * X. */
/*     ADJUST ARGUMENT TO CORRESPOND TO VALUES IN FIRST */
/*     QUADRANT AND DETERMINE SIGN */
/* --------------------------------------------------------------------- */
    n = nq / 2;
    if (n + n != nq) {
	w = 1. - w;
    }
    z__ = piov4 * w;
    m = n / 2;
    if (m + m != n) {
	sgn = -sgn;
    }
/* --------------------------------------------------------------------- */
/*     DETERMINE FINAL VALUE FOR  -PI*COTAN(PI*X) */
/* --------------------------------------------------------------------- */
    n = (nq + 1) / 2;
    m = n / 2;
    m += m;
    if (m != n) {
	goto L30;
    }
/* --------------------------------------------------------------------- */
/*     CHECK FOR SINGULARITY */
/* --------------------------------------------------------------------- */
    if (z__ == 0.) {
	goto L100;
    }
/* --------------------------------------------------------------------- */
/*     USE COS/SIN AS A SUBSTITUTE FOR COTAN, AND */
/*     SIN/COS AS A SUBSTITUTE FOR TAN */
/* --------------------------------------------------------------------- */
    aug = sgn * (cos(z__) / sin(z__) * 4.);
    goto L40;
L30:
    aug = sgn * (sin(z__) / cos(z__) * 4.);
L40:
    x = 1. - x;
L50:
    if (x > 3.) {
	goto L70;
    }
/* --------------------------------------------------------------------- */
/*     0.5 .LE. X .LE. 3.0 */
/* --------------------------------------------------------------------- */
    den = x;
    upper = p1[0] * x;

    for (i__ = 1; i__ <= 5; ++i__) {
	den = (den + q1[i__ - 1]) * x;
	upper = (upper + p1[i__]) * x;
/* L60: */
    }

    den = (upper + p1[6]) / (den + q1[5]);
    xmx0 = x - dx0;
    ret_val = den * xmx0 + aug;
    return ret_val;
/* --------------------------------------------------------------------- */
/*     IF X .GE. XMAX1, PSI = LN(X) */
/* --------------------------------------------------------------------- */
L70:
    if (x >= xmax1) {
	goto L90;
    }
/* --------------------------------------------------------------------- */
/*     3.0 .LT. X .LT. XMAX1 */
/* --------------------------------------------------------------------- */
    w = 1. / (x * x);
    den = w;
    upper = p2[0] * w;

    for (i__ = 1; i__ <= 3; ++i__) {
	den = (den + q2[i__ - 1]) * w;
	upper = (upper + p2[i__]) * w;
/* L80: */
    }

    aug = upper / (den + q2[3]) - .5 / x + aug;
L90:
    ret_val = aug + log(x);
    return ret_val;
/* --------------------------------------------------------------------- */
/*     ERROR RETURN */
/* --------------------------------------------------------------------- */
L100:
    ret_val = 0.;
    return ret_val;
} /* psi_ */

