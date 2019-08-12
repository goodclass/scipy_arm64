/* brcmp1.f -- translated by f2c (version 20190311).
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

doublereal brcmp1_(integer *mu, doublereal *a, doublereal *b, doublereal *x, 
	doublereal *y)
{
    /* Initialized data */

    static doublereal const__ = .398942280401433;

    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double log(doublereal), sqrt(doublereal), exp(doublereal);

    /* Local variables */
    static doublereal c__, e, h__;
    static integer i__, n;
    static doublereal t, u, v, z__, a0, b0, x0, y0, apb, lnx, lny;
    extern doublereal gam1_(doublereal *), esum_(integer *, doublereal *), 
	    rlog1_(doublereal *), bcorr_(doublereal *, doublereal *), gamln1_(
	    doublereal *);
    static doublereal lambda;
    extern doublereal betaln_(doublereal *, doublereal *), algdiv_(doublereal 
	    *, doublereal *), alnrel_(doublereal *);

/* ----------------------------------------------------------------------- */
/*          EVALUATION OF  EXP(MU) * (X**A*Y**B/BETA(A,B)) */
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
/* ----------------- */
/*     CONST = 1/SQRT(2*PI) */
/* ----------------- */
/*     .. */
/*     .. Executable Statements .. */

    a0 = min(*a,*b);
    if (a0 >= 8.) {
	goto L130;
    }

    if (*x > .375) {
	goto L10;
    }
    lnx = log(*x);
    d__1 = -(*x);
    lny = alnrel_(&d__1);
    goto L30;
L10:
    if (*y > .375) {
	goto L20;
    }
    d__1 = -(*y);
    lnx = alnrel_(&d__1);
    lny = log(*y);
    goto L30;
L20:
    lnx = log(*x);
    lny = log(*y);

L30:
    z__ = *a * lnx + *b * lny;
    if (a0 < 1.) {
	goto L40;
    }
    z__ -= betaln_(a, b);
    ret_val = esum_(mu, &z__);
    return ret_val;
/* ----------------------------------------------------------------------- */
/*              PROCEDURE FOR A .LT. 1 OR B .LT. 1 */
/* ----------------------------------------------------------------------- */
L40:
    b0 = max(*a,*b);
    if (b0 >= 8.) {
	goto L120;
    }
    if (b0 > 1.) {
	goto L70;
    }

/*                   ALGORITHM FOR B0 .LE. 1 */

    ret_val = esum_(mu, &z__);
    if (ret_val == 0.) {
	return ret_val;
    }

    apb = *a + *b;
    if (apb > 1.) {
	goto L50;
    }
    z__ = gam1_(&apb) + 1.;
    goto L60;
L50:
    u = *a + *b - 1.;
    z__ = (gam1_(&u) + 1.) / apb;

L60:
    c__ = (gam1_(a) + 1.) * (gam1_(b) + 1.) / z__;
    ret_val = ret_val * (a0 * c__) / (a0 / b0 + 1.);
    return ret_val;

/*                ALGORITHM FOR 1 .LT. B0 .LT. 8 */

L70:
    u = gamln1_(&a0);
    n = (integer) (b0 - 1.);
    if (n < 1) {
	goto L90;
    }
    c__ = 1.;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	b0 += -1.;
	c__ *= b0 / (a0 + b0);
/* L80: */
    }
    u = log(c__) + u;

L90:
    z__ -= u;
    b0 += -1.;
    apb = a0 + b0;
    if (apb > 1.) {
	goto L100;
    }
    t = gam1_(&apb) + 1.;
    goto L110;
L100:
    u = a0 + b0 - 1.;
    t = (gam1_(&u) + 1.) / apb;
L110:
    ret_val = a0 * esum_(mu, &z__) * (gam1_(&b0) + 1.) / t;
    return ret_val;

/*                   ALGORITHM FOR B0 .GE. 8 */

L120:
    u = gamln1_(&a0) + algdiv_(&a0, &b0);
    d__1 = z__ - u;
    ret_val = a0 * esum_(mu, &d__1);
    return ret_val;
/* ----------------------------------------------------------------------- */
/*              PROCEDURE FOR A .GE. 8 AND B .GE. 8 */
/* ----------------------------------------------------------------------- */
L130:
    if (*a > *b) {
	goto L140;
    }
    h__ = *a / *b;
    x0 = h__ / (h__ + 1.);
    y0 = 1. / (h__ + 1.);
    lambda = *a - (*a + *b) * *x;
    goto L150;
L140:
    h__ = *b / *a;
    x0 = 1. / (h__ + 1.);
    y0 = h__ / (h__ + 1.);
    lambda = (*a + *b) * *y - *b;

L150:
    e = -lambda / *a;
    if (abs(e) > .6) {
	goto L160;
    }
    u = rlog1_(&e);
    goto L170;
L160:
    u = e - log(*x / x0);

L170:
    e = lambda / *b;
    if (abs(e) > .6) {
	goto L180;
    }
    v = rlog1_(&e);
    goto L190;
L180:
    v = e - log(*y / y0);

L190:
    d__1 = -(*a * u + *b * v);
    z__ = esum_(mu, &d__1);
    ret_val = const__ * sqrt(*b * x0) * z__ * exp(-bcorr_(a, b));
    return ret_val;
} /* brcmp1_ */

