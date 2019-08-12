/* betaln.f -- translated by f2c (version 20190311).
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

doublereal betaln_(doublereal *a0, doublereal *b0)
{
    /* Initialized data */

    static doublereal e = .918938533204673;

    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static doublereal a, b, c__, h__;
    static integer i__, n;
    static doublereal u, v, w, z__;
    extern doublereal gamln_(doublereal *), bcorr_(doublereal *, doublereal *)
	    , algdiv_(doublereal *, doublereal *), alnrel_(doublereal *), 
	    gsumln_(doublereal *, doublereal *);

/* ----------------------------------------------------------------------- */
/*     EVALUATION OF THE LOGARITHM OF THE BETA FUNCTION */
/* ----------------------------------------------------------------------- */
/*     E = 0.5*LN(2*PI) */
/* -------------------------- */
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
/* -------------------------- */
    a = min(*a0,*b0);
    b = max(*a0,*b0);
    if (a >= 8.) {
	goto L100;
    }
    if (a >= 1.) {
	goto L20;
    }
/* ----------------------------------------------------------------------- */
/*                   PROCEDURE WHEN A .LT. 1 */
/* ----------------------------------------------------------------------- */
    if (b >= 8.) {
	goto L10;
    }
    d__1 = a + b;
    ret_val = gamln_(&a) + (gamln_(&b) - gamln_(&d__1));
    return ret_val;
L10:
    ret_val = gamln_(&a) + algdiv_(&a, &b);
    return ret_val;
/* ----------------------------------------------------------------------- */
/*                PROCEDURE WHEN 1 .LE. A .LT. 8 */
/* ----------------------------------------------------------------------- */
L20:
    if (a > 2.) {
	goto L40;
    }
    if (b > 2.) {
	goto L30;
    }
    ret_val = gamln_(&a) + gamln_(&b) - gsumln_(&a, &b);
    return ret_val;
L30:
    w = 0.;
    if (b < 8.) {
	goto L60;
    }
    ret_val = gamln_(&a) + algdiv_(&a, &b);
    return ret_val;

/*                REDUCTION OF A WHEN B .LE. 1000 */

L40:
    if (b > 1e3) {
	goto L80;
    }
    n = (integer) (a - 1.);
    w = 1.;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a += -1.;
	h__ = a / b;
	w *= h__ / (h__ + 1.);
/* L50: */
    }
    w = log(w);
    if (b < 8.) {
	goto L60;
    }
    ret_val = w + gamln_(&a) + algdiv_(&a, &b);
    return ret_val;

/*                 REDUCTION OF B WHEN B .LT. 8 */

L60:
    n = (integer) (b - 1.);
    z__ = 1.;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	b += -1.;
	z__ *= b / (a + b);
/* L70: */
    }
    ret_val = w + log(z__) + (gamln_(&a) + (gamln_(&b) - gsumln_(&a, &b)));
    return ret_val;

/*                REDUCTION OF A WHEN B .GT. 1000 */

L80:
    n = (integer) (a - 1.);
    w = 1.;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a += -1.;
	w *= a / (a / b + 1.);
/* L90: */
    }
    ret_val = log(w) - n * log(b) + (gamln_(&a) + algdiv_(&a, &b));
    return ret_val;
/* ----------------------------------------------------------------------- */
/*                   PROCEDURE WHEN A .GE. 8 */
/* ----------------------------------------------------------------------- */
L100:
    w = bcorr_(&a, &b);
    h__ = a / b;
    c__ = h__ / (h__ + 1.);
    u = -(a - .5) * log(c__);
    v = b * alnrel_(&h__);
    if (u <= v) {
	goto L110;
    }
    ret_val = log(b) * -.5 + e + w - v - u;
    return ret_val;
L110:
    ret_val = log(b) * -.5 + e + w - u - v;
    return ret_val;
} /* betaln_ */

