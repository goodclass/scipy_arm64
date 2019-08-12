/* bfrac.f -- translated by f2c (version 20190311).
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

doublereal bfrac_(doublereal *a, doublereal *b, doublereal *x, doublereal *y, 
	doublereal *lambda, doublereal *eps)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Local variables */
    static doublereal c__, e, n, p, r__, s, t, w, c0, c1, r0, an, bn, yp1, 
	    anp1, bnp1, beta, alpha;
    extern doublereal brcomp_(doublereal *, doublereal *, doublereal *, 
	    doublereal *);

/* ----------------------------------------------------------------------- */
/*     CONTINUED FRACTION EXPANSION FOR IX(A,B) WHEN A,B .GT. 1. */
/*     IT IS ASSUMED THAT  LAMBDA = (A + B)*Y - B. */
/* ----------------------------------------------------------------------- */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */
/* -------------------- */
    ret_val = brcomp_(a, b, x, y);
    if (ret_val == 0.) {
	return ret_val;
    }

    c__ = *lambda + 1.;
    c0 = *b / *a;
    c1 = 1. / *a + 1.;
    yp1 = *y + 1.;

    n = 0.;
    p = 1.;
    s = *a + 1.;
    an = 0.;
    bn = 1.;
    anp1 = 1.;
    bnp1 = c__ / c1;
    r__ = c1 / c__;

/*        CONTINUED FRACTION CALCULATION */

L10:
    n += 1.;
    t = n / *a;
    w = n * (*b - n) * *x;
    e = *a / s;
    alpha = p * (p + c0) * e * e * (w * *x);
    e = (t + 1.) / (c1 + t + t);
    beta = n + w / s + e * (c__ + n * yp1);
    p = t + 1.;
    s += 2.;

/*        UPDATE AN, BN, ANP1, AND BNP1 */

    t = alpha * an + beta * anp1;
    an = anp1;
    anp1 = t;
    t = alpha * bn + beta * bnp1;
    bn = bnp1;
    bnp1 = t;

    r0 = r__;
    r__ = anp1 / bnp1;
    if (! ((d__1 = r__ - r0, abs(d__1)) > *eps * r__)) {
	goto L20;
    }

/*        RESCALE AN, BN, ANP1, AND BNP1 */

    an /= bnp1;
    bn /= bnp1;
    anp1 = r__;
    bnp1 = 1.;
    goto L10;

/*                 TERMINATION */

L20:
    ret_val *= r__;
    return ret_val;
} /* bfrac_ */

