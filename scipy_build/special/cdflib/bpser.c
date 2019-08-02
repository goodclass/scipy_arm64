/* bpser.f -- translated by f2c (version 20190311).
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

doublereal bpser_(doublereal *a, doublereal *b, doublereal *x, doublereal *
	eps)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Builtin functions */
    double log(doublereal), exp(doublereal), pow_dd(doublereal *, doublereal *
	    );

    /* Local variables */
    static doublereal c__;
    static integer i__, m;
    static doublereal n, t, u, w, z__, a0, b0, apb, tol, sum;
    extern doublereal gam1_(doublereal *), gamln1_(doublereal *), betaln_(
	    doublereal *, doublereal *), algdiv_(doublereal *, doublereal *);

/* ----------------------------------------------------------------------- */
/*     POWER SERIES EXPANSION FOR EVALUATING IX(A,B) WHEN B .LE. 1 */
/*     OR B*X .LE. 0.7.  EPS IS THE TOLERANCE USED. */
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

    ret_val = 0.;
    if (*x == 0.) {
	return ret_val;
    }
/* ----------------------------------------------------------------------- */
/*            COMPUTE THE FACTOR X**A/(A*BETA(A,B)) */
/* ----------------------------------------------------------------------- */
    a0 = min(*a,*b);
    if (a0 < 1.) {
	goto L10;
    }
    z__ = *a * log(*x) - betaln_(a, b);
    ret_val = exp(z__) / *a;
    goto L100;
L10:
    b0 = max(*a,*b);
    if (b0 >= 8.) {
	goto L90;
    }
    if (b0 > 1.) {
	goto L40;
    }

/*            PROCEDURE FOR A0 .LT. 1 AND B0 .LE. 1 */

    ret_val = pow_dd(x, a);
    if (ret_val == 0.) {
	return ret_val;
    }

    apb = *a + *b;
    if (apb > 1.) {
	goto L20;
    }
    z__ = gam1_(&apb) + 1.;
    goto L30;
L20:
    u = *a + *b - 1.;
    z__ = (gam1_(&u) + 1.) / apb;

L30:
    c__ = (gam1_(a) + 1.) * (gam1_(b) + 1.) / z__;
    ret_val = ret_val * c__ * (*b / apb);
    goto L100;

/*         PROCEDURE FOR A0 .LT. 1 AND 1 .LT. B0 .LT. 8 */

L40:
    u = gamln1_(&a0);
    m = (integer) (b0 - 1.);
    if (m < 1) {
	goto L60;
    }
    c__ = 1.;
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	b0 += -1.;
	c__ *= b0 / (a0 + b0);
/* L50: */
    }
    u = log(c__) + u;

L60:
    z__ = *a * log(*x) - u;
    b0 += -1.;
    apb = a0 + b0;
    if (apb > 1.) {
	goto L70;
    }
    t = gam1_(&apb) + 1.;
    goto L80;
L70:
    u = a0 + b0 - 1.;
    t = (gam1_(&u) + 1.) / apb;
L80:
    ret_val = exp(z__) * (a0 / *a) * (gam1_(&b0) + 1.) / t;
    goto L100;

/*            PROCEDURE FOR A0 .LT. 1 AND B0 .GE. 8 */

L90:
    u = gamln1_(&a0) + algdiv_(&a0, &b0);
    z__ = *a * log(*x) - u;
    ret_val = a0 / *a * exp(z__);
L100:
    if (ret_val == 0. || *a <= *eps * .1) {
	return ret_val;
    }
/* ----------------------------------------------------------------------- */
/*                     COMPUTE THE SERIES */
/* ----------------------------------------------------------------------- */
    sum = 0.;
    n = 0.;
    c__ = 1.;
    tol = *eps / *a;
L110:
    n += 1.;
    c__ = c__ * (.5 - *b / n + .5) * *x;
    w = c__ / (*a + n);
    sum += w;
    if (abs(w) > tol) {
	goto L110;
    }
    ret_val *= *a * sum + 1.;
    return ret_val;
} /* bpser_ */

