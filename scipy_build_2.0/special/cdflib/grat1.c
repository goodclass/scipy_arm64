/* grat1.f -- translated by f2c (version 20190311).
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

static integer c__0 = 0;

/* Subroutine */ int grat1_(doublereal *a, doublereal *x, doublereal *r__, 
	doublereal *p, doublereal *q, doublereal *eps)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double log(doublereal), exp(doublereal), sqrt(doublereal);

    /* Local variables */
    static doublereal c__, g, h__, j, l, t, w, z__, an, am0, an0, a2n, b2n, 
	    cma;
    extern doublereal erf_(doublereal *);
    static doublereal tol, sum;
    extern doublereal gam1_(doublereal *);
    static doublereal a2nm1, b2nm1;
    extern doublereal rexp_(doublereal *), erfc1_(integer *, doublereal *);

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */
/* ----------------------------------------------------------------------- */
/*        EVALUATION OF THE INCOMPLETE GAMMA RATIO FUNCTIONS */
/*                      P(A,X) AND Q(A,X) */

/*     IT IS ASSUMED THAT A .LE. 1.  EPS IS THE TOLERANCE TO BE USED. */
/*     THE INPUT ARGUMENT R HAS THE VALUE E**(-X)*X**A/GAMMA(A). */
/* ----------------------------------------------------------------------- */
    if (*a * *x == 0.) {
	goto L120;
    }
    if (*a == .5) {
	goto L100;
    }
    if (*x < 1.1) {
	goto L10;
    }
    goto L60;

/*             TAYLOR SERIES FOR P(A,X)/X**A */

L10:
    an = 3.;
    c__ = *x;
    sum = *x / (*a + 3.);
    tol = *eps * .1 / (*a + 1.);
L20:
    an += 1.;
    c__ = -c__ * (*x / an);
    t = c__ / (*a + an);
    sum += t;
    if (abs(t) > tol) {
	goto L20;
    }
    j = *a * *x * ((sum / 6. - .5 / (*a + 2.)) * *x + 1. / (*a + 1.));

    z__ = *a * log(*x);
    h__ = gam1_(a);
    g = h__ + 1.;
    if (*x < .25) {
	goto L30;
    }
    if (*a < *x / 2.59) {
	goto L50;
    }
    goto L40;
L30:
    if (z__ > -.13394) {
	goto L50;
    }

L40:
    w = exp(z__);
    *p = w * g * (.5 - j + .5);
    *q = .5 - *p + .5;
    return 0;

L50:
    l = rexp_(&z__);
    w = l + .5 + .5;
    *q = (w * j - l) * g - h__;
    if (*q < 0.) {
	goto L90;
    }
    *p = .5 - *q + .5;
    return 0;

/*              CONTINUED FRACTION EXPANSION */

L60:
    a2nm1 = 1.;
    a2n = 1.;
    b2nm1 = *x;
    b2n = *x + (1. - *a);
    c__ = 1.;
L70:
    a2nm1 = *x * a2n + c__ * a2nm1;
    b2nm1 = *x * b2n + c__ * b2nm1;
    am0 = a2nm1 / b2nm1;
    c__ += 1.;
    cma = c__ - *a;
    a2n = a2nm1 + cma * a2n;
    b2n = b2nm1 + cma * b2n;
    an0 = a2n / b2n;
    if ((d__1 = an0 - am0, abs(d__1)) >= *eps * an0) {
	goto L70;
    }
    *q = *r__ * an0;
    *p = .5 - *q + .5;
    return 0;

/*                SPECIAL CASES */

L80:
    *p = 0.;
    *q = 1.;
    return 0;

L90:
    *p = 1.;
    *q = 0.;
    return 0;

L100:
    if (*x >= .25) {
	goto L110;
    }
    d__1 = sqrt(*x);
    *p = erf_(&d__1);
    *q = .5 - *p + .5;
    return 0;
L110:
    d__1 = sqrt(*x);
    *q = erfc1_(&c__0, &d__1);
    *p = .5 - *q + .5;
    return 0;

L120:
    if (*x <= *a) {
	goto L80;
    }
    goto L90;
} /* grat1_ */

