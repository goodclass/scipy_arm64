/* bgrat.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int bgrat_(doublereal *a, doublereal *b, doublereal *x, 
	doublereal *y, doublereal *w, doublereal *eps, integer *ierr)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double log(doublereal), exp(doublereal);

    /* Local variables */
    static doublereal c__[30], d__[30];
    static integer i__;
    static doublereal j, l;
    static integer n;
    static doublereal p, q, r__, s, t, u, v, z__, n2, t2, dj, cn, nu, bm1;
    static integer nm1;
    static doublereal lnx, sum;
    extern doublereal gam1_(doublereal *);
    static doublereal bp2n, coef;
    extern /* Subroutine */ int grat1_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *);
    extern doublereal algdiv_(doublereal *, doublereal *), alnrel_(doublereal 
	    *);

/* ----------------------------------------------------------------------- */
/*     ASYMPTOTIC EXPANSION FOR IX(A,B) WHEN A IS LARGER THAN B. */
/*     THE RESULT OF THE EXPANSION IS ADDED TO W. IT IS ASSUMED */
/*     THAT A .GE. 15 AND B .LE. 1.  EPS IS THE TOLERANCE USED. */
/*     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS. */
/* ----------------------------------------------------------------------- */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

    bm1 = *b - .5 - .5;
    nu = *a + bm1 * .5;
    if (*y > .375) {
	goto L10;
    }
    d__1 = -(*y);
    lnx = alnrel_(&d__1);
    goto L20;
L10:
    lnx = log(*x);
L20:
    z__ = -nu * lnx;
    if (*b * z__ == 0.) {
	goto L70;
    }

/*                 COMPUTATION OF THE EXPANSION */
/*                 SET R = EXP(-Z)*Z**B/GAMMA(B) */

    r__ = *b * (gam1_(b) + 1.) * exp(*b * log(z__));
    r__ = r__ * exp(*a * lnx) * exp(bm1 * .5 * lnx);
    u = algdiv_(b, a) + *b * log(nu);
    u = r__ * exp(-u);
    if (u == 0.) {
	goto L70;
    }
    grat1_(b, &z__, &r__, &p, &q, eps);

/* Computing 2nd power */
    d__1 = 1. / nu;
    v = d__1 * d__1 * .25;
    t2 = lnx * .25 * lnx;
    l = *w / u;
    j = q / r__;
    sum = j;
    t = 1.;
    cn = 1.;
    n2 = 0.;
    for (n = 1; n <= 30; ++n) {
	bp2n = *b + n2;
	j = (bp2n * (bp2n + 1.) * j + (z__ + bp2n + 1.) * t) * v;
	n2 += 2.;
	t *= t2;
	cn /= n2 * (n2 + 1.);
	c__[n - 1] = cn;
	s = 0.;
	if (n == 1) {
	    goto L40;
	}
	nm1 = n - 1;
	coef = *b - n;
	i__1 = nm1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    s += coef * c__[i__ - 1] * d__[n - i__ - 1];
	    coef += *b;
/* L30: */
	}
L40:
	d__[n - 1] = bm1 * cn + s / n;
	dj = d__[n - 1] * j;
	sum += dj;
	if (sum <= 0.) {
	    goto L70;
	}
	if (abs(dj) <= *eps * (sum + l)) {
	    goto L60;
	}
/* L50: */
    }

/*                    ADD THE RESULTS TO W */

L60:
    *ierr = 0;
    *w += u * sum;
    return 0;

/*               THE EXPANSION CANNOT BE COMPUTED */

L70:
    *ierr = 1;
    return 0;
} /* bgrat_ */

