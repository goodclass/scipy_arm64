/* fpcuro.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int fpcuro_(doublereal *a, doublereal *b, doublereal *c__, 
	doublereal *d__, doublereal *x, integer *n)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double atan(doublereal), sqrt(doublereal), atan2(doublereal, doublereal), 
	    cos(doublereal), pow_dd(doublereal *, doublereal *), d_sign(
	    doublereal *, doublereal *);

    /* Local variables */
    static doublereal f;
    static integer i__;
    static doublereal q, r__, u, y, a1, b1, c1, d1, e3, p3, u1, u2, df, pi3, 
	    two, half, disc, ovfl, tent, four, step, three;

/*  subroutine fpcuro finds the real zeros of a cubic polynomial */
/*  p(x) = a*x**3+b*x**2+c*x+d. */

/*  calling sequence: */
/*     call fpcuro(a,b,c,d,x,n) */

/*  input parameters: */
/*    a,b,c,d: real values, containing the coefficients of p(x). */

/*  output parameters: */
/*    x      : real array,length 3, which contains the real zeros of p(x) */
/*    n      : integer, giving the number of real zeros of p(x). */
/*  .. */
/*  ..scalar arguments.. */
/*  ..array argument.. */
/*  ..local scalars.. */
/*  ..function references.. */
/*  set constants */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    two = 2.;
    three = 3.;
    four = 4.;
    ovfl = 1e4;
    half = .5;
    tent = .1;
    e3 = tent / .3;
    pi3 = atan(1.) / .75;
    a1 = abs(*a);
    b1 = abs(*b);
    c1 = abs(*c__);
    d1 = abs(*d__);
/*  test whether p(x) is a third degree polynomial. */
/* Computing MAX */
    d__1 = max(b1,c1);
    if (max(d__1,d1) < a1 * ovfl) {
	goto L300;
    }
/*  test whether p(x) is a second degree polynomial. */
    if (max(c1,d1) < b1 * ovfl) {
	goto L200;
    }
/*  test whether p(x) is a first degree polynomial. */
    if (d1 < c1 * ovfl) {
	goto L100;
    }
/*  p(x) is a constant function. */
    *n = 0;
    goto L800;
/*  p(x) is a first degree polynomial. */
L100:
    *n = 1;
    x[1] = -(*d__) / *c__;
    goto L500;
/*  p(x) is a second degree polynomial. */
L200:
    disc = *c__ * *c__ - four * *b * *d__;
    *n = 0;
    if (disc < 0.f) {
	goto L800;
    }
    *n = 2;
    u = sqrt(disc);
    b1 = *b + *b;
    x[1] = (-(*c__) + u) / b1;
    x[2] = (-(*c__) - u) / b1;
    goto L500;
/*  p(x) is a third degree polynomial. */
L300:
    b1 = *b / *a * e3;
    c1 = *c__ / *a;
    d1 = *d__ / *a;
    q = c1 * e3 - b1 * b1;
    r__ = b1 * b1 * b1 + (d1 - b1 * c1) * half;
    disc = q * q * q + r__ * r__;
    if (disc > 0.f) {
	goto L400;
    }
    u = sqrt((abs(q)));
    if (r__ < 0.f) {
	u = -u;
    }
    p3 = atan2(sqrt(-disc), (abs(r__))) * e3;
    u2 = u + u;
    *n = 3;
    x[1] = -u2 * cos(p3) - b1;
    x[2] = u2 * cos(pi3 - p3) - b1;
    x[3] = u2 * cos(pi3 + p3) - b1;
    goto L500;
L400:
    u = sqrt(disc);
    u1 = -r__ + u;
    u2 = -r__ - u;
    *n = 1;
    d__2 = abs(u1);
    d__1 = pow_dd(&d__2, &e3);
    d__4 = abs(u2);
    d__3 = pow_dd(&d__4, &e3);
    x[1] = d_sign(&d__1, &u1) + d_sign(&d__3, &u2) - b1;
/*  apply a newton iteration to improve the accuracy of the roots. */
L500:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y = x[i__];
	f = ((*a * y + *b) * y + *c__) * y + *d__;
	df = (three * *a * y + two * *b) * y + *c__;
	step = 0.f;
	if (abs(f) < abs(df) * tent) {
	    step = f / df;
	}
	x[i__] = y - step;
/* L700: */
    }
L800:
    return 0;
} /* fpcuro_ */

