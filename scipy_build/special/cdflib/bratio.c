/* bratio.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int bratio_(doublereal *a, doublereal *b, doublereal *x, 
	doublereal *y, doublereal *w, doublereal *w1, integer *ierr)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static integer n;
    static doublereal t, z__, a0, b0, x0, y0;
    static integer ind;
    extern doublereal bup_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *);
    static doublereal eps;
    static integer ierr1;
    extern doublereal bfrac_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    extern /* Subroutine */ int bgrat_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, doublereal *, integer *);
    extern doublereal apser_(doublereal *, doublereal *, doublereal *, 
	    doublereal *), basym_(doublereal *, doublereal *, doublereal *, 
	    doublereal *), bpser_(doublereal *, doublereal *, doublereal *, 
	    doublereal *), fpser_(doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal lambda;
    extern doublereal spmpar_(integer *);

/* ----------------------------------------------------------------------- */

/*            EVALUATION OF THE INCOMPLETE BETA FUNCTION IX(A,B) */

/*                     -------------------- */

/*     IT IS ASSUMED THAT A AND B ARE NONNEGATIVE, AND THAT X .LE. 1 */
/*     AND Y = 1 - X.  BRATIO ASSIGNS W AND W1 THE VALUES */

/*                      W  = IX(A,B) */
/*                      W1 = 1 - IX(A,B) */

/*     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS. */
/*     IF NO INPUT ERRORS ARE DETECTED THEN IERR IS SET TO 0 AND */
/*     W AND W1 ARE COMPUTED. OTHERWISE, IF AN ERROR IS DETECTED, */
/*     THEN W AND W1 ARE ASSIGNED THE VALUE 0 AND IERR IS SET TO */
/*     ONE OF THE FOLLOWING VALUES ... */

/*        IERR = 1  IF A OR B IS NEGATIVE */
/*        IERR = 2  IF A = B = 0 */
/*        IERR = 3  IF X .LT. 0 OR X .GT. 1 */
/*        IERR = 4  IF Y .LT. 0 OR Y .GT. 1 */
/*        IERR = 5  IF X + Y .NE. 1 */
/*        IERR = 6  IF X = A = 0 */
/*        IERR = 7  IF Y = B = 0 */

/* -------------------- */
/*     WRITTEN BY ALFRED H. MORRIS, JR. */
/*        NAVAL SURFACE WARFARE CENTER */
/*        DAHLGREN, VIRGINIA */
/*     REVISED ... NOV 1991 */
/* ----------------------------------------------------------------------- */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */
/* ----------------------------------------------------------------------- */

/*     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE SMALLEST */
/*            FLOATING POINT NUMBER FOR WHICH 1.0 + EPS .GT. 1.0 */

    eps = spmpar_(&c__1);

/* ----------------------------------------------------------------------- */
    *w = 0.;
    *w1 = 0.;
    if (*a < 0. || *b < 0.) {
	goto L270;
    }
    if (*a == 0. && *b == 0.) {
	goto L280;
    }
    if (*x < 0. || *x > 1.) {
	goto L290;
    }
    if (*y < 0. || *y > 1.) {
	goto L300;
    }
    z__ = *x + *y - .5 - .5;
    if (abs(z__) > eps * 3.) {
	goto L310;
    }

    *ierr = 0;
    if (*x == 0.) {
	goto L210;
    }
    if (*y == 0.) {
	goto L230;
    }
    if (*a == 0.) {
	goto L240;
    }
    if (*b == 0.) {
	goto L220;
    }

    eps = max(eps,1e-15);
    if (max(*a,*b) < eps * .001) {
	goto L260;
    }

    ind = 0;
    a0 = *a;
    b0 = *b;
    x0 = *x;
    y0 = *y;
    if (min(a0,b0) > 1.) {
	goto L40;
    }

/*             PROCEDURE FOR A0 .LE. 1 OR B0 .LE. 1 */

    if (*x <= .5) {
	goto L10;
    }
    ind = 1;
    a0 = *b;
    b0 = *a;
    x0 = *y;
    y0 = *x;

L10:
/* Computing MIN */
    d__1 = eps, d__2 = eps * a0;
    if (b0 < min(d__1,d__2)) {
	goto L90;
    }
/* Computing MIN */
    d__1 = eps, d__2 = eps * b0;
    if (a0 < min(d__1,d__2) && b0 * x0 <= 1.) {
	goto L100;
    }
    if (max(a0,b0) > 1.) {
	goto L20;
    }
    if (a0 >= min(.2,b0)) {
	goto L110;
    }
    if (pow_dd(&x0, &a0) <= .9) {
	goto L110;
    }
    if (x0 >= .3) {
	goto L120;
    }
    n = 20;
    goto L140;

L20:
    if (b0 <= 1.) {
	goto L110;
    }
    if (x0 >= .3) {
	goto L120;
    }
    if (x0 >= .1) {
	goto L30;
    }
    d__1 = x0 * b0;
    if (pow_dd(&d__1, &a0) <= .7) {
	goto L110;
    }
L30:
    if (b0 > 15.) {
	goto L150;
    }
    n = 20;
    goto L140;

/*             PROCEDURE FOR A0 .GT. 1 AND B0 .GT. 1 */

L40:
    if (*a > *b) {
	goto L50;
    }
    lambda = *a - (*a + *b) * *x;
    goto L60;
L50:
    lambda = (*a + *b) * *y - *b;
L60:
    if (lambda >= 0.) {
	goto L70;
    }
    ind = 1;
    a0 = *b;
    b0 = *a;
    x0 = *y;
    y0 = *x;
    lambda = abs(lambda);

L70:
    if (b0 < 40. && b0 * x0 <= .7) {
	goto L110;
    }
    if (b0 < 40.) {
	goto L160;
    }
    if (a0 > b0) {
	goto L80;
    }
    if (a0 <= 100.) {
	goto L130;
    }
    if (lambda > a0 * .03) {
	goto L130;
    }
    goto L200;
L80:
    if (b0 <= 100.) {
	goto L130;
    }
    if (lambda > b0 * .03) {
	goto L130;
    }
    goto L200;

/*            EVALUATION OF THE APPROPRIATE ALGORITHM */

L90:
    *w = fpser_(&a0, &b0, &x0, &eps);
    *w1 = .5 - *w + .5;
    goto L250;

L100:
    *w1 = apser_(&a0, &b0, &x0, &eps);
    *w = .5 - *w1 + .5;
    goto L250;

L110:
    *w = bpser_(&a0, &b0, &x0, &eps);
    *w1 = .5 - *w + .5;
    goto L250;

L120:
    *w1 = bpser_(&b0, &a0, &y0, &eps);
    *w = .5 - *w1 + .5;
    goto L250;

L130:
    d__1 = eps * 15.;
    *w = bfrac_(&a0, &b0, &x0, &y0, &lambda, &d__1);
    *w1 = .5 - *w + .5;
    goto L250;

L140:
    *w1 = bup_(&b0, &a0, &y0, &x0, &n, &eps);
    b0 += n;
L150:
    d__1 = eps * 15.;
    bgrat_(&b0, &a0, &y0, &x0, w1, &d__1, &ierr1);
    *w = .5 - *w1 + .5;
    goto L250;

L160:
    n = (integer) b0;
    b0 -= n;
    if (b0 != 0.) {
	goto L170;
    }
    --n;
    b0 = 1.;
L170:
    *w = bup_(&b0, &a0, &y0, &x0, &n, &eps);
    if (x0 > .7) {
	goto L180;
    }
    *w += bpser_(&a0, &b0, &x0, &eps);
    *w1 = .5 - *w + .5;
    goto L250;

L180:
    if (a0 > 15.) {
	goto L190;
    }
    n = 20;
    *w += bup_(&a0, &b0, &x0, &y0, &n, &eps);
    a0 += n;
L190:
    d__1 = eps * 15.;
    bgrat_(&a0, &b0, &x0, &y0, w, &d__1, &ierr1);
    *w1 = .5 - *w + .5;
    goto L250;

L200:
    d__1 = eps * 100.;
    *w = basym_(&a0, &b0, &lambda, &d__1);
    *w1 = .5 - *w + .5;
    goto L250;

/*               TERMINATION OF THE PROCEDURE */

L210:
    if (*a == 0.) {
	goto L320;
    }
L220:
    *w = 0.;
    *w1 = 1.;
    return 0;

L230:
    if (*b == 0.) {
	goto L330;
    }
L240:
    *w = 1.;
    *w1 = 0.;
    return 0;

L250:
    if (ind == 0) {
	return 0;
    }
    t = *w;
    *w = *w1;
    *w1 = t;
    return 0;

/*           PROCEDURE FOR A AND B .LT. 1.E-3*EPS */

L260:
    *w = *b / (*a + *b);
    *w1 = *a / (*a + *b);
    return 0;

/*                       ERROR RETURN */

L270:
    *ierr = 1;
    return 0;
L280:
    *ierr = 2;
    return 0;
L290:
    *ierr = 3;
    return 0;
L300:
    *ierr = 4;
    return 0;
L310:
    *ierr = 5;
    return 0;
L320:
    *ierr = 6;
    return 0;
L330:
    *ierr = 7;
    return 0;
} /* bratio_ */

