/* gaminv.f -- translated by f2c (version 20190311).
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
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__0 = 0;

/* Subroutine */ int gaminv_(doublereal *a, doublereal *x, doublereal *x0, 
	doublereal *p, doublereal *q, integer *ierr)
{
    /* Initialized data */

    static doublereal ln10 = 2.302585;
    static doublereal b4 = .036117081018842;
    static doublereal eps0[2] = { 1e-10,1e-8 };
    static doublereal amin[2] = { 500.,100. };
    static doublereal bmin[2] = { 1e-28,1e-13 };
    static doublereal dmin__[2] = { 1e-6,1e-4 };
    static doublereal emin[2] = { .002,.006 };
    static doublereal tol = 1e-5;
    static doublereal c__ = .577215664901533;
    static doublereal a0 = 3.31125922108741;
    static doublereal a1 = 11.6616720288968;
    static doublereal a2 = 4.28342155967104;
    static doublereal a3 = .213623493715853;
    static doublereal b1 = 6.61053765625462;
    static doublereal b2 = 6.40691597760039;
    static doublereal b3 = 1.27364489782223;

    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double exp(doublereal), log(doublereal), sqrt(doublereal);

    /* Local variables */
    static doublereal b, d__, e, g, h__, r__, s, t, u, w, y, z__, c1, c2, c3, 
	    c4, c5, e2, s2, qg, pn, qn, xn, am1, ap1, ap2, ap3, apn, rta, eps;
    static integer iop;
    static doublereal sum, amax, xmin, xmax;
    extern doublereal gamma_(doublereal *), gamln_(doublereal *), rcomp_(
	    doublereal *, doublereal *), gamln1_(doublereal *), alnrel_(
	    doublereal *);
    extern /* Subroutine */ int gratio_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *);
    extern doublereal spmpar_(integer *);

/* ---------------------------------------------------------------------- */
/*            INVERSE INCOMPLETE GAMMA RATIO FUNCTION */

/*     GIVEN POSITIVE A, AND NONEGATIVE P AND Q WHERE P + Q = 1. */
/*     THEN X IS COMPUTED WHERE P(A,X) = P AND Q(A,X) = Q. SCHRODER */
/*     ITERATION IS EMPLOYED. THE ROUTINE ATTEMPTS TO COMPUTE X */
/*     TO 10 SIGNIFICANT DIGITS IF THIS IS POSSIBLE FOR THE */
/*     PARTICULAR COMPUTER ARITHMETIC BEING USED. */

/*                      ------------ */

/*     X IS A VARIABLE. IF P = 0 THEN X IS ASSIGNED THE VALUE 0, */
/*     AND IF Q = 0 THEN X IS SET TO THE LARGEST FLOATING POINT */
/*     NUMBER AVAILABLE. OTHERWISE, GAMINV ATTEMPTS TO OBTAIN */
/*     A SOLUTION FOR P(A,X) = P AND Q(A,X) = Q. IF THE ROUTINE */
/*     IS SUCCESSFUL THEN THE SOLUTION IS STORED IN X. */

/*     X0 IS AN OPTIONAL INITIAL APPROXIMATION FOR X. IF THE USER */
/*     DOES NOT WISH TO SUPPLY AN INITIAL APPROXIMATION, THEN SET */
/*     X0 .LE. 0. */

/*     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS. */
/*     WHEN THE ROUTINE TERMINATES, IERR HAS ONE OF THE FOLLOWING */
/*     VALUES ... */

/*       IERR =  0    THE SOLUTION WAS OBTAINED. ITERATION WAS */
/*                    NOT USED. */
/*       IERR.GT.0    THE SOLUTION WAS OBTAINED. IERR ITERATIONS */
/*                    WERE PERFORMED. */
/*       IERR = -2    (INPUT ERROR) A .LE. 0 */
/*       IERR = -3    NO SOLUTION WAS OBTAINED. THE RATIO Q/A */
/*                    IS TOO LARGE. */
/*       IERR = -4    (INPUT ERROR) P + Q .NE. 1 */
/*       IERR = -6    20 ITERATIONS WERE PERFORMED. THE MOST */
/*                    RECENT VALUE OBTAINED FOR X IS GIVEN. */
/*                    THIS CANNOT OCCUR IF X0 .LE. 0. */
/*       IERR = -7    ITERATION FAILED. NO VALUE IS GIVEN FOR X. */
/*                    THIS MAY OCCUR WHEN X IS APPROXIMATELY 0. */
/*       IERR = -8    A VALUE FOR X HAS BEEN OBTAINED, BUT THE */
/*                    ROUTINE IS NOT CERTAIN OF ITS ACCURACY. */
/*                    ITERATION CANNOT BE PERFORMED IN THIS */
/*                    CASE. IF X0 .LE. 0, THIS CAN OCCUR ONLY */
/*                    WHEN P OR Q IS APPROXIMATELY 0. IF X0 IS */
/*                    POSITIVE THEN THIS CAN OCCUR WHEN A IS */
/*                    EXCEEDINGLY CLOSE TO X AND A IS EXTREMELY */
/*                    LARGE (SAY A .GE. 1.E20). */
/* ---------------------------------------------------------------------- */
/*     WRITTEN BY ALFRED H. MORRIS, JR. */
/*        NAVAL SURFACE WEAPONS CENTER */
/*        DAHLGREN, VIRGINIA */
/*     ------------------- */
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
/*     .. Data statements .. */
/*     ------------------- */
/*     LN10 = LN(10) */
/*     C = EULER CONSTANT */
/*     ------------------- */
/*     ------------------- */
/*     ------------------- */
/*     ------------------- */
/*     .. */
/*     .. Executable Statements .. */
/*     ------------------- */
/*     ****** E, XMIN, AND XMAX ARE MACHINE DEPENDENT CONSTANTS. */
/*            E IS THE SMALLEST NUMBER FOR WHICH 1.0 + E .GT. 1.0. */
/*            XMIN IS THE SMALLEST POSITIVE NUMBER AND XMAX IS THE */
/*            LARGEST POSITIVE NUMBER. */

    e = spmpar_(&c__1);
    xmin = spmpar_(&c__2);
    xmax = spmpar_(&c__3);
/*     ------------------- */
    *x = 0.;
    if (*a <= 0.) {
	goto L300;
    }
    t = *p + *q - 1.;
    if (abs(t) > e) {
	goto L320;
    }

    *ierr = 0;
    if (*p == 0.) {
	return 0;
    }
    if (*q == 0.) {
	goto L270;
    }
    if (*a == 1.) {
	goto L280;
    }

    e2 = e * 2.;
    amax = 4e-11 / (e * e);
    iop = 1;
    if (e > 1e-10) {
	iop = 2;
    }
    eps = eps0[iop - 1];
    xn = *x0;
    if (*x0 > 0.) {
	goto L160;
    }

/*        SELECTION OF THE INITIAL APPROXIMATION XN OF X */
/*                       WHEN A .LT. 1 */

    if (*a > 1.) {
	goto L80;
    }
    d__1 = *a + 1.;
    g = gamma_(&d__1);
    qg = *q * g;
    if (qg == 0.) {
	goto L360;
    }
    b = qg / *a;
    if (qg > *a * .6) {
	goto L40;
    }
    if (*a >= .3 || b < .35) {
	goto L10;
    }
    t = exp(-(b + c__));
    u = t * exp(t);
    xn = t * exp(u);
    goto L160;

L10:
    if (b >= .45) {
	goto L40;
    }
    if (b == 0.) {
	goto L360;
    }
    y = -log(b);
    s = .5 - *a + .5;
    z__ = log(y);
    t = y - s * z__;
    if (b < .15) {
	goto L20;
    }
    xn = y - s * log(t) - log(s / (t + 1.) + 1.);
    goto L220;
L20:
    if (b <= .01) {
	goto L30;
    }
    u = ((t + (3. - *a) * 2.) * t + (2. - *a) * (3. - *a)) / ((t + (5. - *a)) 
	    * t + 2.);
    xn = y - s * log(t) - log(u);
    goto L220;
L30:
    c1 = -s * z__;
    c2 = -s * (c1 + 1.);
    c3 = s * ((c1 * .5 + (2. - *a)) * c1 + (2.5 - *a * 1.5));
    c4 = -s * (((c1 / 3. + (2.5 - *a * 1.5)) * c1 + ((*a - 6.) * *a + 7.)) * 
	    c1 + ((*a * 11. - 46) * *a + 47.) / 6.);
    c5 = -s * ((((-c1 / 4. + (*a * 11. - 17.) / 6.) * c1 + ((*a * -3. + 13.) *
	     *a - 13.)) * c1 + (((*a * 2. - 25.) * *a + 72.) * *a - 61.) * .5)
	     * c1 + (((*a * 25. - 195.) * *a + 477.) * *a - 379.) / 12.);
    xn = (((c5 / y + c4) / y + c3) / y + c2) / y + c1 + y;
    if (*a > 1.) {
	goto L220;
    }
    if (b > bmin[iop - 1]) {
	goto L220;
    }
    *x = xn;
    return 0;

L40:
    if (b * *q > 1e-8) {
	goto L50;
    }
    xn = exp(-(*q / *a + c__));
    goto L70;
L50:
    if (*p <= .9) {
	goto L60;
    }
    d__1 = -(*q);
    xn = exp((alnrel_(&d__1) + gamln1_(a)) / *a);
    goto L70;
L60:
    xn = exp(log(*p * g) / *a);
L70:
    if (xn == 0.) {
	goto L310;
    }
    t = .5 - xn / (*a + 1.) + .5;
    xn /= t;
    goto L160;

/*        SELECTION OF THE INITIAL APPROXIMATION XN OF X */
/*                       WHEN A .GT. 1 */

L80:
    if (*q <= .5) {
	goto L90;
    }
    w = log(*p);
    goto L100;
L90:
    w = log(*q);
L100:
    t = sqrt(w * -2.);
    s = t - (((a3 * t + a2) * t + a1) * t + a0) / ((((b4 * t + b3) * t + b2) *
	     t + b1) * t + 1.);
    if (*q > .5) {
	s = -s;
    }

    rta = sqrt(*a);
    s2 = s * s;
    xn = *a + s * rta + (s2 - 1.) / 3. + s * (s2 - 7.) / (rta * 36.) - ((s2 * 
	    3. + 7.) * s2 - 16.) / (*a * 810.) + s * ((s2 * 9. + 256.) * s2 - 
	    433.) / (*a * 38880. * rta);
    xn = max(xn,0.);
    if (*a < amin[iop - 1]) {
	goto L110;
    }
    *x = xn;
    d__ = .5 - *x / *a + .5;
    if (abs(d__) <= dmin__[iop - 1]) {
	return 0;
    }

L110:
    if (*p <= .5) {
	goto L130;
    }
    if (xn < *a * 3.) {
	goto L220;
    }
    y = -(w + gamln_(a));
/* Computing MAX */
    d__1 = 2., d__2 = *a * (*a - 1.);
    d__ = max(d__1,d__2);
    if (y < ln10 * d__) {
	goto L120;
    }
    s = 1. - *a;
    z__ = log(y);
    goto L30;
L120:
    t = *a - 1.;
    d__1 = -t / (xn + 1.);
    xn = y + t * log(xn) - alnrel_(&d__1);
    d__1 = -t / (xn + 1.);
    xn = y + t * log(xn) - alnrel_(&d__1);
    goto L220;

L130:
    ap1 = *a + 1.;
    if (xn > ap1 * .7) {
	goto L170;
    }
    w += gamln_(&ap1);
    if (xn > ap1 * .15) {
	goto L140;
    }
    ap2 = *a + 2.;
    ap3 = *a + 3.;
    *x = exp((w + *x) / *a);
    *x = exp((w + *x - log(*x / ap1 * (*x / ap2 + 1.) + 1.)) / *a);
    *x = exp((w + *x - log(*x / ap1 * (*x / ap2 + 1.) + 1.)) / *a);
    *x = exp((w + *x - log(*x / ap1 * (*x / ap2 * (*x / ap3 + 1.) + 1.) + 1.))
	     / *a);
    xn = *x;
    if (xn > ap1 * .01) {
	goto L140;
    }
    if (xn <= emin[iop - 1] * ap1) {
	return 0;
    }
    goto L170;

L140:
    apn = ap1;
    t = xn / apn;
    sum = t + 1.;
L150:
    apn += 1.;
    t *= xn / apn;
    sum += t;
    if (t > 1e-4) {
	goto L150;
    }
    t = w - log(sum);
    xn = exp((xn + t) / *a);
    xn *= 1. - (*a * log(xn) - xn - t) / (*a - xn);
    goto L170;

/*                 SCHRODER ITERATION USING P */

L160:
    if (*p > .5) {
	goto L220;
    }
L170:
    if (*p <= xmin * 1e10) {
	goto L350;
    }
    am1 = *a - .5 - .5;
L180:
    if (*a <= amax) {
	goto L190;
    }
    d__ = .5 - xn / *a + .5;
    if (abs(d__) <= e2) {
	goto L350;
    }

L190:
    if (*ierr >= 20) {
	goto L330;
    }
    ++(*ierr);
    gratio_(a, &xn, &pn, &qn, &c__0);
    if (pn == 0. || qn == 0.) {
	goto L350;
    }
    r__ = rcomp_(a, &xn);
    if (r__ == 0.) {
	goto L350;
    }
    t = (pn - *p) / r__;
    w = (am1 - xn) * .5;
    if (abs(t) <= .1 && (d__1 = w * t, abs(d__1)) <= .1) {
	goto L200;
    }
    *x = xn * (1. - t);
    if (*x <= 0.) {
	goto L340;
    }
    d__ = abs(t);
    goto L210;

L200:
    h__ = t * (w * t + 1.);
    *x = xn * (1. - h__);
    if (*x <= 0.) {
	goto L340;
    }
    if (abs(w) >= 1. && abs(w) * t * t <= eps) {
	return 0;
    }
    d__ = abs(h__);
L210:
    xn = *x;
    if (d__ > tol) {
	goto L180;
    }
    if (d__ <= eps) {
	return 0;
    }
    if ((d__1 = *p - pn, abs(d__1)) <= tol * *p) {
	return 0;
    }
    goto L180;

/*                 SCHRODER ITERATION USING Q */

L220:
    if (*q <= xmin * 1e10) {
	goto L350;
    }
    am1 = *a - .5 - .5;
L230:
    if (*a <= amax) {
	goto L240;
    }
    d__ = .5 - xn / *a + .5;
    if (abs(d__) <= e2) {
	goto L350;
    }

L240:
    if (*ierr >= 20) {
	goto L330;
    }
    ++(*ierr);
    gratio_(a, &xn, &pn, &qn, &c__0);
    if (pn == 0. || qn == 0.) {
	goto L350;
    }
    r__ = rcomp_(a, &xn);
    if (r__ == 0.) {
	goto L350;
    }
    t = (*q - qn) / r__;
    w = (am1 - xn) * .5;
    if (abs(t) <= .1 && (d__1 = w * t, abs(d__1)) <= .1) {
	goto L250;
    }
    *x = xn * (1. - t);
    if (*x <= 0.) {
	goto L340;
    }
    d__ = abs(t);
    goto L260;

L250:
    h__ = t * (w * t + 1.);
    *x = xn * (1. - h__);
    if (*x <= 0.) {
	goto L340;
    }
    if (abs(w) >= 1. && abs(w) * t * t <= eps) {
	return 0;
    }
    d__ = abs(h__);
L260:
    xn = *x;
    if (d__ > tol) {
	goto L230;
    }
    if (d__ <= eps) {
	return 0;
    }
    if ((d__1 = *q - qn, abs(d__1)) <= tol * *q) {
	return 0;
    }
    goto L230;

/*                       SPECIAL CASES */

L270:
    *x = xmax;
    return 0;

L280:
    if (*q < .9) {
	goto L290;
    }
    d__1 = -(*p);
    *x = -alnrel_(&d__1);
    return 0;
L290:
    *x = -log(*q);
    return 0;

/*                       ERROR RETURN */

L300:
    *ierr = -2;
    return 0;

L310:
    *ierr = -3;
    return 0;

L320:
    *ierr = -4;
    return 0;

L330:
    *ierr = -6;
    return 0;

L340:
    *ierr = -7;
    return 0;

L350:
    *x = xn;
    *ierr = -8;
    return 0;

L360:
    *x = xmax;
    *ierr = -8;
    return 0;
} /* gaminv_ */

