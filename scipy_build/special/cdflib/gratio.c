/* gratio.f -- translated by f2c (version 20190311).
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
static integer c__0 = 0;

/* Subroutine */ int gratio_(doublereal *a, doublereal *x, doublereal *ans, 
	doublereal *qans, integer *ind)
{
    /* Initialized data */

    static doublereal acc0[3] = { 5e-15,5e-7,5e-4 };
    static doublereal d10 = -.00185185185185185;
    static doublereal d1[12] = { -.00347222222222222,.00264550264550265,
	    -9.9022633744856e-4,2.05761316872428e-4,-4.01877572016461e-7,
	    -1.809855033449e-5,7.64916091608111e-6,-1.61209008945634e-6,
	    4.64712780280743e-9,1.37863344691572e-7,-5.7525456035177e-8,
	    1.19516285997781e-8 };
    static doublereal d20 = .00413359788359788;
    static doublereal d2[10] = { -.00268132716049383,7.71604938271605e-4,
	    2.0093878600823e-6,-1.07366532263652e-4,5.29234488291201e-5,
	    -1.27606351886187e-5,3.42357873409614e-8,1.37219573090629e-6,
	    -6.29899213838006e-7,1.42806142060642e-7 };
    static doublereal d30 = 6.49434156378601e-4;
    static doublereal d3[8] = { 2.29472093621399e-4,-4.69189494395256e-4,
	    2.67720632062839e-4,-7.56180167188398e-5,-2.3965051138673e-7,
	    1.10826541153473e-5,-5.6749528269916e-6,1.42309007324359e-6 };
    static doublereal d40 = -8.61888290916712e-4;
    static doublereal d4[6] = { 7.84039221720067e-4,-2.9907248030319e-4,
	    -1.46384525788434e-6,6.64149821546512e-5,-3.96836504717943e-5,
	    1.13757269706784e-5 };
    static doublereal d50 = -3.36798553366358e-4;
    static doublereal d5[4] = { -6.97281375836586e-5,2.77275324495939e-4,
	    -1.99325705161888e-4,6.79778047793721e-5 };
    static doublereal big[3] = { 20.,14.,10. };
    static doublereal d60 = 5.31307936463992e-4;
    static doublereal d6[2] = { -5.92166437353694e-4,2.70878209671804e-4 };
    static doublereal d70 = 3.44367606892378e-4;
    static doublereal e00[3] = { 2.5e-4,.025,.14 };
    static doublereal x00[3] = { 31.,17.,9.7 };
    static doublereal alog10 = 2.30258509299405;
    static doublereal rt2pin = .398942280401433;
    static doublereal rtpi = 1.77245385090552;
    static doublereal third = .333333333333333;
    static doublereal d0[13] = { .0833333333333333,-.0148148148148148,
	    .00115740740740741,3.52733686067019e-4,-1.78755144032922e-4,
	    3.91926317852244e-5,-2.18544851067999e-6,-1.85406221071516e-6,
	    8.29671134095309e-7,-1.76659527368261e-7,6.7078535434015e-9,
	    1.02618097842403e-8,-4.38203601845335e-9 };

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double log(doublereal), exp(doublereal), sqrt(doublereal);

    /* Local variables */
    static doublereal c__, e, g, h__;
    static integer i__;
    static doublereal j, l;
    static integer m, n;
    static doublereal r__, s, t, u, w, y, z__, c0, c1, c2, c3, c4, c5, c6, e0,
	     t1, x0, an, wk[20], am0, an0, a2n, b2n, acc, cma, amn;
    extern doublereal erf_(doublereal *);
    static doublereal apn;
    static integer max__;
    static doublereal rta;
    static integer iop;
    static doublereal tol, sum, rtx;
    extern doublereal gam1_(doublereal *);
    static doublereal a2nm1, b2nm1;
    extern doublereal rlog_(doublereal *);
    static doublereal twoa;
    extern doublereal rexp_(doublereal *), erfc1_(integer *, doublereal *), 
	    gamma_(doublereal *), spmpar_(integer *);

/* ---------------------------------------------------------------------- */
/*        EVALUATION OF THE INCOMPLETE GAMMA RATIO FUNCTIONS */
/*                      P(A,X) AND Q(A,X) */

/*                        ---------- */

/*     IT IS ASSUMED THAT A AND X ARE NONNEGATIVE, WHERE A AND X */
/*     ARE NOT BOTH 0. */

/*     ANS AND QANS ARE VARIABLES. GRATIO ASSIGNS ANS THE VALUE */
/*     P(A,X) AND QANS THE VALUE Q(A,X). IND MAY BE ANY INTEGER. */
/*     IF IND = 0 THEN THE USER IS REQUESTING AS MUCH ACCURACY AS */
/*     POSSIBLE (UP TO 14 SIGNIFICANT DIGITS). OTHERWISE, IF */
/*     IND = 1 THEN ACCURACY IS REQUESTED TO WITHIN 1 UNIT OF THE */
/*     6-TH SIGNIFICANT DIGIT, AND IF IND .NE. 0,1 THEN ACCURACY */
/*     IS REQUESTED TO WITHIN 1 UNIT OF THE 3RD SIGNIFICANT DIGIT. */

/*     ERROR RETURN ... */
/*        ANS IS ASSIGNED THE VALUE 2 WHEN A OR X IS NEGATIVE, */
/*     WHEN A*X = 0, OR WHEN P(A,X) AND Q(A,X) ARE INDETERMINANT. */
/*     P(A,X) AND Q(A,X) ARE COMPUTATIONALLY INDETERMINANT WHEN */
/*     X IS EXCEEDINGLY CLOSE TO A AND A IS EXTREMELY LARGE. */
/* ---------------------------------------------------------------------- */
/*     WRITTEN BY ALFRED H. MORRIS, JR. */
/*        NAVAL SURFACE WEAPONS CENTER */
/*        DAHLGREN, VIRGINIA */
/*     -------------------- */
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
/*     -------------------- */
/*     -------------------- */
/*     ALOG10 = LN(10) */
/*     RT2PIN = 1/SQRT(2*PI) */
/*     RTPI   = SQRT(PI) */
/*     -------------------- */
/*     -------------------- */
/*     -------------------- */
/*     -------------------- */
/*     -------------------- */
/*     -------------------- */
/*     -------------------- */
/*     -------------------- */
/*     -------------------- */
/*     .. */
/*     .. Executable Statements .. */
/*     -------------------- */
/*     ****** E IS A MACHINE DEPENDENT CONSTANT. E IS THE SMALLEST */
/*            FLOATING POINT NUMBER FOR WHICH 1.0 + E .GT. 1.0 . */

    e = spmpar_(&c__1);

/*     -------------------- */
    if (*a < 0. || *x < 0.) {
	goto L430;
    }
    if (*a == 0. && *x == 0.) {
	goto L430;
    }
    if (*a * *x == 0.) {
	goto L420;
    }

    iop = *ind + 1;
    if (iop != 1 && iop != 2) {
	iop = 3;
    }
/* Computing MAX */
    d__1 = acc0[iop - 1];
    acc = max(d__1,e);
    e0 = e00[iop - 1];
    x0 = x00[iop - 1];

/*            SELECT THE APPROPRIATE ALGORITHM */

    if (*a >= 1.) {
	goto L10;
    }
    if (*a == .5) {
	goto L390;
    }
    if (*x < 1.1) {
	goto L160;
    }
    t1 = *a * log(*x) - *x;
    u = *a * exp(t1);
    if (u == 0.) {
	goto L380;
    }
    r__ = u * (gam1_(a) + 1.);
    goto L250;

L10:
    if (*a >= big[iop - 1]) {
	goto L30;
    }
    if (*a > *x || *x >= x0) {
	goto L20;
    }
    twoa = *a + *a;
    m = (integer) twoa;
    if (twoa != (doublereal) m) {
	goto L20;
    }
    i__ = m / 2;
    if (*a == (doublereal) i__) {
	goto L210;
    }
    goto L220;
L20:
    t1 = *a * log(*x) - *x;
    r__ = exp(t1) / gamma_(a);
    goto L40;

L30:
    l = *x / *a;
    if (l == 0.) {
	goto L370;
    }
    s = .5 - l + .5;
    z__ = rlog_(&l);
    if (z__ >= 700. / *a) {
	goto L410;
    }
    y = *a * z__;
    rta = sqrt(*a);
    if (abs(s) <= e0 / rta) {
	goto L330;
    }
    if (abs(s) <= .4) {
	goto L270;
    }

/* Computing 2nd power */
    d__1 = 1. / *a;
    t = d__1 * d__1;
    t1 = (((t * .75 - 1.) * t + 3.5) * t - 105.) / (*a * 1260.);
    t1 -= y;
    r__ = rt2pin * rta * exp(t1);

L40:
    if (r__ == 0.) {
	goto L420;
    }
    if (*x <= max(*a,alog10)) {
	goto L50;
    }
    if (*x < x0) {
	goto L250;
    }
    goto L100;

/*                 TAYLOR SERIES FOR P/R */

L50:
    apn = *a + 1.;
    t = *x / apn;
    wk[0] = t;
    for (n = 2; n <= 20; ++n) {
	apn += 1.;
	t *= *x / apn;
	if (t <= .001) {
	    goto L70;
	}
	wk[n - 1] = t;
/* L60: */
    }
    n = 20;

L70:
    sum = t;
    tol = acc * .5;
L80:
    apn += 1.;
    t *= *x / apn;
    sum += t;
    if (t > tol) {
	goto L80;
    }

    max__ = n - 1;
    i__1 = max__;
    for (m = 1; m <= i__1; ++m) {
	--n;
	sum += wk[n - 1];
/* L90: */
    }
    *ans = r__ / *a * (sum + 1.);
    *qans = .5 - *ans + .5;
    return 0;

/*                 ASYMPTOTIC EXPANSION */

L100:
    amn = *a - 1.;
    t = amn / *x;
    wk[0] = t;
    for (n = 2; n <= 20; ++n) {
	amn += -1.;
	t *= amn / *x;
	if (abs(t) <= .001) {
	    goto L120;
	}
	wk[n - 1] = t;
/* L110: */
    }
    n = 20;

L120:
    sum = t;
L130:
    if (! (abs(t) > acc)) {
	goto L140;
    }
    amn += -1.;
    t *= amn / *x;
    sum += t;
    goto L130;

L140:
    max__ = n - 1;
    i__1 = max__;
    for (m = 1; m <= i__1; ++m) {
	--n;
	sum += wk[n - 1];
/* L150: */
    }
    *qans = r__ / *x * (sum + 1.);
    *ans = .5 - *qans + .5;
    return 0;

/*             TAYLOR SERIES FOR P(A,X)/X**A */

L160:
    an = 3.;
    c__ = *x;
    sum = *x / (*a + 3.);
    tol = acc * 3. / (*a + 1.);
L170:
    an += 1.;
    c__ = -c__ * (*x / an);
    t = c__ / (*a + an);
    sum += t;
    if (abs(t) > tol) {
	goto L170;
    }
    j = *a * *x * ((sum / 6. - .5 / (*a + 2.)) * *x + 1. / (*a + 1.));

    z__ = *a * log(*x);
    h__ = gam1_(a);
    g = h__ + 1.;
    if (*x < .25) {
	goto L180;
    }
    if (*a < *x / 2.59) {
	goto L200;
    }
    goto L190;
L180:
    if (z__ > -.13394) {
	goto L200;
    }

L190:
    w = exp(z__);
    *ans = w * g * (.5 - j + .5);
    *qans = .5 - *ans + .5;
    return 0;

L200:
    l = rexp_(&z__);
    w = l + .5 + .5;
    *qans = (w * j - l) * g - h__;
    if (*qans < 0.) {
	goto L380;
    }
    *ans = .5 - *qans + .5;
    return 0;

/*             FINITE SUMS FOR Q WHEN A .GE. 1 */
/*                 AND 2*A IS AN INTEGER */

L210:
    sum = exp(-(*x));
    t = sum;
    n = 1;
    c__ = 0.;
    goto L230;

L220:
    rtx = sqrt(*x);
    sum = erfc1_(&c__0, &rtx);
    t = exp(-(*x)) / (rtpi * rtx);
    n = 0;
    c__ = -.5;

L230:
    if (n == i__) {
	goto L240;
    }
    ++n;
    c__ += 1.;
    t = *x * t / c__;
    sum += t;
    goto L230;
L240:
    *qans = sum;
    *ans = .5 - *qans + .5;
    return 0;

/*              CONTINUED FRACTION EXPANSION */

L250:
/* Computing MAX */
    d__1 = e * 5.;
    tol = max(d__1,acc);
    a2nm1 = 1.;
    a2n = 1.;
    b2nm1 = *x;
    b2n = *x + (1. - *a);
    c__ = 1.;
L260:
    a2nm1 = *x * a2n + c__ * a2nm1;
    b2nm1 = *x * b2n + c__ * b2nm1;
    am0 = a2nm1 / b2nm1;
    c__ += 1.;
    cma = c__ - *a;
    a2n = a2nm1 + cma * a2n;
    b2n = b2nm1 + cma * b2n;
    an0 = a2n / b2n;
    if ((d__1 = an0 - am0, abs(d__1)) >= tol * an0) {
	goto L260;
    }

    *qans = r__ * an0;
    *ans = .5 - *qans + .5;
    return 0;

/*                GENERAL TEMME EXPANSION */

L270:
    if (abs(s) <= e * 2. && *a * e * e > .00328) {
	goto L430;
    }
    c__ = exp(-y);
    d__1 = sqrt(y);
    w = erfc1_(&c__1, &d__1) * .5;
    u = 1. / *a;
    z__ = sqrt(z__ + z__);
    if (l < 1.) {
	z__ = -z__;
    }
    if (iop < 2) {
	goto L280;
    }
    if (iop == 2) {
	goto L290;
    }
    goto L300;

L280:
    if (abs(s) <= .001) {
	goto L340;
    }
    c0 = ((((((((((((d0[12] * z__ + d0[11]) * z__ + d0[10]) * z__ + d0[9]) * 
	    z__ + d0[8]) * z__ + d0[7]) * z__ + d0[6]) * z__ + d0[5]) * z__ + 
	    d0[4]) * z__ + d0[3]) * z__ + d0[2]) * z__ + d0[1]) * z__ + d0[0])
	     * z__ - third;
    c1 = (((((((((((d1[11] * z__ + d1[10]) * z__ + d1[9]) * z__ + d1[8]) * 
	    z__ + d1[7]) * z__ + d1[6]) * z__ + d1[5]) * z__ + d1[4]) * z__ + 
	    d1[3]) * z__ + d1[2]) * z__ + d1[1]) * z__ + d1[0]) * z__ + d10;
    c2 = (((((((((d2[9] * z__ + d2[8]) * z__ + d2[7]) * z__ + d2[6]) * z__ + 
	    d2[5]) * z__ + d2[4]) * z__ + d2[3]) * z__ + d2[2]) * z__ + d2[1])
	     * z__ + d2[0]) * z__ + d20;
    c3 = (((((((d3[7] * z__ + d3[6]) * z__ + d3[5]) * z__ + d3[4]) * z__ + d3[
	    3]) * z__ + d3[2]) * z__ + d3[1]) * z__ + d3[0]) * z__ + d30;
    c4 = (((((d4[5] * z__ + d4[4]) * z__ + d4[3]) * z__ + d4[2]) * z__ + d4[1]
	    ) * z__ + d4[0]) * z__ + d40;
    c5 = (((d5[3] * z__ + d5[2]) * z__ + d5[1]) * z__ + d5[0]) * z__ + d50;
    c6 = (d6[1] * z__ + d6[0]) * z__ + d60;
    t = ((((((d70 * u + c6) * u + c5) * u + c4) * u + c3) * u + c2) * u + c1) 
	    * u + c0;
    goto L310;

L290:
    c0 = (((((d0[5] * z__ + d0[4]) * z__ + d0[3]) * z__ + d0[2]) * z__ + d0[1]
	    ) * z__ + d0[0]) * z__ - third;
    c1 = (((d1[3] * z__ + d1[2]) * z__ + d1[1]) * z__ + d1[0]) * z__ + d10;
    c2 = d2[0] * z__ + d20;
    t = (c2 * u + c1) * u + c0;
    goto L310;

L300:
    t = ((d0[2] * z__ + d0[1]) * z__ + d0[0]) * z__ - third;

L310:
    if (l < 1.) {
	goto L320;
    }
    *qans = c__ * (w + rt2pin * t / rta);
    *ans = .5 - *qans + .5;
    return 0;
L320:
    *ans = c__ * (w - rt2pin * t / rta);
    *qans = .5 - *ans + .5;
    return 0;

/*               TEMME EXPANSION FOR L = 1 */

L330:
    if (*a * e * e > .00328) {
	goto L430;
    }
    c__ = .5 - y + .5;
    w = (.5 - sqrt(y) * (.5 - y / 3. + .5) / rtpi) / c__;
    u = 1. / *a;
    z__ = sqrt(z__ + z__);
    if (l < 1.) {
	z__ = -z__;
    }
    if (iop < 2) {
	goto L340;
    }
    if (iop == 2) {
	goto L350;
    }
    goto L360;

L340:
    c0 = ((((((d0[6] * z__ + d0[5]) * z__ + d0[4]) * z__ + d0[3]) * z__ + d0[
	    2]) * z__ + d0[1]) * z__ + d0[0]) * z__ - third;
    c1 = (((((d1[5] * z__ + d1[4]) * z__ + d1[3]) * z__ + d1[2]) * z__ + d1[1]
	    ) * z__ + d1[0]) * z__ + d10;
    c2 = ((((d2[4] * z__ + d2[3]) * z__ + d2[2]) * z__ + d2[1]) * z__ + d2[0])
	     * z__ + d20;
    c3 = (((d3[3] * z__ + d3[2]) * z__ + d3[1]) * z__ + d3[0]) * z__ + d30;
    c4 = (d4[1] * z__ + d4[0]) * z__ + d40;
    c5 = (d5[1] * z__ + d5[0]) * z__ + d50;
    c6 = d6[0] * z__ + d60;
    t = ((((((d70 * u + c6) * u + c5) * u + c4) * u + c3) * u + c2) * u + c1) 
	    * u + c0;
    goto L310;

L350:
    c0 = (d0[1] * z__ + d0[0]) * z__ - third;
    c1 = d1[0] * z__ + d10;
    t = (d20 * u + c1) * u + c0;
    goto L310;

L360:
    t = d0[0] * z__ - third;
    goto L310;

/*                     SPECIAL CASES */

L370:
    *ans = 0.;
    *qans = 1.;
    return 0;

L380:
    *ans = 1.;
    *qans = 0.;
    return 0;

L390:
    if (*x >= .25) {
	goto L400;
    }
    d__1 = sqrt(*x);
    *ans = erf_(&d__1);
    *qans = .5 - *ans + .5;
    return 0;
L400:
    d__1 = sqrt(*x);
    *qans = erfc1_(&c__0, &d__1);
    *ans = .5 - *qans + .5;
    return 0;

L410:
    if (abs(s) <= e * 2.) {
	goto L430;
    }
L420:
    if (*x <= *a) {
	goto L370;
    }
    goto L380;

/*                     ERROR RETURN */

L430:
    *ans = 2.;
    return 0;
} /* gratio_ */

