/* cdfbet.f -- translated by f2c (version 20190311).
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
static doublereal c_b35 = 0.;
static doublereal c_b36 = 1.;
static doublereal c_b37 = 1e-50;
static doublereal c_b38 = 1e-8;
static doublereal c_b48 = 1e-100;
static doublereal c_b49 = 1e100;
static doublereal c_b50 = .5;
static doublereal c_b52 = 5.;

/* Subroutine */ int cdfbet_(integer *which, doublereal *p, doublereal *q, 
	doublereal *x, doublereal *y, doublereal *a, doublereal *b, integer *
	status, doublereal *bound)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static doublereal fx, pq, xy;
    static logical qhi;
    static doublereal cum, xhi, xlo, ccum;
    static logical qleft;
    extern /* Subroutine */ int dinvr_(integer *, doublereal *, doublereal *, 
	    logical *, logical *), dzror_(integer *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, logical *, logical *);
    static logical qporq;
    extern /* Subroutine */ int dstzr_(doublereal *, doublereal *, doublereal 
	    *, doublereal *), cumbet_(doublereal *, doublereal *, doublereal *
	    , doublereal *, doublereal *, doublereal *);
    extern doublereal spmpar_(integer *);
    extern /* Subroutine */ int dstinv_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);

/* ********************************************************************** */

/*      SUBROUTINE CDFBET( WHICH, P, Q, X, Y, A, B, STATUS, BOUND ) */
/*               Cumulative Distribution Function */
/*                         BETa Distribution */


/*                              Function */


/*     Calculates any one parameter of the beta distribution given */
/*     values for the others. */


/*                              Arguments */


/*     WHICH --> Integer indicating which of the next four argument */
/*               values is to be calculated from the others. */
/*               Legal range: 1..4 */
/*               iwhich = 1 : Calculate P and Q from X,Y,A and B */
/*               iwhich = 2 : Calculate X and Y from P,Q,A and B */
/*               iwhich = 3 : Calculate A from P,Q,X,Y and B */
/*               iwhich = 4 : Calculate B from P,Q,X,Y and A */

/*                    INTEGER WHICH */

/*     P <--> The integral from 0 to X of the chi-square */
/*            distribution. */
/*            Input range: [0, 1]. */
/*                    DOUBLE PRECISION P */

/*     Q <--> 1-P. */
/*            Input range: [0, 1]. */
/*            P + Q = 1.0. */
/*                    DOUBLE PRECISION Q */

/*     X <--> Upper limit of integration of beta density. */
/*            Input range: [0,1]. */
/*            Search range: [0,1] */
/*                    DOUBLE PRECISION X */

/*     Y <--> 1-X. */
/*            Input range: [0,1]. */
/*            Search range: [0,1] */
/*            X + Y = 1.0. */
/*                    DOUBLE PRECISION Y */

/*     A <--> The first parameter of the beta density. */
/*            Input range: (0, +infinity). */
/*            Search range: [1D-100,1D100] */
/*                    DOUBLE PRECISION A */

/*     B <--> The second parameter of the beta density. */
/*            Input range: (0, +infinity). */
/*            Search range: [1D-100,1D100] */
/*                    DOUBLE PRECISION B */

/*     STATUS <-- 0 if calculation completed correctly */
/*               -I if input parameter number I is out of range */
/*                1 if answer appears to be lower than lowest */
/*                  search bound */
/*                2 if answer appears to be higher than greatest */
/*                  search bound */
/*                3 if P + Q .ne. 1 */
/*                4 if X + Y .ne. 1 */
/*                    INTEGER STATUS */

/*     BOUND <-- Undefined if STATUS is 0 */

/*               Bound exceeded by parameter number I if STATUS */
/*               is negative. */

/*               Lower search bound if STATUS is 1. */

/*               Upper search bound if STATUS is 2. */


/*                              Method */


/*     Cumulative distribution function  (P)  is calculated directly by */
/*     code associated with the following reference. */

/*     DiDinato, A. R. and Morris,  A.   H.  Algorithm 708: Significant */
/*     Digit Computation of the Incomplete  Beta  Function Ratios.  ACM */
/*     Trans. Math.  Softw. 18 (1993), 360-373. */

/*     Computation of other parameters involve a search for a value that */
/*     produces  the desired  value  of P.   The search relies  on  the */
/*     monotinicity of P with the other parameter. */


/*                              Note */


/*     The beta density is proportional to */
/*               t^(A-1) * (1-t)^(B-1) */

/* ********************************************************************** */
/*     .. Parameters .. */
/*     .. */
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
    if (! (*which < 1 || *which > 4)) {
	goto L30;
    }
    if (! (*which < 1)) {
	goto L10;
    }
    *bound = 1.;
    goto L20;
L10:
    *bound = 4.;
L20:
    *status = -1;
    return 0;
L30:
    if (*which == 1) {
	goto L70;
    }
    if (! (*p < 0. || *p > 1.)) {
	goto L60;
    }
    if (! (*p < 0.)) {
	goto L40;
    }
    *bound = 0.;
    goto L50;
L40:
    *bound = 1.;
L50:
    *status = -2;
    return 0;
L60:
L70:
    if (*which == 1) {
	goto L110;
    }
    if (! (*q < 0. || *q > 1.)) {
	goto L100;
    }
    if (! (*q < 0.)) {
	goto L80;
    }
    *bound = 0.;
    goto L90;
L80:
    *bound = 1.;
L90:
    *status = -3;
    return 0;
L100:
L110:
    if (*which == 2) {
	goto L150;
    }
    if (! (*x < 0. || *x > 1.)) {
	goto L140;
    }
    if (! (*x < 0.)) {
	goto L120;
    }
    *bound = 0.;
    goto L130;
L120:
    *bound = 1.;
L130:
    *status = -4;
    return 0;
L140:
L150:
    if (*which == 2) {
	goto L190;
    }
    if (! (*y < 0. || *y > 1.)) {
	goto L180;
    }
    if (! (*y < 0.)) {
	goto L160;
    }
    *bound = 0.;
    goto L170;
L160:
    *bound = 1.;
L170:
    *status = -5;
    return 0;
L180:
L190:
    if (*which == 3) {
	goto L210;
    }
    if (! (*a <= 0.)) {
	goto L200;
    }
    *bound = 0.;
    *status = -6;
    return 0;
L200:
L210:
    if (*which == 4) {
	goto L230;
    }
    if (! (*b <= 0.)) {
	goto L220;
    }
    *bound = 0.;
    *status = -7;
    return 0;
L220:
L230:
    if (*which == 1) {
	goto L270;
    }
    pq = *p + *q;
    if (! ((d__1 = pq - .5 - .5, abs(d__1)) > spmpar_(&c__1) * 3.)) {
	goto L260;
    }
    if (! (pq < 0.)) {
	goto L240;
    }
    *bound = 0.;
    goto L250;
L240:
    *bound = 1.;
L250:
    *status = 3;
    return 0;
L260:
L270:
    if (*which == 2) {
	goto L310;
    }
    xy = *x + *y;
    if (! ((d__1 = xy - .5 - .5, abs(d__1)) > spmpar_(&c__1) * 3.)) {
	goto L300;
    }
    if (! (xy < 0.)) {
	goto L280;
    }
    *bound = 0.;
    goto L290;
L280:
    *bound = 1.;
L290:
    *status = 4;
    return 0;
L300:
L310:
    if (! (*which == 1)) {
	qporq = *p <= *q;
    }
    if (1 == *which) {
	cumbet_(x, y, a, b, p, q);
	*status = 0;
    } else if (2 == *which) {
	dstzr_(&c_b35, &c_b36, &c_b37, &c_b38);
	if (! qporq) {
	    goto L340;
	}
	*status = 0;
	dzror_(status, x, &fx, &xlo, &xhi, &qleft, &qhi);
	*y = 1. - *x;
L320:
	if (! (*status == 1)) {
	    goto L330;
	}
	cumbet_(x, y, a, b, &cum, &ccum);
	fx = cum - *p;
	dzror_(status, x, &fx, &xlo, &xhi, &qleft, &qhi);
	*y = 1. - *x;
	goto L320;
L330:
	goto L370;
L340:
	*status = 0;
	dzror_(status, y, &fx, &xlo, &xhi, &qleft, &qhi);
	*x = 1. - *y;
L350:
	if (! (*status == 1)) {
	    goto L360;
	}
	cumbet_(x, y, a, b, &cum, &ccum);
	fx = ccum - *q;
	dzror_(status, y, &fx, &xlo, &xhi, &qleft, &qhi);
	*x = 1. - *y;
	goto L350;
L360:
L370:
	if (! (*status == -1)) {
	    goto L400;
	}
	if (! qleft) {
	    goto L380;
	}
	*status = 1;
	*bound = 0.;
	goto L390;
L380:
	*status = 2;
	*bound = 1.;
L390:
L400:
	;
    } else if (3 == *which) {
	*a = 5.;
	dstinv_(&c_b48, &c_b49, &c_b50, &c_b50, &c_b52, &c_b37, &c_b38);
	*status = 0;
	dinvr_(status, a, &fx, &qleft, &qhi);
L410:
	if (! (*status == 1)) {
	    goto L440;
	}
	cumbet_(x, y, a, b, &cum, &ccum);
	if (! qporq) {
	    goto L420;
	}
	fx = cum - *p;
	goto L430;
L420:
	fx = ccum - *q;
L430:
	dinvr_(status, a, &fx, &qleft, &qhi);
	goto L410;
L440:
	if (! (*status == -1)) {
	    goto L470;
	}
	if (! qleft) {
	    goto L450;
	}
	*status = 1;
	*bound = 1e-100;
	goto L460;
L450:
	*status = 2;
	*bound = 1e100;
L460:
L470:
	;
    } else if (4 == *which) {
	*b = 5.;
	dstinv_(&c_b48, &c_b49, &c_b50, &c_b50, &c_b52, &c_b37, &c_b38);
	*status = 0;
	dinvr_(status, b, &fx, &qleft, &qhi);
L480:
	if (! (*status == 1)) {
	    goto L510;
	}
	cumbet_(x, y, a, b, &cum, &ccum);
	if (! qporq) {
	    goto L490;
	}
	fx = cum - *p;
	goto L500;
L490:
	fx = ccum - *q;
L500:
	dinvr_(status, b, &fx, &qleft, &qhi);
	goto L480;
L510:
	if (! (*status == -1)) {
	    goto L540;
	}
	if (! qleft) {
	    goto L520;
	}
	*status = 1;
	*bound = 1e-100;
	goto L530;
L520:
	*status = 2;
	*bound = 1e100;
L530:
L540:
	;
    }
    return 0;
} /* cdfbet_ */

