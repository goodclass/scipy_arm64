/* cdfpoi.f -- translated by f2c (version 20190311).
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
static doublereal c_b23 = 0.;
static doublereal c_b24 = 1e100;
static doublereal c_b25 = .5;
static doublereal c_b27 = 5.;
static doublereal c_b28 = 1e-50;
static doublereal c_b29 = 1e-8;

/* Subroutine */ int cdfpoi_(integer *which, doublereal *p, doublereal *q, 
	doublereal *s, doublereal *xlam, integer *status, doublereal *bound)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static doublereal fx, pq;
    static logical qhi;
    static doublereal cum, ccum;
    static logical qleft;
    extern /* Subroutine */ int dinvr_(integer *, doublereal *, doublereal *, 
	    logical *, logical *);
    static logical qporq;
    extern /* Subroutine */ int cumpoi_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    extern doublereal spmpar_(integer *);
    extern /* Subroutine */ int dstinv_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);

/* ********************************************************************** */

/*      SUBROUTINE CDFPOI( WHICH, P, Q, S, XLAM, STATUS, BOUND ) */
/*               Cumulative Distribution Function */
/*               POIsson distribution */


/*                              Function */


/*     Calculates any one parameter of the Poisson */
/*     distribution given values for the others. */


/*                              Arguments */


/*     WHICH --> Integer indicating which  argument */
/*               value is to be calculated from the others. */
/*               Legal range: 1..3 */
/*               iwhich = 1 : Calculate P and Q from S and XLAM */
/*               iwhich = 2 : Calculate S from P,Q and XLAM */
/*               iwhich = 3 : Calculate XLAM from P,Q and S */
/*                    INTEGER WHICH */

/*        P <--> The cumulation from 0 to S of the poisson density. */
/*               Input range: [0,1]. */
/*                    DOUBLE PRECISION P */

/*        Q <--> 1-P. */
/*               Input range: (0, 1]. */
/*               P + Q = 1.0. */
/*                    DOUBLE PRECISION Q */

/*        S <--> Upper limit of cumulation of the Poisson. */
/*               Input range: [0, +infinity). */
/*               Search range: [0,1E100] */
/*                    DOUBLE PRECISION S */

/*     XLAM <--> Mean of the Poisson distribution. */
/*               Input range: [0, +infinity). */
/*               Search range: [0,1E100] */
/*                    DOUBLE PRECISION XLAM */

/*     STATUS <-- 0 if calculation completed correctly */
/*               -I if input parameter number I is out of range */
/*                1 if answer appears to be lower than lowest */
/*                  search bound */
/*                2 if answer appears to be higher than greatest */
/*                  search bound */
/*                3 if P + Q .ne. 1 */
/*                    INTEGER STATUS */

/*     BOUND <-- Undefined if STATUS is 0 */

/*               Bound exceeded by parameter number I if STATUS */
/*               is negative. */

/*               Lower search bound if STATUS is 1. */

/*               Upper search bound if STATUS is 2. */


/*                              Method */


/*     Formula   26.4.21  of   Abramowitz  and   Stegun,   Handbook  of */
/*     Mathematical Functions (1966) is used  to reduce the computation */
/*     of  the cumulative distribution function to that  of computing a */
/*     chi-square, hence an incomplete gamma function. */

/*     Cumulative  distribution function  (P) is  calculated  directly. */
/*     Computation of other parameters involve a search for a value that */
/*     produces  the desired value of  P.   The  search relies  on  the */
/*     monotinicity of P with the other parameter. */


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
    if (! (*which < 1 || *which > 3)) {
	goto L30;
    }
    if (! (*which < 1)) {
	goto L10;
    }
    *bound = 1.;
    goto L20;
L10:
    *bound = 3.;
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
    if (! (*q <= 0. || *q > 1.)) {
	goto L100;
    }
    if (! (*q <= 0.)) {
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
	goto L130;
    }
    if (! (*s < 0.)) {
	goto L120;
    }
    *bound = 0.;
    *status = -4;
    return 0;
L120:
L130:
    if (*which == 3) {
	goto L150;
    }
    if (! (*xlam < 0.)) {
	goto L140;
    }
    *bound = 0.;
    *status = -5;
    return 0;
L140:
L150:
    if (*which == 1) {
	goto L190;
    }
    pq = *p + *q;
    if (! ((d__1 = pq - .5 - .5, abs(d__1)) > spmpar_(&c__1) * 3.)) {
	goto L180;
    }
    if (! (pq < 0.)) {
	goto L160;
    }
    *bound = 0.;
    goto L170;
L160:
    *bound = 1.;
L170:
    *status = 3;
    return 0;
L180:
L190:
    if (! (*which == 1)) {
	qporq = *p <= *q;
    }
    if (1 == *which) {
	cumpoi_(s, xlam, p, q);
	*status = 0;
    } else if (2 == *which) {
	if (*xlam < .01 && *p < .975) {
/*             For sufficiently small xlam and p, the result is 0.0. */
	    *s = 0.;
	    *status = 0;
	    goto L260;
	}
	*s = 5.;
	dstinv_(&c_b23, &c_b24, &c_b25, &c_b25, &c_b27, &c_b28, &c_b29);
	*status = 0;
	dinvr_(status, s, &fx, &qleft, &qhi);
L200:
	if (! (*status == 1)) {
	    goto L230;
	}
	cumpoi_(s, xlam, &cum, &ccum);
	if (! qporq) {
	    goto L210;
	}
	fx = cum - *p;
	goto L220;
L210:
	fx = ccum - *q;
L220:
	dinvr_(status, s, &fx, &qleft, &qhi);
	goto L200;
L230:
	if (! (*status == -1)) {
	    goto L260;
	}
	if (! qleft) {
	    goto L240;
	}
	*status = 1;
	*bound = 0.;
	goto L250;
L240:
	*status = 2;
	*bound = 1e100;
L250:
L260:
	;
    } else if (3 == *which) {
	*xlam = 5.;
	dstinv_(&c_b23, &c_b24, &c_b25, &c_b25, &c_b27, &c_b28, &c_b29);
	*status = 0;
	dinvr_(status, xlam, &fx, &qleft, &qhi);
L270:
	if (! (*status == 1)) {
	    goto L300;
	}
	cumpoi_(s, xlam, &cum, &ccum);
	if (! qporq) {
	    goto L280;
	}
	fx = cum - *p;
	goto L290;
L280:
	fx = ccum - *q;
L290:
	dinvr_(status, xlam, &fx, &qleft, &qhi);
	goto L270;
L300:
	if (! (*status == -1)) {
	    goto L330;
	}
	if (! qleft) {
	    goto L310;
	}
	*status = 1;
	*bound = 0.;
	goto L320;
L310:
	*status = 2;
	*bound = 1e100;
L320:
L330:
	;
    }
    return 0;
} /* cdfpoi_ */

