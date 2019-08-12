/* cdfnbn.f -- translated by f2c (version 20190311).
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
static doublereal c_b36 = 1e100;
static doublereal c_b37 = .5;
static doublereal c_b39 = 5.;
static doublereal c_b40 = 1e-50;
static doublereal c_b41 = 1e-8;
static doublereal c_b64 = 1.;

/* Subroutine */ int cdfnbn_(integer *which, doublereal *p, doublereal *q, 
	doublereal *s, doublereal *xn, doublereal *pr, doublereal *ompr, 
	integer *status, doublereal *bound)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static doublereal fx, pq;
    static logical qhi;
    static doublereal cum, xhi, xlo, ccum;
    static logical qleft;
    extern /* Subroutine */ int dinvr_(integer *, doublereal *, doublereal *, 
	    logical *, logical *), dzror_(integer *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, logical *, logical *);
    static logical qporq;
    extern /* Subroutine */ int dstzr_(doublereal *, doublereal *, doublereal 
	    *, doublereal *), cumnbn_(doublereal *, doublereal *, doublereal *
	    , doublereal *, doublereal *, doublereal *);
    extern doublereal spmpar_(integer *);
    extern /* Subroutine */ int dstinv_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal prompr;

/* ********************************************************************** */

/*      SUBROUTINE CDFNBN ( WHICH, P, S, XN, PR, STATUS, BOUND ) */
/*               Cumulative Distribution Function */
/*               Negative BiNomial distribution */


/*                              Function */


/*     Calculates any one parameter of the negative binomial */
/*     distribution given values for the others. */

/*     The  cumulative  negative   binomial  distribution  returns  the */
/*     probability that there  will be  F or fewer failures before  the */
/*     XNth success in binomial trials each of which has probability of */
/*     success PR. */

/*     The individual term of the negative binomial is the probability of */
/*     S failures before XN successes and is */
/*          Choose( S, XN+S-1 ) * PR^(XN) * (1-PR)^S */


/*                              Arguments */


/*     WHICH --> Integer indicating which of the next four argument */
/*               values is to be calculated from the others. */
/*               Legal range: 1..4 */
/*               iwhich = 1 : Calculate P and Q from S,XN,PR and OMPR */
/*               iwhich = 2 : Calculate S from P,Q,XN,PR and OMPR */
/*               iwhich = 3 : Calculate XN from P,Q,S,PR and OMPR */
/*               iwhich = 4 : Calculate PR and OMPR from P,Q,S and XN */
/*                    INTEGER WHICH */

/*     P <--> The cumulation from 0 to S of the  negative */
/*            binomial distribution. */
/*            Input range: [0,1]. */
/*                    DOUBLE PRECISION P */

/*     Q <--> 1-P. */
/*            Input range: (0, 1]. */
/*            P + Q = 1.0. */
/*                    DOUBLE PRECISION Q */

/*     S <--> The upper limit of cumulation of the binomial distribution. */
/*            There are F or fewer failures before the XNth success. */
/*            Input range: [0, +infinity). */
/*            Search range: [0, 1E100] */
/*                    DOUBLE PRECISION S */

/*     XN  <--> The number of successes. */
/*              Input range: [0, +infinity). */
/*              Search range: [0, 1E100] */
/*                    DOUBLE PRECISION XN */

/*     PR  <--> The probability of success in each binomial trial. */
/*              Input range: [0,1]. */
/*              Search range: [0,1]. */
/*                    DOUBLE PRECISION PR */

/*     OMPR  <--> 1-PR */
/*              Input range: [0,1]. */
/*              Search range: [0,1] */
/*              PR + OMPR = 1.0 */
/*                    DOUBLE PRECISION OMPR */

/*     STATUS <-- 0 if calculation completed correctly */
/*               -I if input parameter number I is out of range */
/*                1 if answer appears to be lower than lowest */
/*                  search bound */
/*                2 if answer appears to be higher than greatest */
/*                  search bound */
/*                3 if P + Q .ne. 1 */
/*                4 if PR + OMPR .ne. 1 */
/*                    INTEGER STATUS */

/*     BOUND <-- Undefined if STATUS is 0 */

/*               Bound exceeded by parameter number I if STATUS */
/*               is negative. */

/*               Lower search bound if STATUS is 1. */

/*               Upper search bound if STATUS is 2. */


/*                              Method */


/*     Formula   26.5.26   of   Abramowitz  and  Stegun,  Handbook   of */
/*     Mathematical Functions (1966) is used  to  reduce calculation of */
/*     the cumulative distribution  function to that of  an  incomplete */
/*     beta. */

/*     Computation of other parameters involve a search for a value that */
/*     produces  the desired  value  of P.   The search relies  on  the */
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
    if (! (*xn < 0.)) {
	goto L140;
    }
    *bound = 0.;
    *status = -5;
    return 0;
L140:
L150:
    if (*which == 4) {
	goto L190;
    }
    if (! (*pr < 0. || *pr > 1.)) {
	goto L180;
    }
    if (! (*pr < 0.)) {
	goto L160;
    }
    *bound = 0.;
    goto L170;
L160:
    *bound = 1.;
L170:
    *status = -6;
    return 0;
L180:
L190:
    if (*which == 4) {
	goto L230;
    }
    if (! (*ompr < 0. || *ompr > 1.)) {
	goto L220;
    }
    if (! (*ompr < 0.)) {
	goto L200;
    }
    *bound = 0.;
    goto L210;
L200:
    *bound = 1.;
L210:
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
    if (*which == 4) {
	goto L310;
    }
    prompr = *pr + *ompr;
    if (! ((d__1 = prompr - .5 - .5, abs(d__1)) > spmpar_(&c__1) * 3.)) {
	goto L300;
    }
    if (! (prompr < 0.)) {
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
	cumnbn_(s, xn, pr, ompr, p, q);
	*status = 0;
    } else if (2 == *which) {
	*s = 5.;
	dstinv_(&c_b35, &c_b36, &c_b37, &c_b37, &c_b39, &c_b40, &c_b41);
	*status = 0;
	dinvr_(status, s, &fx, &qleft, &qhi);
L320:
	if (! (*status == 1)) {
	    goto L350;
	}
	cumnbn_(s, xn, pr, ompr, &cum, &ccum);
	if (! qporq) {
	    goto L330;
	}
	fx = cum - *p;
	goto L340;
L330:
	fx = ccum - *q;
L340:
	dinvr_(status, s, &fx, &qleft, &qhi);
	goto L320;
L350:
	if (! (*status == -1)) {
	    goto L380;
	}
	if (! qleft) {
	    goto L360;
	}
	*status = 1;
	*bound = 0.;
	goto L370;
L360:
	*status = 2;
	*bound = 1e100;
L370:
L380:
	;
    } else if (3 == *which) {
	*xn = 5.;
	dstinv_(&c_b35, &c_b36, &c_b37, &c_b37, &c_b39, &c_b40, &c_b41);
	*status = 0;
	dinvr_(status, xn, &fx, &qleft, &qhi);
L390:
	if (! (*status == 1)) {
	    goto L420;
	}
	cumnbn_(s, xn, pr, ompr, &cum, &ccum);
	if (! qporq) {
	    goto L400;
	}
	fx = cum - *p;
	goto L410;
L400:
	fx = ccum - *q;
L410:
	dinvr_(status, xn, &fx, &qleft, &qhi);
	goto L390;
L420:
	if (! (*status == -1)) {
	    goto L450;
	}
	if (! qleft) {
	    goto L430;
	}
	*status = 1;
	*bound = 0.;
	goto L440;
L430:
	*status = 2;
	*bound = 1e100;
L440:
L450:
	;
    } else if (4 == *which) {
	dstzr_(&c_b35, &c_b64, &c_b40, &c_b41);
	if (! qporq) {
	    goto L480;
	}
	*status = 0;
	dzror_(status, pr, &fx, &xlo, &xhi, &qleft, &qhi);
	*ompr = 1. - *pr;
L460:
	if (! (*status == 1)) {
	    goto L470;
	}
	cumnbn_(s, xn, pr, ompr, &cum, &ccum);
	fx = cum - *p;
	dzror_(status, pr, &fx, &xlo, &xhi, &qleft, &qhi);
	*ompr = 1. - *pr;
	goto L460;
L470:
	goto L510;
L480:
	*status = 0;
	dzror_(status, ompr, &fx, &xlo, &xhi, &qleft, &qhi);
	*pr = 1. - *ompr;
L490:
	if (! (*status == 1)) {
	    goto L500;
	}
	cumnbn_(s, xn, pr, ompr, &cum, &ccum);
	fx = ccum - *q;
	dzror_(status, ompr, &fx, &xlo, &xhi, &qleft, &qhi);
	*pr = 1. - *ompr;
	goto L490;
L500:
L510:
	if (! (*status == -1)) {
	    goto L540;
	}
	if (! qleft) {
	    goto L520;
	}
	*status = 1;
	*bound = 0.;
	goto L530;
L520:
	*status = 2;
	*bound = 1.;
L530:
L540:
	;
    }
    return 0;
} /* cdfnbn_ */

