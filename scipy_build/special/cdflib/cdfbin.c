/* cdfbin.f -- translated by f2c (version 20190311).
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
static doublereal c_b37 = 0.;
static doublereal c_b38 = .5;
static doublereal c_b40 = 5.;
static doublereal c_b41 = 1e-50;
static doublereal c_b42 = 1e-8;
static doublereal c_b50 = 1e-100;
static doublereal c_b51 = 1e100;
static doublereal c_b65 = 1.;

/* Subroutine */ int cdfbin_(integer *which, doublereal *p, doublereal *q, 
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
	    *, doublereal *), cumbin_(doublereal *, doublereal *, doublereal *
	    , doublereal *, doublereal *, doublereal *);
    extern doublereal spmpar_(integer *);
    extern /* Subroutine */ int dstinv_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal prompr;

/* ********************************************************************** */

/*      SUBROUTINE CDFBIN ( WHICH, P, Q, S, XN, PR, OMPR, STATUS, BOUND ) */
/*               Cumulative Distribution Function */
/*                         BINomial distribution */


/*                              Function */


/*     Calculates any one parameter of the binomial */
/*     distribution given values for the others. */


/*                              Arguments */


/*     WHICH --> Integer indicating which of the next four argument */
/*               values is to be calculated from the others. */
/*               Legal range: 1..4 */
/*               iwhich = 1 : Calculate P and Q from S,XN,PR and OMPR */
/*               iwhich = 2 : Calculate S from P,Q,XN,PR and OMPR */
/*               iwhich = 3 : Calculate XN from P,Q,S,PR and OMPR */
/*               iwhich = 4 : Calculate PR and OMPR from P,Q,S and XN */
/*                    INTEGER WHICH */

/*     P <--> The cumulation from 0 to S of the binomial distribution. */
/*            (Probablility of S or fewer successes in XN trials each */
/*            with probability of success PR.) */
/*            Input range: [0,1]. */
/*                    DOUBLE PRECISION P */

/*     Q <--> 1-P. */
/*            Input range: [0, 1]. */
/*            P + Q = 1.0. */
/*                    DOUBLE PRECISION Q */

/*     S <--> The number of successes observed. */
/*            Input range: [0, XN] */
/*            Search range: [0, XN] */
/*                    DOUBLE PRECISION S */

/*     XN  <--> The number of binomial trials. */
/*              Input range: (0, +infinity). */
/*              Search range: [1E-100, 1E100] */
/*                    DOUBLE PRECISION XN */

/*     PR  <--> The probability of success in each binomial trial. */
/*              Input range: [0,1]. */
/*              Search range: [0,1] */
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


/*     Formula  26.5.24    of   Abramowitz  and    Stegun,  Handbook   of */
/*     Mathematical   Functions (1966) is   used  to reduce the  binomial */
/*     distribution  to  the  cumulative incomplete    beta distribution. */

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
    if (! (*which < 1 && *which > 4)) {
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
    if (*which == 3) {
	goto L130;
    }
    if (! (*xn <= 0.)) {
	goto L120;
    }
    *bound = 0.;
    *status = -5;
    return 0;
L120:
L130:
    if (*which == 2) {
	goto L170;
    }
    if (! (*s < 0. || *which != 3 && *s > *xn)) {
	goto L160;
    }
    if (! (*s < 0.)) {
	goto L140;
    }
    *bound = 0.;
    goto L150;
L140:
    *bound = *xn;
L150:
    *status = -4;
    return 0;
L160:
L170:
    if (*which == 4) {
	goto L210;
    }
    if (! (*pr < 0. || *pr > 1.)) {
	goto L200;
    }
    if (! (*pr < 0.)) {
	goto L180;
    }
    *bound = 0.;
    goto L190;
L180:
    *bound = 1.;
L190:
    *status = -6;
    return 0;
L200:
L210:
    if (*which == 4) {
	goto L250;
    }
    if (! (*ompr < 0. || *ompr > 1.)) {
	goto L240;
    }
    if (! (*ompr < 0.)) {
	goto L220;
    }
    *bound = 0.;
    goto L230;
L220:
    *bound = 1.;
L230:
    *status = -7;
    return 0;
L240:
L250:
    if (*which == 1) {
	goto L290;
    }
    pq = *p + *q;
    if (! ((d__1 = pq - .5 - .5, abs(d__1)) > spmpar_(&c__1) * 3.)) {
	goto L280;
    }
    if (! (pq < 0.)) {
	goto L260;
    }
    *bound = 0.;
    goto L270;
L260:
    *bound = 1.;
L270:
    *status = 3;
    return 0;
L280:
L290:
    if (*which == 4) {
	goto L330;
    }
    prompr = *pr + *ompr;
    if (! ((d__1 = prompr - .5 - .5, abs(d__1)) > spmpar_(&c__1) * 3.)) {
	goto L320;
    }
    if (! (prompr < 0.)) {
	goto L300;
    }
    *bound = 0.;
    goto L310;
L300:
    *bound = 1.;
L310:
    *status = 4;
    return 0;
L320:
L330:
    if (! (*which == 1)) {
	qporq = *p <= *q;
    }
    if (1 == *which) {
	cumbin_(s, xn, pr, ompr, p, q);
	*status = 0;
    } else if (2 == *which) {
	*s = *xn / 2.;
	dstinv_(&c_b37, xn, &c_b38, &c_b38, &c_b40, &c_b41, &c_b42);
	*status = 0;
	dinvr_(status, s, &fx, &qleft, &qhi);
L340:
	if (! (*status == 1)) {
	    goto L370;
	}
	cumbin_(s, xn, pr, ompr, &cum, &ccum);
	if (! qporq) {
	    goto L350;
	}
	fx = cum - *p;
	goto L360;
L350:
	fx = ccum - *q;
L360:
	dinvr_(status, s, &fx, &qleft, &qhi);
	goto L340;
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
	*bound = *xn;
L390:
L400:
	;
    } else if (3 == *which) {
	*xn = 5.;
	dstinv_(&c_b50, &c_b51, &c_b38, &c_b38, &c_b40, &c_b41, &c_b42);
	*status = 0;
	dinvr_(status, xn, &fx, &qleft, &qhi);
L410:
	if (! (*status == 1)) {
	    goto L440;
	}
	cumbin_(s, xn, pr, ompr, &cum, &ccum);
	if (! qporq) {
	    goto L420;
	}
	fx = cum - *p;
	goto L430;
L420:
	fx = ccum - *q;
L430:
	dinvr_(status, xn, &fx, &qleft, &qhi);
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
	dstzr_(&c_b37, &c_b65, &c_b41, &c_b42);
	if (! qporq) {
	    goto L500;
	}
	*status = 0;
	dzror_(status, pr, &fx, &xlo, &xhi, &qleft, &qhi);
	*ompr = 1. - *pr;
L480:
	if (! (*status == 1)) {
	    goto L490;
	}
	cumbin_(s, xn, pr, ompr, &cum, &ccum);
	fx = cum - *p;
	dzror_(status, pr, &fx, &xlo, &xhi, &qleft, &qhi);
	*ompr = 1. - *pr;
	goto L480;
L490:
	goto L530;
L500:
	*status = 0;
	dzror_(status, ompr, &fx, &xlo, &xhi, &qleft, &qhi);
	*pr = 1. - *ompr;
L510:
	if (! (*status == 1)) {
	    goto L520;
	}
	cumbin_(s, xn, pr, ompr, &cum, &ccum);
	fx = ccum - *q;
	dzror_(status, ompr, &fx, &xlo, &xhi, &qleft, &qhi);
	*pr = 1. - *ompr;
	goto L510;
L520:
L530:
	if (! (*status == -1)) {
	    goto L560;
	}
	if (! qleft) {
	    goto L540;
	}
	*status = 1;
	*bound = 0.;
	goto L550;
L540:
	*status = 2;
	*bound = 1.;
L550:
L560:
	;
    }
    return 0;
} /* cdfbin_ */

