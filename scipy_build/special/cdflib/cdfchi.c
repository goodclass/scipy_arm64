/* cdfchi.f -- translated by f2c (version 20190311).
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
static doublereal c_b25 = 0.;
static doublereal c_b26 = 1e100;
static doublereal c_b27 = .5;
static doublereal c_b29 = 5.;
static doublereal c_b30 = 1e-50;
static doublereal c_b31 = 1e-8;
static doublereal c_b40 = 1e-100;

/* Subroutine */ int cdfchi_(integer *which, doublereal *p, doublereal *q, 
	doublereal *x, doublereal *df, integer *status, doublereal *bound)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static doublereal fx, pq;
    static logical qhi;
    static doublereal cum, ccum, porq;
    static logical qleft;
    extern /* Subroutine */ int dinvr_(integer *, doublereal *, doublereal *, 
	    logical *, logical *);
    static logical qporq;
    extern /* Subroutine */ int cumchi_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    extern doublereal spmpar_(integer *);
    extern /* Subroutine */ int dstinv_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);

/* ********************************************************************** */

/*      SUBROUTINE CDFCHI( WHICH, P, Q, X, DF, STATUS, BOUND ) */
/*               Cumulative Distribution Function */
/*               CHI-Square distribution */


/*                              Function */


/*     Calculates any one parameter of the chi-square */
/*     distribution given values for the others. */


/*                              Arguments */


/*     WHICH --> Integer indicating which of the next three argument */
/*               values is to be calculated from the others. */
/*               Legal range: 1..3 */
/*               iwhich = 1 : Calculate P and Q from X and DF */
/*               iwhich = 2 : Calculate X from P,Q and DF */
/*               iwhich = 3 : Calculate DF from P,Q and X */
/*                    INTEGER WHICH */

/*     P <--> The integral from 0 to X of the chi-square */
/*            distribution. */
/*            Input range: [0, 1]. */
/*                    DOUBLE PRECISION P */

/*     Q <--> 1-P. */
/*            Input range: (0, 1]. */
/*            P + Q = 1.0. */
/*                    DOUBLE PRECISION Q */

/*     X <--> Upper limit of integration of the non-central */
/*            chi-square distribution. */
/*            Input range: [0, +infinity). */
/*            Search range: [0,1E100] */
/*                    DOUBLE PRECISION X */

/*     DF <--> Degrees of freedom of the */
/*             chi-square distribution. */
/*             Input range: (0, +infinity). */
/*             Search range: [ 1E-100, 1E100] */
/*                    DOUBLE PRECISION DF */

/*     STATUS <-- 0 if calculation completed correctly */
/*               -I if input parameter number I is out of range */
/*                1 if answer appears to be lower than lowest */
/*                  search bound */
/*                2 if answer appears to be higher than greatest */
/*                  search bound */
/*                3 if P + Q .ne. 1 */
/*               10 indicates error returned from cumgam.  See */
/*                  references in cdfgam */
/*                    INTEGER STATUS */

/*     BOUND <-- Undefined if STATUS is 0 */

/*               Bound exceeded by parameter number I if STATUS */
/*               is negative. */

/*               Lower search bound if STATUS is 1. */

/*               Upper search bound if STATUS is 2. */


/*                              Method */


/*     Formula    26.4.19   of Abramowitz  and     Stegun, Handbook  of */
/*     Mathematical Functions   (1966) is used   to reduce the chisqure */
/*     distribution to the incomplete distribution. */

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
    if (! (*x < 0.)) {
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
    if (! (*df <= 0.)) {
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
    if (*which == 1) {
	goto L220;
    }
    qporq = *p <= *q;
    if (! qporq) {
	goto L200;
    }
    porq = *p;
    goto L210;
L200:
    porq = *q;
L210:
L220:
    if (1 == *which) {
	*status = 0;
	cumchi_(x, df, p, q);
	if (porq > 1.5) {
	    *status = 10;
	    return 0;
	}
    } else if (2 == *which) {
	*x = 5.;
	dstinv_(&c_b25, &c_b26, &c_b27, &c_b27, &c_b29, &c_b30, &c_b31);
	*status = 0;
	dinvr_(status, x, &fx, &qleft, &qhi);
L230:
	if (! (*status == 1)) {
	    goto L270;
	}
	cumchi_(x, df, &cum, &ccum);
	if (! qporq) {
	    goto L240;
	}
	fx = cum - *p;
	goto L250;
L240:
	fx = ccum - *q;
L250:
	if (! (fx + porq > 1.5)) {
	    goto L260;
	}
	*status = 10;
	return 0;
L260:
	dinvr_(status, x, &fx, &qleft, &qhi);
	goto L230;
L270:
	if (! (*status == -1)) {
	    goto L300;
	}
	if (! qleft) {
	    goto L280;
	}
	*status = 1;
	*bound = 0.;
	goto L290;
L280:
	*status = 2;
	*bound = 1e100;
L290:
L300:
	;
    } else if (3 == *which) {
	*df = 5.;
	dstinv_(&c_b40, &c_b26, &c_b27, &c_b27, &c_b29, &c_b30, &c_b31);
	*status = 0;
	dinvr_(status, df, &fx, &qleft, &qhi);
L310:
	if (! (*status == 1)) {
	    goto L350;
	}
	cumchi_(x, df, &cum, &ccum);
	if (! qporq) {
	    goto L320;
	}
	fx = cum - *p;
	goto L330;
L320:
	fx = ccum - *q;
L330:
	if (! (fx + porq > 1.5)) {
	    goto L340;
	}
	*status = 10;
	return 0;
L340:
	dinvr_(status, df, &fx, &qleft, &qhi);
	goto L310;
L350:
	if (! (*status == -1)) {
	    goto L380;
	}
	if (! qleft) {
	    goto L360;
	}
	*status = 1;
	*bound = 1e-100;
	goto L370;
L360:
	*status = 2;
	*bound = 1e100;
L370:
L380:
	;
    }
    return 0;
} /* cdfchi_ */

