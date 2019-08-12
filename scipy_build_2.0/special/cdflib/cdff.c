/* cdff.f -- translated by f2c (version 20190311).
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
static doublereal c_b24 = 0.;
static doublereal c_b25 = 1e100;
static doublereal c_b26 = .5;
static doublereal c_b28 = 5.;
static doublereal c_b29 = 1e-50;
static doublereal c_b30 = 1e-8;
static doublereal c_b38 = 1e-100;

/* Subroutine */ int cdff_(integer *which, doublereal *p, doublereal *q, 
	doublereal *f, doublereal *dfn, doublereal *dfd, integer *status, 
	doublereal *bound)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static doublereal fx, pq;
    static logical qhi;
    static doublereal cum, ccum;
    extern /* Subroutine */ int cumf_(doublereal *, doublereal *, doublereal *
	    , doublereal *, doublereal *);
    static logical qleft;
    extern /* Subroutine */ int dinvr_(integer *, doublereal *, doublereal *, 
	    logical *, logical *);
    static logical qporq;
    extern doublereal spmpar_(integer *);
    extern /* Subroutine */ int dstinv_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);

/* ********************************************************************** */

/*      SUBROUTINE CDFF( WHICH, P, Q, F, DFN, DFD, STATUS, BOUND ) */
/*               Cumulative Distribution Function */
/*               F distribution */


/*                              Function */


/*     Calculates any one parameter of the F distribution */
/*     given values for the others. */


/*                              Arguments */


/*     WHICH --> Integer indicating which of the next four argument */
/*               values is to be calculated from the others. */
/*               Legal range: 1..4 */
/*               iwhich = 1 : Calculate P and Q from F,DFN and DFD */
/*               iwhich = 2 : Calculate F from P,Q,DFN and DFD */
/*               iwhich = 3 : Calculate DFN from P,Q,F and DFD */
/*               iwhich = 4 : Calculate DFD from P,Q,F and DFN */
/*                    INTEGER WHICH */

/*       P <--> The integral from 0 to F of the f-density. */
/*              Input range: [0,1]. */
/*                    DOUBLE PRECISION P */

/*       Q <--> 1-P. */
/*              Input range: (0, 1]. */
/*              P + Q = 1.0. */
/*                    DOUBLE PRECISION Q */

/*       F <--> Upper limit of integration of the f-density. */
/*              Input range: [0, +infinity). */
/*              Search range: [0,1E100] */
/*                    DOUBLE PRECISION F */

/*     DFN < --> Degrees of freedom of the numerator sum of squares. */
/*               Input range: (0, +infinity). */
/*               Search range: [ 1E-100, 1E100] */
/*                    DOUBLE PRECISION DFN */

/*     DFD < --> Degrees of freedom of the denominator sum of squares. */
/*               Input range: (0, +infinity). */
/*               Search range: [ 1E-100, 1E100] */
/*                    DOUBLE PRECISION DFD */

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


/*     Formula   26.6.2   of   Abramowitz   and   Stegun,  Handbook  of */
/*     Mathematical  Functions (1966) is used to reduce the computation */
/*     of the  cumulative  distribution function for the  F  variate to */
/*     that of an incomplete beta. */

/*     Computation of other parameters involve a search for a value that */
/*     produces  the desired  value  of P.   The search relies  on  the */
/*     monotinicity of P with the other parameter. */

/*                              WARNING */

/*     The value of the  cumulative  F distribution is  not necessarily */
/*     monotone in  either degrees of freedom.  There  thus may  be two */
/*     values  that  provide a given CDF  value.   This routine assumes */
/*     monotonicity and will find an arbitrary one of the two values. */

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
    if (! (*f < 0.)) {
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
    if (! (*dfn <= 0.)) {
	goto L140;
    }
    *bound = 0.;
    *status = -5;
    return 0;
L140:
L150:
    if (*which == 4) {
	goto L170;
    }
    if (! (*dfd <= 0.)) {
	goto L160;
    }
    *bound = 0.;
    *status = -6;
    return 0;
L160:
L170:
    if (*which == 1) {
	goto L210;
    }
    pq = *p + *q;
    if (! ((d__1 = pq - .5 - .5, abs(d__1)) > spmpar_(&c__1) * 3.)) {
	goto L200;
    }
    if (! (pq < 0.)) {
	goto L180;
    }
    *bound = 0.;
    goto L190;
L180:
    *bound = 1.;
L190:
    *status = 3;
    return 0;
L200:
L210:
    if (! (*which == 1)) {
	qporq = *p <= *q;
    }
    if (1 == *which) {
	cumf_(f, dfn, dfd, p, q);
	*status = 0;
    } else if (2 == *which) {
	*f = 5.;
	dstinv_(&c_b24, &c_b25, &c_b26, &c_b26, &c_b28, &c_b29, &c_b30);
	*status = 0;
	dinvr_(status, f, &fx, &qleft, &qhi);
L220:
	if (! (*status == 1)) {
	    goto L250;
	}
	cumf_(f, dfn, dfd, &cum, &ccum);
	if (! qporq) {
	    goto L230;
	}
	fx = cum - *p;
	goto L240;
L230:
	fx = ccum - *q;
L240:
	dinvr_(status, f, &fx, &qleft, &qhi);
	goto L220;
L250:
	if (! (*status == -1)) {
	    goto L280;
	}
	if (! qleft) {
	    goto L260;
	}
	*status = 1;
	*bound = 0.;
	goto L270;
L260:
	*status = 2;
	*bound = 1e100;
L270:
L280:
	;
    } else if (3 == *which) {
	*dfn = 5.;
	dstinv_(&c_b38, &c_b25, &c_b26, &c_b26, &c_b28, &c_b29, &c_b30);
	*status = 0;
	dinvr_(status, dfn, &fx, &qleft, &qhi);
L290:
	if (! (*status == 1)) {
	    goto L320;
	}
	cumf_(f, dfn, dfd, &cum, &ccum);
	if (! qporq) {
	    goto L300;
	}
	fx = cum - *p;
	goto L310;
L300:
	fx = ccum - *q;
L310:
	dinvr_(status, dfn, &fx, &qleft, &qhi);
	goto L290;
L320:
	if (! (*status == -1)) {
	    goto L350;
	}
	if (! qleft) {
	    goto L330;
	}
	*status = 1;
	*bound = 1e-100;
	goto L340;
L330:
	*status = 2;
	*bound = 1e100;
L340:
L350:
	;
    } else if (4 == *which) {
	*dfd = 5.;
	dstinv_(&c_b38, &c_b25, &c_b26, &c_b26, &c_b28, &c_b29, &c_b30);
	*status = 0;
	dinvr_(status, dfd, &fx, &qleft, &qhi);
L360:
	if (! (*status == 1)) {
	    goto L390;
	}
	cumf_(f, dfn, dfd, &cum, &ccum);
	if (! qporq) {
	    goto L370;
	}
	fx = cum - *p;
	goto L380;
L370:
	fx = ccum - *q;
L380:
	dinvr_(status, dfd, &fx, &qleft, &qhi);
	goto L360;
L390:
	if (! (*status == -1)) {
	    goto L420;
	}
	if (! qleft) {
	    goto L400;
	}
	*status = 1;
	*bound = 1e-100;
	goto L410;
L400:
	*status = 2;
	*bound = 1e100;
L410:
L420:
	;
    }
    return 0;
} /* cdff_ */

