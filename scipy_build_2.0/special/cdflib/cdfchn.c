/* cdfchn.f -- translated by f2c (version 20190311).
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

static doublereal c_b15 = 0.;
static doublereal c_b16 = 1e100;
static doublereal c_b17 = .5;
static doublereal c_b19 = 5.;
static doublereal c_b20 = 1e-50;
static doublereal c_b21 = 1e-8;
static doublereal c_b27 = 1e-100;
static doublereal c_b40 = 1e4;

/* Subroutine */ int cdfchn_(integer *which, doublereal *p, doublereal *q, 
	doublereal *x, doublereal *df, doublereal *pnonc, integer *status, 
	doublereal *bound)
{
    static doublereal fx;
    static logical qhi;
    static doublereal cum, ccum;
    static logical qleft;
    extern /* Subroutine */ int dinvr_(integer *, doublereal *, doublereal *, 
	    logical *, logical *), cumchn_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), dstinv_(doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);

/* ********************************************************************** */

/*      SUBROUTINE CDFCHN( WHICH, P, Q, X, DF, PNONC, STATUS, BOUND ) */
/*               Cumulative Distribution Function */
/*               Non-central Chi-Square */


/*                              Function */


/*     Calculates any one parameter of the non-central chi-square */
/*     distribution given values for the others. */


/*                              Arguments */


/*     WHICH --> Integer indicating which of the next three argument */
/*               values is to be calculated from the others. */
/*               Input range: 1..4 */
/*               iwhich = 1 : Calculate P and Q from X and DF */
/*               iwhich = 2 : Calculate X from P,DF and PNONC */
/*               iwhich = 3 : Calculate DF from P,X and PNONC */
/*               iwhich = 3 : Calculate PNONC from P,X and DF */
/*                    INTEGER WHICH */

/*     P <--> The integral from 0 to X of the non-central chi-square */
/*            distribution. */
/*            Input range: [0, 1-1E-16). */
/*                    DOUBLE PRECISION P */

/*     Q <--> 1-P. */
/*            Q is not used by this subroutine and is only included */
/*            for similarity with other cdf* routines. */
/*                    DOUBLE PRECISION Q */

/*     X <--> Upper limit of integration of the non-central */
/*            chi-square distribution. */
/*            Input range: [0, +infinity). */
/*            Search range: [0,1E100] */
/*                    DOUBLE PRECISION X */

/*     DF <--> Degrees of freedom of the non-central */
/*             chi-square distribution. */
/*             Input range: (0, +infinity). */
/*             Search range: [ 1E-100, 1E100] */
/*                    DOUBLE PRECISION DF */

/*     PNONC <--> Non-centrality parameter of the non-central */
/*                chi-square distribution. */
/*                Input range: [0, +infinity). */
/*                Search range: [0,1E4] */
/*                    DOUBLE PRECISION PNONC */

/*     STATUS <-- 0 if calculation completed correctly */
/*               -I if input parameter number I is out of range */
/*                1 if answer appears to be lower than lowest */
/*                  search bound */
/*                2 if answer appears to be higher than greatest */
/*                  search bound */
/*                    INTEGER STATUS */

/*     BOUND <-- Undefined if STATUS is 0 */

/*               Bound exceeded by parameter number I if STATUS */
/*               is negative. */

/*               Lower search bound if STATUS is 1. */

/*               Upper search bound if STATUS is 2. */


/*                              Method */


/*     Formula  26.4.25   of   Abramowitz   and   Stegun,  Handbook  of */
/*     Mathematical  Functions (1966) is used to compute the cumulative */
/*     distribution function. */

/*     Computation of other parameters involve a search for a value that */
/*     produces  the desired  value  of P.   The search relies  on  the */
/*     monotinicity of P with the other parameter. */


/*                            WARNING */

/*     The computation time  required for this  routine is proportional */
/*     to the noncentrality  parameter  (PNONC).  Very large  values of */
/*     this parameter can consume immense  computer resources.  This is */
/*     why the search range is bounded by 10,000. */

/* ********************************************************************** */
/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
    if (*x > 1e100) {
	*x = 1e100;
    }
    if (*df > 1e100) {
	*df = 1e100;
    }
    if (*pnonc > 1e4) {
	*pnonc = 1e4;
    }
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
    if (! (*p < 0. || *p > .99999999999999989)) {
	goto L60;
    }
    if (! (*p < 0.)) {
	goto L40;
    }
    *bound = 0.;
    goto L50;
L40:
    *bound = .99999999999999989;
L50:
    *status = -2;
    return 0;
L60:
L70:
    if (*which == 2) {
	goto L90;
    }
    if (*x >= 0.) {
	goto L80;
    }
    *bound = 0.;
    *status = -4;
    return 0;
L80:
L90:
    if (*which == 3) {
	goto L110;
    }
    if (*df > 0.) {
	goto L100;
    }
    *bound = 0.;
    *status = -5;
    return 0;
L100:
L110:
    if (*which == 4) {
	goto L130;
    }
    if (*pnonc >= 0.) {
	goto L120;
    }
    *bound = 0.;
    *status = -6;
    return 0;
L120:
L130:
    if (1 == *which) {
	cumchn_(x, df, pnonc, p, q);
	*status = 0;
    } else if (2 == *which) {
	*x = 5.;
	dstinv_(&c_b15, &c_b16, &c_b17, &c_b17, &c_b19, &c_b20, &c_b21);
	*status = 0;
	dinvr_(status, x, &fx, &qleft, &qhi);
L140:
	if (! (*status == 1)) {
	    goto L150;
	}
	cumchn_(x, df, pnonc, &cum, &ccum);
	fx = cum - *p;
	dinvr_(status, x, &fx, &qleft, &qhi);
	goto L140;
L150:
	if (! (*status == -1)) {
	    goto L180;
	}
	if (! qleft) {
	    goto L160;
	}
	*status = 1;
	*bound = 0.;
	goto L170;
L160:
	*status = 2;
	*bound = 1e100;
L170:
L180:
	;
    } else if (3 == *which) {
	*df = 5.;
	dstinv_(&c_b27, &c_b16, &c_b17, &c_b17, &c_b19, &c_b20, &c_b21);
	*status = 0;
	dinvr_(status, df, &fx, &qleft, &qhi);
L190:
	if (! (*status == 1)) {
	    goto L200;
	}
	cumchn_(x, df, pnonc, &cum, &ccum);
	fx = cum - *p;
	dinvr_(status, df, &fx, &qleft, &qhi);
	goto L190;
L200:
	if (! (*status == -1)) {
	    goto L230;
	}
	if (! qleft) {
	    goto L210;
	}
	*status = 1;
	*bound = 1e-100;
	goto L220;
L210:
	*status = 2;
	*bound = 1e100;
L220:
L230:
	;
    } else if (4 == *which) {
	*pnonc = 5.;
	dstinv_(&c_b15, &c_b40, &c_b17, &c_b17, &c_b19, &c_b20, &c_b21);
	*status = 0;
	dinvr_(status, pnonc, &fx, &qleft, &qhi);
L240:
	if (! (*status == 1)) {
	    goto L250;
	}
	cumchn_(x, df, pnonc, &cum, &ccum);
	fx = cum - *p;
	dinvr_(status, pnonc, &fx, &qleft, &qhi);
	goto L240;
L250:
	if (! (*status == -1)) {
	    goto L280;
	}
	if (! qleft) {
	    goto L260;
	}
	*status = 1;
	*bound = 1e-100;
	goto L270;
L260:
	*status = 2;
	*bound = 1e4;
L270:
L280:
	;
    }
    return 0;
} /* cdfchn_ */

