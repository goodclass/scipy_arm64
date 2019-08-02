/* cdffnc.f -- translated by f2c (version 20190311).
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

static doublereal c_b17 = 0.;
static doublereal c_b18 = 1e100;
static doublereal c_b19 = .5;
static doublereal c_b21 = 5.;
static doublereal c_b22 = 1e-50;
static doublereal c_b23 = 1e-8;
static doublereal c_b29 = 1e-100;
static doublereal c_b54 = 1e4;

/* Subroutine */ int cdffnc_(integer *which, doublereal *p, doublereal *q, 
	doublereal *f, doublereal *dfn, doublereal *dfd, doublereal *phonc, 
	integer *status, doublereal *bound)
{
    static doublereal fx;
    static logical qhi;
    static doublereal cum, ccum;
    static logical qleft;
    extern /* Subroutine */ int dinvr_(integer *, doublereal *, doublereal *, 
	    logical *, logical *), cumfnc_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    , dstinv_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static integer errflag;

/* ********************************************************************** */

/*      SUBROUTINE CDFFNC( WHICH, P, Q, F, DFN, DFD, PNONC, STATUS, BOUND ) */
/*               Cumulative Distribution Function */
/*               Non-central F distribution */


/*                              Function */


/*     Calculates any one parameter of the Non-central F */
/*     distribution given values for the others. */


/*                              Arguments */


/*     WHICH --> Integer indicating which of the next five argument */
/*               values is to be calculated from the others. */
/*               Legal range: 1..5 */
/*               iwhich = 1 : Calculate P and Q from F,DFN,DFD and PNONC */
/*               iwhich = 2 : Calculate F from P,Q,DFN,DFD and PNONC */
/*               iwhich = 3 : Calculate DFN from P,Q,F,DFD and PNONC */
/*               iwhich = 4 : Calculate DFD from P,Q,F,DFN and PNONC */
/*               iwhich = 5 : Calculate PNONC from P,Q,F,DFN and DFD */
/*                    INTEGER WHICH */

/*       P <--> The integral from 0 to F of the non-central f-density. */
/*              Input range: [0,1-1E-16). */
/*                    DOUBLE PRECISION P */

/*       Q <--> 1-P. */
/*            Q is not used by this subroutine and is only included */
/*            for similarity with other cdf* routines. */
/*                    DOUBLE PRECISION Q */

/*       F <--> Upper limit of integration of the non-central f-density. */
/*              Input range: [0, +infinity). */
/*              Search range: [0,1E100] */
/*                    DOUBLE PRECISION F */

/*     DFN < --> Degrees of freedom of the numerator sum of squares. */
/*               Input range: (0, +infinity). */
/*               Search range: [ 1E-100, 1E100] */
/*                    DOUBLE PRECISION DFN */

/*     DFD < --> Degrees of freedom of the denominator sum of squares. */
/*               Must be in range: (0, +infinity). */
/*               Input range: (0, +infinity). */
/*               Search range: [ 1E-100, 1E100] */
/*                    DOUBLE PRECISION DFD */

/*     PNONC <-> The non-centrality parameter */
/*               Input range: [0,infinity) */
/*               Search range: [0,1E4] */
/*                    DOUBLE PRECISION PHONC */

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


/*     Formula  26.6.20   of   Abramowitz   and   Stegun,  Handbook  of */
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

/*                              WARNING */

/*     The  value  of the  cumulative  noncentral F distribution is not */
/*     necessarily monotone in either degrees  of freedom.  There  thus */
/*     may be two values that provide a given  CDF value.  This routine */
/*     assumes monotonicity  and will find  an arbitrary one of the two */
/*     values. */

/* ********************************************************************** */
/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
    if (! (*which < 1 || *which > 5)) {
	goto L30;
    }
    if (! (*which < 1)) {
	goto L10;
    }
    *bound = 1.;
    goto L20;
L10:
    *bound = 5.;
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
    if (! (*f < 0.)) {
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
    if (! (*dfn <= 0.)) {
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
    if (! (*dfd <= 0.)) {
	goto L120;
    }
    *bound = 0.;
    *status = -6;
    return 0;
L120:
L130:
    if (*which == 5) {
	goto L150;
    }
    if (! (*phonc < 0.)) {
	goto L140;
    }
    *bound = 0.;
    *status = -7;
    return 0;
L140:
L150:
    if (1 == *which) {
	cumfnc_(f, dfn, dfd, phonc, p, q, &errflag);
	if (errflag != 0) {
	    *status = 10;
	    return 0;
	}
	*status = 0;
    } else if (2 == *which) {
	*f = 5.;
	dstinv_(&c_b17, &c_b18, &c_b19, &c_b19, &c_b21, &c_b22, &c_b23);
	*status = 0;
	dinvr_(status, f, &fx, &qleft, &qhi);
L160:
	if (! (*status == 1)) {
	    goto L170;
	}
	cumfnc_(f, dfn, dfd, phonc, &cum, &ccum, &errflag);
	if (errflag != 0) {
	    *status = 10;
	    return 0;
	}
	fx = cum - *p;
	dinvr_(status, f, &fx, &qleft, &qhi);
	goto L160;
L170:
	if (! (*status == -1)) {
	    goto L200;
	}
	if (! qleft) {
	    goto L180;
	}
	*status = 1;
	*bound = 0.;
	goto L190;
L180:
	*status = 2;
	*bound = 1e100;
L190:
L200:
	;
    } else if (3 == *which) {
	*dfn = 5.;
	dstinv_(&c_b29, &c_b18, &c_b19, &c_b19, &c_b21, &c_b22, &c_b23);
	*status = 0;
	dinvr_(status, dfn, &fx, &qleft, &qhi);
L210:
	if (! (*status == 1)) {
	    goto L220;
	}
	cumfnc_(f, dfn, dfd, phonc, &cum, &ccum, &errflag);
	if (errflag != 0) {
	    *status = 10;
	    return 0;
	}
	fx = cum - *p;
	dinvr_(status, dfn, &fx, &qleft, &qhi);
	goto L210;
L220:
	if (! (*status == -1)) {
	    goto L250;
	}
	if (! qleft) {
	    goto L230;
	}
	*status = 1;
	*bound = 1e-100;
	goto L240;
L230:
	*status = 2;
	*bound = 1e100;
L240:
L250:
	;
    } else if (4 == *which) {
	*dfd = 5.;
	dstinv_(&c_b29, &c_b18, &c_b19, &c_b19, &c_b21, &c_b22, &c_b23);
	*status = 0;
	dinvr_(status, dfd, &fx, &qleft, &qhi);
L260:
	if (! (*status == 1)) {
	    goto L270;
	}
	cumfnc_(f, dfn, dfd, phonc, &cum, &ccum, &errflag);
	if (errflag != 0) {
	    *status = 10;
	    return 0;
	}
	fx = cum - *p;
	dinvr_(status, dfd, &fx, &qleft, &qhi);
	goto L260;
L270:
	if (! (*status == -1)) {
	    goto L300;
	}
	if (! qleft) {
	    goto L280;
	}
	*status = 1;
	*bound = 1e-100;
	goto L290;
L280:
	*status = 2;
	*bound = 1e100;
L290:
L300:
	;
    } else if (5 == *which) {
	*phonc = 5.;
	dstinv_(&c_b17, &c_b54, &c_b19, &c_b19, &c_b21, &c_b22, &c_b23);
	*status = 0;
	dinvr_(status, phonc, &fx, &qleft, &qhi);
L310:
	if (! (*status == 1)) {
	    goto L320;
	}
	cumfnc_(f, dfn, dfd, phonc, &cum, &ccum, &errflag);
	if (errflag != 0) {
	    *status = 10;
	    return 0;
	}
	fx = cum - *p;
	dinvr_(status, phonc, &fx, &qleft, &qhi);
	goto L310;
L320:
	if (! (*status == -1)) {
	    goto L350;
	}
	if (! qleft) {
	    goto L330;
	}
	*status = 1;
	*bound = 0.;
	goto L340;
L330:
	*status = 2;
	*bound = 1e4;
L340:
L350:
	;
    }
    return 0;
} /* cdffnc_ */

