/* cdft.f -- translated by f2c (version 20190311).
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
static doublereal c_b20 = -1e100;
static doublereal c_b21 = 1e100;
static doublereal c_b22 = .5;
static doublereal c_b24 = 5.;
static doublereal c_b25 = 1e-50;
static doublereal c_b26 = 1e-8;
static doublereal c_b34 = 1e-100;
static doublereal c_b35 = 1e10;

/* Subroutine */ int cdft_(integer *which, doublereal *p, doublereal *q, 
	doublereal *t, doublereal *df, integer *status, doublereal *bound)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static doublereal fx, pq;
    extern doublereal dt1_(doublereal *, doublereal *, doublereal *);
    static logical qhi;
    static doublereal cum, ccum;
    extern /* Subroutine */ int cumt_(doublereal *, doublereal *, doublereal *
	    , doublereal *);
    static logical qleft;
    extern /* Subroutine */ int dinvr_(integer *, doublereal *, doublereal *, 
	    logical *, logical *);
    static logical qporq;
    extern doublereal spmpar_(integer *);
    extern /* Subroutine */ int dstinv_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);

/* ********************************************************************** */

/*      SUBROUTINE CDFT( WHICH, P, Q, T, DF, STATUS, BOUND ) */
/*               Cumulative Distribution Function */
/*                         T distribution */


/*                              Function */


/*     Calculates any one parameter of the t distribution given */
/*     values for the others. */


/*                              Arguments */


/*     WHICH --> Integer indicating which  argument */
/*               values is to be calculated from the others. */
/*               Legal range: 1..3 */
/*               iwhich = 1 : Calculate P and Q from T and DF */
/*               iwhich = 2 : Calculate T from P,Q and DF */
/*               iwhich = 3 : Calculate DF from P,Q and T */
/*                    INTEGER WHICH */

/*        P <--> The integral from -infinity to t of the t-density. */
/*              Input range: (0,1]. */
/*                    DOUBLE PRECISION P */

/*     Q <--> 1-P. */
/*            Input range: (0, 1]. */
/*            P + Q = 1.0. */
/*                    DOUBLE PRECISION Q */

/*        T <--> Upper limit of integration of the t-density. */
/*               Input range: ( -infinity, +infinity). */
/*               Search range: [ -1E100, 1E100 ] */
/*                    DOUBLE PRECISION T */

/*        DF <--> Degrees of freedom of the t-distribution. */
/*                Input range: (0 , +infinity). */
/*                Search range: [1e-100, 1E10] */
/*                    DOUBLE PRECISION DF */

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


/*     Formula  26.5.27  of   Abramowitz   and  Stegun,   Handbook   of */
/*     Mathematical Functions  (1966) is used to reduce the computation */
/*     of the cumulative distribution function to that of an incomplete */
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
    if (! (*p <= 0. || *p > 1.)) {
	goto L60;
    }
    if (! (*p <= 0.)) {
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
    if (*which == 3) {
	goto L130;
    }
    if (! (*df <= 0.)) {
	goto L120;
    }
    *bound = 0.;
    *status = -5;
    return 0;
L120:
L130:
    if (*which == 1) {
	goto L170;
    }
    pq = *p + *q;
    if (! ((d__1 = pq - .5 - .5, abs(d__1)) > spmpar_(&c__1) * 3.)) {
	goto L160;
    }
    if (! (pq < 0.)) {
	goto L140;
    }
    *bound = 0.;
    goto L150;
L140:
    *bound = 1.;
L150:
    *status = 3;
    return 0;
L160:
L170:
    if (! (*which == 1)) {
	qporq = *p <= *q;
    }
    if (1 == *which) {
	cumt_(t, df, p, q);
	*status = 0;
    } else if (2 == *which) {
	*t = dt1_(p, q, df);
	dstinv_(&c_b20, &c_b21, &c_b22, &c_b22, &c_b24, &c_b25, &c_b26);
	*status = 0;
	dinvr_(status, t, &fx, &qleft, &qhi);
L180:
	if (! (*status == 1)) {
	    goto L210;
	}
	cumt_(t, df, &cum, &ccum);
	if (! qporq) {
	    goto L190;
	}
	fx = cum - *p;
	goto L200;
L190:
	fx = ccum - *q;
L200:
	dinvr_(status, t, &fx, &qleft, &qhi);
	goto L180;
L210:
	if (! (*status == -1)) {
	    goto L240;
	}
	if (! qleft) {
	    goto L220;
	}
	*status = 1;
	*bound = -1e100;
	goto L230;
L220:
	*status = 2;
	*bound = 1e100;
L230:
L240:
	;
    } else if (3 == *which) {
	*df = 5.;
	dstinv_(&c_b34, &c_b35, &c_b22, &c_b22, &c_b24, &c_b25, &c_b26);
	*status = 0;
	dinvr_(status, df, &fx, &qleft, &qhi);
L250:
	if (! (*status == 1)) {
	    goto L280;
	}
	cumt_(t, df, &cum, &ccum);
	if (! qporq) {
	    goto L260;
	}
	fx = cum - *p;
	goto L270;
L260:
	fx = ccum - *q;
L270:
	dinvr_(status, df, &fx, &qleft, &qhi);
	goto L250;
L280:
	if (! (*status == -1)) {
	    goto L310;
	}
	if (! qleft) {
	    goto L290;
	}
	*status = 1;
	*bound = 1e-100;
	goto L300;
L290:
	*status = 2;
	*bound = 1e10;
L300:
L310:
	;
    }
    return 0;
} /* cdft_ */

