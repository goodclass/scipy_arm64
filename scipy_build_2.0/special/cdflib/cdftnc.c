/* cdftnc.f -- translated by f2c (version 20190311).
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

static doublereal c_b12 = -1e100;
static doublereal c_b13 = 1e100;
static doublereal c_b14 = .5;
static doublereal c_b16 = 5.;
static doublereal c_b17 = 1e-50;
static doublereal c_b18 = 1e-8;
static doublereal c_b24 = 1e-100;
static doublereal c_b25 = 1e6;
static doublereal c_b36 = -1e6;

/* Subroutine */ int cdftnc_(integer *which, doublereal *p, doublereal *q, 
	doublereal *t, doublereal *df, doublereal *pnonc, integer *status, 
	doublereal *bound)
{
    static doublereal fx;
    static logical qhi;
    static doublereal cum, ccum;
    static logical qleft;
    extern /* Subroutine */ int dinvr_(integer *, doublereal *, doublereal *, 
	    logical *, logical *), cumtnc_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), dstinv_(doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);

/* *********************************************************************** */

/*      SUBROUTINE CDFTNC( WHICH, P, Q, T, DF, PNONC, STATUS, BOUND ) */
/*               Cumulative Distribution Function */
/*                  Non-Central T distribution */

/*                               Function */

/*     Calculates any one parameter of the noncentral t distribution give */
/*     values for the others. */

/*                               Arguments */

/*     WHICH --> Integer indicating which  argument */
/*               values is to be calculated from the others. */
/*               Legal range: 1..3 */
/*               iwhich = 1 : Calculate P and Q from T,DF,PNONC */
/*               iwhich = 2 : Calculate T from P,Q,DF,PNONC */
/*               iwhich = 3 : Calculate DF from P,Q,T */
/*               iwhich = 4 : Calculate PNONC from P,Q,DF,T */
/*                    INTEGER WHICH */

/*        P <--> The integral from -infinity to t of the noncentral t-den */
/*              Input range: (0,1]. */
/*                    DOUBLE PRECISION P */

/*     Q <--> 1-P. */
/*            Input range: (0, 1]. */
/*            P + Q = 1.0. */
/*                    DOUBLE PRECISION Q */

/*        T <--> Upper limit of integration of the noncentral t-density. */
/*               Input range: ( -infinity, +infinity). */
/*               Search range: [ -1E100, 1E100 ] */
/*                    DOUBLE PRECISION T */

/*        DF <--> Degrees of freedom of the noncentral t-distribution. */
/*                Input range: (0 , +infinity). */
/*                Search range: [1e-100, 1E10] */
/*                    DOUBLE PRECISION DF */

/*     PNONC <--> Noncentrality parameter of the noncentral t-distributio */
/*                Input range: [-1e6, 1E6]. */

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

/*                                Method */

/*     Upper tail    of  the  cumulative  noncentral t is calculated usin */
/*     formulae  from page 532  of Johnson, Kotz,  Balakrishnan, Coninuou */
/*     Univariate Distributions, Vol 2, 2nd Edition.  Wiley (1995) */

/*     Computation of other parameters involve a search for a value that */
/*     produces  the desired  value  of P.   The search relies  on  the */
/*     monotinicity of P with the other parameter. */

/* *********************************************************************** */
/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
    if (*t > 1e100) {
	*t = 1e100;
    } else if (*t < -1e100) {
	*t = -1e100;
    }
    if (*df > 1e10) {
	*df = 1e10;
    }
    if (*t != *t) {
	*status = -4;
	return 0;
    }
    if (*which != 4) {
	if (! (*pnonc >= -1e6)) {
	    *status = -6;
	    *bound = -1e6;
	    return 0;
	} else if (! (*pnonc <= 1e6)) {
	    *status = -6;
	    *bound = 1e6;
	    return 0;
	}
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
    if (*which == 3) {
	goto L90;
    }
    if (*df > 0.) {
	goto L80;
    }
    *bound = 0.;
    *status = -5;
    return 0;
L80:
L90:
    if (*which == 4) {
	goto L100;
    }
L100:
    if (1 == *which) {
	cumtnc_(t, df, pnonc, p, q);
	*status = 0;
    } else if (2 == *which) {
	*t = 5.;
	dstinv_(&c_b12, &c_b13, &c_b14, &c_b14, &c_b16, &c_b17, &c_b18);
	*status = 0;
	dinvr_(status, t, &fx, &qleft, &qhi);
L110:
	if (! (*status == 1)) {
	    goto L120;
	}
	cumtnc_(t, df, pnonc, &cum, &ccum);
	fx = cum - *p;
	dinvr_(status, t, &fx, &qleft, &qhi);
	goto L110;
L120:
	if (! (*status == -1)) {
	    goto L150;
	}
	if (! qleft) {
	    goto L130;
	}
	*status = 1;
	*bound = -1e100;
	goto L140;
L130:
	*status = 2;
	*bound = 1e100;
L140:
L150:
	;
    } else if (3 == *which) {
	*df = 5.;
	dstinv_(&c_b24, &c_b25, &c_b14, &c_b14, &c_b16, &c_b17, &c_b18);
	*status = 0;
	dinvr_(status, df, &fx, &qleft, &qhi);
L160:
	if (! (*status == 1)) {
	    goto L170;
	}
	cumtnc_(t, df, pnonc, &cum, &ccum);
	fx = cum - *p;
	dinvr_(status, df, &fx, &qleft, &qhi);
	goto L160;
L170:
	if (! (*status == -1)) {
	    goto L200;
	}
	if (! qleft) {
	    goto L180;
	}
	*status = 1;
	*bound = 1e-100;
	goto L190;
L180:
	*status = 2;
	*bound = 1e100;
L190:
L200:
	;
    } else if (4 == *which) {
	*pnonc = 5.;
	dstinv_(&c_b36, &c_b25, &c_b14, &c_b14, &c_b16, &c_b17, &c_b18);
	*status = 0;
	dinvr_(status, pnonc, &fx, &qleft, &qhi);
L210:
	if (! (*status == 1)) {
	    goto L220;
	}
	cumtnc_(t, df, pnonc, &cum, &ccum);
	fx = cum - *p;
	dinvr_(status, pnonc, &fx, &qleft, &qhi);
	goto L210;
L220:
	if (! (*status == -1)) {
	    goto L250;
	}
	if (! qleft) {
	    goto L230;
	}
	*status = 1;
	*bound = 0.;
	goto L240;
L230:
	*status = 2;
	*bound = 1e6;
L240:
L250:
	;
    }
    return 0;
} /* cdftnc_ */

