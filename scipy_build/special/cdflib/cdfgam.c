/* cdfgam.f -- translated by f2c (version 20190311).
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
static doublereal c_b27 = -1.;
static doublereal c_b28 = 1e-100;
static doublereal c_b29 = 1e100;
static doublereal c_b30 = .5;
static doublereal c_b32 = 5.;
static doublereal c_b33 = 1e-50;
static doublereal c_b34 = 1e-8;

/* Subroutine */ int cdfgam_(integer *which, doublereal *p, doublereal *q, 
	doublereal *x, doublereal *shape, doublereal *scale, integer *status, 
	doublereal *bound)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static doublereal fx, pq, xx;
    static logical qhi;
    static doublereal cum, ccum;
    static integer ierr;
    static doublereal porq;
    static logical qleft;
    extern /* Subroutine */ int dinvr_(integer *, doublereal *, doublereal *, 
	    logical *, logical *);
    static logical qporq;
    extern /* Subroutine */ int cumgam_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal xscale;
    extern /* Subroutine */ int gaminv_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *);
    extern doublereal spmpar_(integer *);
    extern /* Subroutine */ int dstinv_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);

/* ********************************************************************** */

/*      SUBROUTINE CDFGAM( WHICH, P, Q, X, SHAPE, SCALE, STATUS, BOUND ) */
/*               Cumulative Distribution Function */
/*                         GAMma Distribution */


/*                              Function */


/*     Calculates any one parameter of the gamma */
/*     distribution given values for the others. */


/*                              Arguments */


/*     WHICH --> Integer indicating which of the next four argument */
/*               values is to be calculated from the others. */
/*               Legal range: 1..4 */
/*               iwhich = 1 : Calculate P and Q from X,SHAPE and SCALE */
/*               iwhich = 2 : Calculate X from P,Q,SHAPE and SCALE */
/*               iwhich = 3 : Calculate SHAPE from P,Q,X and SCALE */
/*               iwhich = 4 : Calculate SCALE from P,Q,X and SHAPE */
/*                    INTEGER WHICH */

/*     P <--> The integral from 0 to X of the gamma density. */
/*            Input range: [0,1]. */
/*                    DOUBLE PRECISION P */

/*     Q <--> 1-P. */
/*            Input range: (0, 1]. */
/*            P + Q = 1.0. */
/*                    DOUBLE PRECISION Q */


/*     X <--> The upper limit of integration of the gamma density. */
/*            Input range: [0, +infinity). */
/*            Search range: [0,1E100] */
/*                    DOUBLE PRECISION X */

/*     SHAPE <--> The shape parameter of the gamma density. */
/*                Input range: (0, +infinity). */
/*                Search range: [1E-100,1E100] */
/*                  DOUBLE PRECISION SHAPE */


/*     SCALE <--> The scale parameter of the gamma density. */
/*                Input range: (0, +infinity). */
/*                Search range: (1E-100,1E100] */
/*                   DOUBLE PRECISION SCALE */

/*     STATUS <-- 0 if calculation completed correctly */
/*               -I if input parameter number I is out of range */
/*                1 if answer appears to be lower than lowest */
/*                  search bound */
/*                2 if answer appears to be higher than greatest */
/*                  search bound */
/*                3 if P + Q .ne. 1 */
/*                10 if the gamma or inverse gamma routine cannot */
/*                   compute the answer.  Usually happens only for */
/*                   X and SHAPE very large (gt 1E10 or more) */
/*                    INTEGER STATUS */

/*     BOUND <-- Undefined if STATUS is 0 */

/*               Bound exceeded by parameter number I if STATUS */
/*               is negative. */

/*               Lower search bound if STATUS is 1. */

/*               Upper search bound if STATUS is 2. */


/*                              Method */


/*     Cumulative distribution function (P) is calculated directly by */
/*     the code associated with: */

/*     DiDinato, A. R. and Morris, A. H. Computation of the  incomplete */
/*     gamma function  ratios  and their  inverse.   ACM  Trans.  Math. */
/*     Softw. 12 (1986), 377-393. */

/*     Computation of other parameters involve a search for a value that */
/*     produces  the desired  value  of P.   The search relies  on  the */
/*     monotinicity of P with the other parameter. */


/*                              Note */



/*     The gamma density is proportional to */
/*       T**(SHAPE - 1) * EXP(- SCALE * T) */


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
    if (! (*shape <= 0.)) {
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
    if (! (*scale <= 0.)) {
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
    if (*which == 1) {
	goto L240;
    }
    qporq = *p <= *q;
    if (! qporq) {
	goto L220;
    }
    porq = *p;
    goto L230;
L220:
    porq = *q;
L230:
L240:
    if (1 == *which) {
	*status = 0;
	xscale = *x * *scale;
	cumgam_(&xscale, shape, p, q);
	if (*p > 1.5) {
	    *status = 10;
	}
    } else if (2 == *which) {
	gaminv_(shape, &xx, &c_b27, p, q, &ierr);
	if ((doublereal) ierr < 0.) {
	    *status = 10;
	    return 0;
	} else {
	    *x = xx / *scale;
	    *status = 0;
	}
    } else if (3 == *which) {
	*shape = 5.;
	xscale = *x * *scale;
	dstinv_(&c_b28, &c_b29, &c_b30, &c_b30, &c_b32, &c_b33, &c_b34);
	*status = 0;
	dinvr_(status, shape, &fx, &qleft, &qhi);
L250:
	if (! (*status == 1)) {
	    goto L290;
	}
	cumgam_(&xscale, shape, &cum, &ccum);
	if (! qporq) {
	    goto L260;
	}
	fx = cum - *p;
	goto L270;
L260:
	fx = ccum - *q;
L270:
	if (! (qporq && cum > 1.5 || ! qporq && ccum > 1.5)) {
	    goto L280;
	}
	*status = 10;
	return 0;
L280:
	dinvr_(status, shape, &fx, &qleft, &qhi);
	goto L250;
L290:
	if (! (*status == -1)) {
	    goto L320;
	}
	if (! qleft) {
	    goto L300;
	}
	*status = 1;
	*bound = 1e-100;
	goto L310;
L300:
	*status = 2;
	*bound = 1e100;
L310:
L320:
	;
    } else if (4 == *which) {
	gaminv_(shape, &xx, &c_b27, p, q, &ierr);
	if ((doublereal) ierr < 0.) {
	    *status = 10;
	    return 0;
	} else {
	    *scale = xx / *x;
	    *status = 0;
	}
    }
    return 0;
} /* cdfgam_ */

