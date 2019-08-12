/* cdfnor.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int cdfnor_(integer *which, doublereal *p, doublereal *q, 
	doublereal *x, doublereal *mean, doublereal *sd, integer *status, 
	doublereal *bound)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static doublereal z__, pq;
    extern doublereal dinvnr_(doublereal *, doublereal *), spmpar_(integer *);
    extern /* Subroutine */ int cumnor_(doublereal *, doublereal *, 
	    doublereal *);

/* ********************************************************************** */

/*      SUBROUTINE CDFNOR( WHICH, P, Q, X, MEAN, SD, STATUS, BOUND ) */
/*               Cumulative Distribution Function */
/*               NORmal distribution */


/*                              Function */


/*     Calculates any one parameter of the normal */
/*     distribution given values for the others. */


/*                              Arguments */


/*     WHICH  --> Integer indicating  which of the  next  parameter */
/*     values is to be calculated using values  of the others. */
/*     Legal range: 1..4 */
/*               iwhich = 1 : Calculate P and Q from X,MEAN and SD */
/*               iwhich = 2 : Calculate X from P,Q,MEAN and SD */
/*               iwhich = 3 : Calculate MEAN from P,Q,X and SD */
/*               iwhich = 4 : Calculate SD from P,Q,X and MEAN */
/*                    INTEGER WHICH */

/*     P <--> The integral from -infinity to X of the normal density. */
/*            Input range: (0,1]. */
/*                    DOUBLE PRECISION P */

/*     Q <--> 1-P. */
/*            Input range: (0, 1]. */
/*            P + Q = 1.0. */
/*                    DOUBLE PRECISION Q */

/*     X < --> Upper limit of integration of the normal-density. */
/*             Input range: ( -infinity, +infinity) */
/*                    DOUBLE PRECISION X */

/*     MEAN <--> The mean of the normal density. */
/*               Input range: (-infinity, +infinity) */
/*                    DOUBLE PRECISION MEAN */

/*     SD <--> Standard Deviation of the normal density. */
/*             Input range: (0, +infinity). */
/*                    DOUBLE PRECISION SD */

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




/*     A slightly modified version of ANORM from */

/*     Cody, W.D. (1993). "ALGORITHM 715: SPECFUN - A Portabel FORTRAN */
/*     Package of Special Function Routines and Test Drivers" */
/*     acm Transactions on Mathematical Software. 19, 22-32. */

/*     is used to calculate the cumulative standard normal distribution. */

/*     The rational functions from pages  90-95  of Kennedy and Gentle, */
/*     Statistical  Computing,  Marcel  Dekker, NY,  1980 are  used  as */
/*     starting values to Newton's Iterations which compute the inverse */
/*     standard normal.  Therefore no  searches  are necessary for  any */
/*     parameter. */

/*     For X < -15, the asymptotic expansion for the normal is used  as */
/*     the starting value in finding the inverse standard normal. */
/*     This is formula 26.2.12 of Abramowitz and Stegun. */


/*                              Note */


/*      The normal density is proportional to */
/*      exp( - 0.5 * (( X - MEAN)/SD)**2) */


/* ********************************************************************** */
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
    *status = 0;
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
    if (*which == 1) {
	goto L150;
    }
    pq = *p + *q;
    if (! ((d__1 = pq - .5 - .5, abs(d__1)) > spmpar_(&c__1) * 3.)) {
	goto L140;
    }
    if (! (pq < 0.)) {
	goto L120;
    }
    *bound = 0.;
    goto L130;
L120:
    *bound = 1.;
L130:
    *status = 3;
    return 0;
L140:
L150:
    if (*which == 4) {
	goto L170;
    }
    if (! (*sd <= 0.)) {
	goto L160;
    }
    *bound = 0.;
    *status = -6;
    return 0;
L160:
L170:
    if (1 == *which) {
	z__ = (*x - *mean) / *sd;
	cumnor_(&z__, p, q);
    } else if (2 == *which) {
	z__ = dinvnr_(p, q);
	*x = *sd * z__ + *mean;
    } else if (3 == *which) {
	z__ = dinvnr_(p, q);
	*mean = *x - *sd * z__;
    } else if (4 == *which) {
	z__ = dinvnr_(p, q);
	*sd = (*x - *mean) / z__;
    }
    return 0;
} /* cdfnor_ */

