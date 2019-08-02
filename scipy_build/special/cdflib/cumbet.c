/* cumbet.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int cumbet_(doublereal *x, doublereal *y, doublereal *a, 
	doublereal *b, doublereal *cum, doublereal *ccum)
{
    static integer ierr;
    extern /* Subroutine */ int bratio_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    ;

/* ********************************************************************** */

/*     SUBROUTINE CUMBET(X,Y,A,B,CUM,CCUM) */
/*          Double precision cUMulative incomplete BETa distribution */


/*                              Function */


/*     Calculates the cdf to X of the incomplete beta distribution */
/*     with parameters a and b.  This is the integral from 0 to x */
/*     of (1/B(a,b))*f(t)) where f(t) = t**(a-1) * (1-t)**(b-1) */


/*                              Arguments */


/*     X --> Upper limit of integration. */
/*                                        X is DOUBLE PRECISION */

/*     Y --> 1 - X. */
/*                                        Y is DOUBLE PRECISION */

/*     A --> First parameter of the beta distribution. */
/*                                        A is DOUBLE PRECISION */

/*     B --> Second parameter of the beta distribution. */
/*                                        B is DOUBLE PRECISION */

/*     CUM <-- Cumulative incomplete beta distribution. */
/*                                        CUM is DOUBLE PRECISION */

/*     CCUM <-- Compliment of Cumulative incomplete beta distribution. */
/*                                        CCUM is DOUBLE PRECISION */


/*                              Method */


/*     Calls the routine BRATIO. */

/*                                   References */

/*     Didonato, Armido R. and Morris, Alfred H. Jr. (1992) Algorithim */
/*     708 Significant Digit Computation of the Incomplete Beta Function */
/*     Ratios. ACM ToMS, Vol.18, No. 3, Sept. 1992, 360-373. */

/* ********************************************************************** */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Routines .. */
/*     .. */
/*     .. Executable Statements .. */
    if (! (*x <= 0.)) {
	goto L10;
    }
    *cum = 0.;
    *ccum = 1.;
    return 0;
L10:
    if (! (*y <= 0.)) {
	goto L20;
    }
    *cum = 1.;
    *ccum = 0.;
    return 0;
L20:
    bratio_(a, b, x, y, cum, ccum, &ierr);
/*     Call bratio routine */
    return 0;
} /* cumbet_ */

