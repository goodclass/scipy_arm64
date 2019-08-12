/* cumbin.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int cumbin_(doublereal *s, doublereal *xn, doublereal *pr, 
	doublereal *ompr, doublereal *cum, doublereal *ccum)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Local variables */
    extern /* Subroutine */ int cumbet_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);

/* ********************************************************************** */

/*     SUBROUTINE CUMBIN(S,XN,PBIN,OMPR,CUM,CCUM) */
/*                    CUmulative BINomial distribution */


/*                              Function */


/*     Returns the probability   of 0  to  S  successes in  XN   binomial */
/*     trials, each of which has a probability of success, PBIN. */


/*                              Arguments */


/*     S --> The upper limit of cumulation of the binomial distribution. */
/*                                                  S is DOUBLE PRECISION */

/*     XN --> The number of binomial trials. */
/*                                                  XN is DOUBLE PRECISIO */

/*     PBIN --> The probability of success in each binomial trial. */
/*                                                  PBIN is DOUBLE PRECIS */

/*     OMPR --> 1 - PBIN */
/*                                                  OMPR is DOUBLE PRECIS */

/*     CUM <-- Cumulative binomial distribution. */
/*                                                  CUM is DOUBLE PRECISI */

/*     CCUM <-- Compliment of Cumulative binomial distribution. */
/*                                                  CCUM is DOUBLE PRECIS */


/*                              Method */


/*     Formula  26.5.24    of   Abramowitz  and    Stegun,  Handbook   of */
/*     Mathematical   Functions (1966) is   used  to reduce the  binomial */
/*     distribution  to  the  cumulative    beta distribution. */

/* ********************************************************************** */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */
    if (! (*s < *xn)) {
	goto L10;
    }
    d__1 = *s + 1.;
    d__2 = *xn - *s;
    cumbet_(pr, ompr, &d__1, &d__2, ccum, cum);
    goto L20;
L10:
    *cum = 1.;
    *ccum = 0.;
L20:
    return 0;
} /* cumbin_ */

