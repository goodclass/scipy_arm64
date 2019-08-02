/* cumnbn.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int cumnbn_(doublereal *s, doublereal *xn, doublereal *pr, 
	doublereal *ompr, doublereal *cum, doublereal *ccum)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    extern /* Subroutine */ int cumbet_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);

/* ********************************************************************** */

/*     SUBROUTINE CUMNNBN(S,XN,PR,OMPR,CUM,CCUM) */
/*                    CUmulative Negative BINomial distribution */


/*                              Function */


/*     Returns the probability that it there will be S or fewer failures */
/*     before there are XN successes, with each binomial trial having */
/*     a probability of success PR. */

/*     Prob(# failures = S | XN successes, PR)  = */
/*                        ( XN + S - 1 ) */
/*                        (            ) * PR^XN * (1-PR)^S */
/*                        (      S     ) */


/*                              Arguments */


/*     S --> The number of failures */
/*                                                  S is DOUBLE PRECISION */

/*     XN --> The number of successes */
/*                                                  XN is DOUBLE PRECISIO */

/*     PR --> The probability of success in each binomial trial. */
/*                                                  PR is DOUBLE PRECISIO */

/*     OMPR --> 1 - PR */
/*                                                  OMPR is DOUBLE PRECIS */

/*     CUM <-- Cumulative negative binomial distribution. */
/*                                                  CUM is DOUBLE PRECISI */

/*     CCUM <-- Compliment of Cumulative negative binomial distribution. */
/*                                                  CCUM is DOUBLE PRECIS */


/*                              Method */


/*     Formula  26.5.26    of   Abramowitz  and    Stegun,  Handbook   of */
/*     Mathematical   Functions (1966) is   used  to reduce the  negative */
/*     binomial distribution to the cumulative beta distribution. */

/* ********************************************************************** */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */
    d__1 = *s + 1.;
    cumbet_(pr, ompr, xn, &d__1, cum, ccum);
    return 0;
} /* cumnbn_ */

