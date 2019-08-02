/* cumf.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int cumf_(doublereal *f, doublereal *dfn, doublereal *dfd, 
	doublereal *cum, doublereal *ccum)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal xx, yy;
    static integer ierr;
    static doublereal prod, dsum;
    extern /* Subroutine */ int bratio_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    ;

/* ********************************************************************** */

/*     SUBROUTINE CUMF(F,DFN,DFD,CUM,CCUM) */
/*                    CUMulative F distribution */


/*                              Function */


/*     Computes  the  integral from  0  to  F of  the f-density  with DFN */
/*     and DFD degrees of freedom. */


/*                              Arguments */


/*     F --> Upper limit of integration of the f-density. */
/*                                                  F is DOUBLE PRECISION */

/*     DFN --> Degrees of freedom of the numerator sum of squares. */
/*                                                  DFN is DOUBLE PRECISI */

/*     DFD --> Degrees of freedom of the denominator sum of squares. */
/*                                                  DFD is DOUBLE PRECISI */

/*     CUM <-- Cumulative f distribution. */
/*                                                  CUM is DOUBLE PRECISI */

/*     CCUM <-- Compliment of Cumulative f distribution. */
/*                                                  CCUM is DOUBLE PRECIS */


/*                              Method */


/*     Formula  26.5.28 of  Abramowitz and   Stegun   is  used to  reduce */
/*     the cumulative F to a cumulative beta distribution. */


/*                              Note */


/*     If F is less than or equal to 0, 0 is returned. */

/* ********************************************************************** */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Parameters .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */
    if (! (*f <= 0.)) {
	goto L10;
    }
    *cum = 0.;
    *ccum = 1.;
    return 0;
L10:
    prod = *dfn * *f;

/*     XX is such that the incomplete beta with parameters */
/*     DFD/2 and DFN/2 evaluated at XX is 1 - CUM or CCUM */

/*     YY is 1 - XX */

/*     Calculate the smaller of XX and YY accurately */

    dsum = *dfd + prod;
    xx = *dfd / dsum;
    if (xx > .5) {
	yy = prod / dsum;
	xx = 1. - yy;
    } else {
	yy = 1. - xx;
    }
    d__1 = *dfd * .5;
    d__2 = *dfn * .5;
    bratio_(&d__1, &d__2, &xx, &yy, ccum, cum, &ierr);
    return 0;
} /* cumf_ */

