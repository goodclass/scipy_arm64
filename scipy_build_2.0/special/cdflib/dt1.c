/* dt1.f -- translated by f2c (version 20190311).
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

doublereal dt1_(doublereal *p, doublereal *q, doublereal *df)
{
    /* Initialized data */

    static doublereal coef[20]	/* was [5][4] */ = { 1.,1.,0.,0.,0.,3.,16.,5.,
	    0.,0.,-15.,17.,19.,3.,0.,-945.,-1920.,1482.,776.,79. };
    static integer ideg[4] = { 2,3,4,5 };
    static doublereal denom[4] = { 4.,96.,384.,92160. };

    /* System generated locals */
    doublereal ret_val, d__1;

    /* Local variables */
    static integer i__;
    static doublereal x, xp, xx, sum, term;
    extern doublereal devlpl_(doublereal *, integer *, doublereal *);
    static doublereal denpow;
    extern doublereal dinvnr_(doublereal *, doublereal *);

/* ********************************************************************** */

/*     DOUBLE PRECISION FUNCTION DT1(P,Q,DF) */
/*     Double precision Initialize Approximation to */
/*           INVerse of the cumulative T distribution */


/*                              Function */


/*     Returns  the  inverse   of  the T   distribution   function, i.e., */
/*     the integral from 0 to INVT of the T density is P. This is an */
/*     initial approximation */


/*                              Arguments */


/*     P --> The p-value whose inverse from the T distribution is */
/*          desired. */
/*                    P is DOUBLE PRECISION */

/*     Q --> 1-P. */
/*                    Q is DOUBLE PRECISION */

/*     DF --> Degrees of freedom of the T distribution. */
/*                    DF is DOUBLE PRECISION */

/* ********************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Data statements .. */
/*     .. */
/*     .. Executable Statements .. */
    x = (d__1 = dinvnr_(p, q), abs(d__1));
    xx = x * x;
    sum = x;
    denpow = 1.;
    for (i__ = 1; i__ <= 4; ++i__) {
	term = devlpl_(&coef[i__ * 5 - 5], &ideg[i__ - 1], &xx) * x;
	denpow *= *df;
	sum += term / (denpow * denom[i__ - 1]);
/* L10: */
    }
    if (! (*p >= .5)) {
	goto L20;
    }
    xp = sum;
    goto L30;
L20:
    xp = -sum;
L30:
    ret_val = xp;
    return ret_val;
} /* dt1_ */

