/* cumt.f -- translated by f2c (version 20190311).
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

static doublereal c_b2 = .5;

/* Subroutine */ int cumt_(doublereal *t, doublereal *df, doublereal *cum, 
	doublereal *ccum)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static doublereal a, tt, xx, yy, oma, dfptt;
    extern /* Subroutine */ int cumbet_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);

/* ********************************************************************** */

/*     SUBROUTINE CUMT(T,DF,CUM,CCUM) */
/*                    CUMulative T-distribution */


/*                              Function */


/*     Computes the integral from -infinity to T of the t-density. */


/*                              Arguments */


/*     T --> Upper limit of integration of the t-density. */
/*                                                  T is DOUBLE PRECISION */

/*     DF --> Degrees of freedom of the t-distribution. */
/*                                                  DF is DOUBLE PRECISIO */

/*     CUM <-- Cumulative t-distribution. */
/*                                                  CCUM is DOUBLE PRECIS */

/*     CCUM <-- Compliment of Cumulative t-distribution. */
/*                                                  CCUM is DOUBLE PRECIS */


/*                              Method */


/*     Formula 26.5.27   of     Abramowitz  and   Stegun,    Handbook  of */
/*     Mathematical Functions  is   used   to  reduce the  t-distribution */
/*     to an incomplete beta. */

/* ********************************************************************** */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */
    tt = *t * *t;
    dfptt = *df + tt;
    xx = *df / dfptt;
    yy = tt / dfptt;
    d__1 = *df * .5;
    cumbet_(&xx, &yy, &d__1, &c_b2, &a, &oma);
    if (! (*t <= 0.)) {
	goto L10;
    }
    *cum = a * .5;
    *ccum = oma + *cum;
    goto L20;
L10:
    *ccum = a * .5;
    *cum = oma + *ccum;
L20:
    return 0;
} /* cumt_ */

