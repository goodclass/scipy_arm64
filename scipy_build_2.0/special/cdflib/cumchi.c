/* cumchi.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int cumchi_(doublereal *x, doublereal *df, doublereal *cum, 
	doublereal *ccum)
{
    static doublereal a, xx;
    extern /* Subroutine */ int cumgam_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);

/* ********************************************************************** */

/*     SUBROUTINE FUNCTION CUMCHI(X,DF,CUM,CCUM) */
/*             CUMulative of the CHi-square distribution */


/*                              Function */


/*     Calculates the cumulative chi-square distribution. */


/*                              Arguments */


/*     X       --> Upper limit of integration of the */
/*                 chi-square distribution. */
/*                                                 X is DOUBLE PRECISION */

/*     DF      --> Degrees of freedom of the */
/*                 chi-square distribution. */
/*                                                 DF is DOUBLE PRECISION */

/*     CUM <-- Cumulative chi-square distribution. */
/*                                                 CUM is DOUBLE PRECISIO */

/*     CCUM <-- Compliment of Cumulative chi-square distribution. */
/*                                                 CCUM is DOUBLE PRECISI */


/*                              Method */


/*     Calls incomplete gamma function (CUMGAM) */

/* ********************************************************************** */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */
    a = *df * .5;
    xx = *x * .5;
    cumgam_(&xx, &a, cum, ccum);
    return 0;
} /* cumchi_ */

