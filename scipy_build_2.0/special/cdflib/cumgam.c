/* cumgam.f -- translated by f2c (version 20190311).
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

static integer c__0 = 0;

/* Subroutine */ int cumgam_(doublereal *x, doublereal *a, doublereal *cum, 
	doublereal *ccum)
{
    extern /* Subroutine */ int gratio_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *);

/* ********************************************************************** */

/*     SUBROUTINE CUMGAM(X,A,CUM,CCUM) */
/*           Double precision cUMulative incomplete GAMma distribution */


/*                              Function */


/*     Computes   the  cumulative        of    the     incomplete   gamma */
/*     distribution, i.e., the integral from 0 to X of */
/*          (1/GAM(A))*EXP(-T)*T**(A-1) DT */
/*     where GAM(A) is the complete gamma function of A, i.e., */
/*          GAM(A) = integral from 0 to infinity of */
/*                    EXP(-T)*T**(A-1) DT */


/*                              Arguments */


/*     X --> The upper limit of integration of the incomplete gamma. */
/*                                                X is DOUBLE PRECISION */

/*     A --> The shape parameter of the incomplete gamma. */
/*                                                A is DOUBLE PRECISION */

/*     CUM <-- Cumulative incomplete gamma distribution. */
/*                                        CUM is DOUBLE PRECISION */

/*     CCUM <-- Compliment of Cumulative incomplete gamma distribution. */
/*                                                CCUM is DOUBLE PRECISIO */


/*                              Method */


/*     Calls the routine GRATIO. */

/* ********************************************************************** */

/*     .. */
/*     .. Scalar Arguments .. */
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
    gratio_(a, x, cum, ccum, &c__0);
/*     Call gratio routine */
    return 0;
} /* cumgam_ */

