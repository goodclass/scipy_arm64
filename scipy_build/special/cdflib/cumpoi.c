/* cumpoi.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int cumpoi_(doublereal *s, doublereal *xlam, doublereal *cum,
	 doublereal *ccum)
{
    static doublereal df, chi;
    extern /* Subroutine */ int cumchi_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);

/* ********************************************************************** */

/*     SUBROUTINE CUMPOI(S,XLAM,CUM,CCUM) */
/*                    CUMulative POIsson distribution */


/*                              Function */


/*     Returns the  probability  of  S   or  fewer events in  a   Poisson */
/*     distribution with mean XLAM. */


/*                              Arguments */


/*     S --> Upper limit of cumulation of the Poisson. */
/*                                                  S is DOUBLE PRECISION */

/*     XLAM --> Mean of the Poisson distribution. */
/*                                                  XLAM is DOUBLE PRECIS */

/*     CUM <-- Cumulative poisson distribution. */
/*                                        CUM is DOUBLE PRECISION */

/*     CCUM <-- Compliment of Cumulative poisson distribution. */
/*                                                  CCUM is DOUBLE PRECIS */


/*                              Method */


/*     Uses formula  26.4.21   of   Abramowitz and  Stegun,  Handbook  of */
/*     Mathematical   Functions  to reduce   the   cumulative Poisson  to */
/*     the cumulative chi-square distribution. */

/* ********************************************************************** */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Executable Statements .. */
    df = (*s + 1.) * 2.;
    chi = *xlam * 2.;
    cumchi_(&chi, &df, ccum, cum);
    return 0;
} /* cumpoi_ */

