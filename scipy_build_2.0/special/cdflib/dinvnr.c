/* dinvnr.f -- translated by f2c (version 20190311).
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

doublereal dinvnr_(doublereal *p, doublereal *q)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal dx, pp, cum, ccum, xcur;
    static logical qporq;
    static doublereal strtx;
    extern /* Subroutine */ int cumnor_(doublereal *, doublereal *, 
	    doublereal *);
    extern doublereal stvaln_(doublereal *);

/* ********************************************************************** */

/*     DOUBLE PRECISION FUNCTION DINVNR(P,Q) */
/*     Double precision NoRmal distribution INVerse */


/*                              Function */


/*     Returns X  such that CUMNOR(X)  =   P,  i.e., the  integral from - */
/*     infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU is P */


/*                              Arguments */


/*     P --> The probability whose normal deviate is sought. */
/*                    P is DOUBLE PRECISION */

/*     Q --> 1-P */
/*                    P is DOUBLE PRECISION */


/*                              Method */


/*     The  rational   function   on  page 95    of Kennedy  and  Gentle, */
/*     Statistical Computing, Marcel Dekker, NY , 1980 is used as a start */
/*     value for the Newton method of finding roots. */


/*                              Note */


/*     If P or Q .lt. machine EPS returns +/- DINVNR(EPS) */

/* ********************************************************************** */
/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Statement Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     FIND MINIMUM OF P AND Q */

    qporq = *p <= *q;
    if (! qporq) {
	goto L10;
    }
    pp = *p;
    goto L20;
L10:
    pp = *q;

/*     INITIALIZATION STEP */

L20:
    strtx = stvaln_(&pp);
    xcur = strtx;

/*     NEWTON ITERATIONS */

    for (i__ = 1; i__ <= 100; ++i__) {
	cumnor_(&xcur, &cum, &ccum);
	dx = (cum - pp) / (.3989422804014326 * exp(-.5 * xcur * xcur));
	xcur -= dx;
	if ((d__1 = dx / xcur, abs(d__1)) < 1e-13) {
	    goto L40;
	}
/* L30: */
    }
    ret_val = strtx;

/*     IF WE GET HERE, NEWTON HAS FAILED */

    if (! qporq) {
	ret_val = -ret_val;
    }
    return ret_val;

/*     IF WE GET HERE, NEWTON HAS SUCCEEDED */

L40:
    ret_val = xcur;
    if (! qporq) {
	ret_val = -ret_val;
    }
    return ret_val;
} /* dinvnr_ */

