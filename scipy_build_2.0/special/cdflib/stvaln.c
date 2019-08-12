/* stvaln.f -- translated by f2c (version 20190311).
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

static integer c__5 = 5;

doublereal stvaln_(doublereal *p)
{
    /* Initialized data */

    static doublereal xnum[5] = { -.322232431088,-1.,-.342242088547,
	    -.0204231210245,-4.53642210148e-5 };
    static doublereal xden[5] = { .099348462606,.588581570495,.531103462366,
	    .10353775285,.0038560700634 };

    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double log(doublereal), sqrt(doublereal);

    /* Local variables */
    static doublereal y, z__, sign;
    extern doublereal devlpl_(doublereal *, integer *, doublereal *);


/* ********************************************************************** */

/*     DOUBLE PRECISION FUNCTION STVALN(P) */
/*                    STarting VALue for Neton-Raphon */
/*                calculation of Normal distribution Inverse */


/*                              Function */


/*     Returns X  such that CUMNOR(X)  =   P,  i.e., the  integral from - */
/*     infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU is P */


/*                              Arguments */


/*     P --> The probability whose normal deviate is sought. */
/*                    P is DOUBLE PRECISION */


/*                              Method */


/*     The  rational   function   on  page 95    of Kennedy  and  Gentle, */
/*     Statistical Computing, Marcel Dekker, NY , 1980. */

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
    if (! (*p <= .5)) {
	goto L10;
    }
    sign = -1.;
    z__ = *p;
    goto L20;
L10:
    sign = 1.;
    z__ = 1. - *p;
L20:
    y = sqrt(log(z__) * -2.);
    ret_val = y + devlpl_(xnum, &c__5, &y) / devlpl_(xden, &c__5, &y);
    ret_val = sign * ret_val;
    return ret_val;
} /* stvaln_ */

