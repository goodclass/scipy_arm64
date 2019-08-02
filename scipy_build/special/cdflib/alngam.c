/* alngam.f -- translated by f2c (version 20190311).
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

static integer c__9 = 9;
static integer c__4 = 4;
static integer c__5 = 5;

doublereal alngam_(doublereal *x)
{
    /* Initialized data */

    static doublereal scoefn[9] = { 62.003838007127258804,
	    36.036772530024836321,20.782472531792126786,6.338067999387272343,
	    2.15994312846059073,.3980671310203570498,.1093115956710439502,
	    .0092381945590275995,.0029737866448101651 };
    static doublereal scoefd[4] = { 62.003838007126989331,
	    9.822521104713994894,-8.906016659497461257,1. };
    static doublereal coef[5] = { .083333333333333023564,
	    -.0027777777768818808,7.9365006754279e-4,-5.94997310889e-4,
	    8.065880899e-4 };

    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static integer i__, n;
    static doublereal xx, prod, offset;
    extern doublereal devlpl_(doublereal *, integer *, doublereal *);

/* ********************************************************************** */

/*     DOUBLE PRECISION FUNCTION ALNGAM(X) */
/*                 double precision LN of the GAMma function */


/*                              Function */


/*     Returns the natural logarithm of GAMMA(X). */


/*                              Arguments */


/*     X --> value at which scaled log gamma is to be returned */
/*                    X is DOUBLE PRECISION */


/*                              Method */


/*     If X .le. 6.0, then use recursion to get X below 3 */
/*     then apply rational approximation number 5236 of */
/*     Hart et al, Computer Approximations, John Wiley and */
/*     Sons, NY, 1968. */

/*     If X .gt. 6.0, then use recursion to get X to at least 12 and */
/*     then use formula 5423 of the same source. */

/* ********************************************************************** */

/*     .. Parameters .. */
/*     .. */
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
    if (! (*x <= 6.)) {
	goto L70;
    }
    prod = 1.;
    xx = *x;
    if (! (*x > 3.)) {
	goto L30;
    }
L10:
    if (! (xx > 3.)) {
	goto L20;
    }
    xx += -1.;
    prod *= xx;
    goto L10;
L20:
L30:
    if (! (*x < 2.)) {
	goto L60;
    }
L40:
    if (! (xx < 2.)) {
	goto L50;
    }
    prod /= xx;
    xx += 1.;
    goto L40;
L50:
L60:
    d__1 = xx - 2.;
    d__2 = xx - 2.;
    ret_val = devlpl_(scoefn, &c__9, &d__1) / devlpl_(scoefd, &c__4, &d__2);


/*     COMPUTE RATIONAL APPROXIMATION TO GAMMA(X) */


    ret_val *= prod;
    ret_val = log(ret_val);
    goto L110;
L70:
    offset = .91893853320467274178;


/*     IF NECESSARY MAKE X AT LEAST 12 AND CARRY CORRECTION IN OFFSET */


    if (*x > 12.) {
	goto L90;
    }
    n = (integer) (12. - *x);
    if (! (n > 0)) {
	goto L90;
    }
    prod = 1.;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	prod *= *x + (doublereal) (i__ - 1);
/* L80: */
    }
    offset -= log(prod);
    xx = *x + (doublereal) n;
    goto L100;
L90:
    xx = *x;


/*     COMPUTE POWER SERIES */


L100:
/* Computing 2nd power */
    d__2 = xx;
    d__1 = 1. / (d__2 * d__2);
    ret_val = devlpl_(coef, &c__5, &d__1) / xx;
    ret_val = ret_val + offset + (xx - .5) * log(xx) - xx;
L110:
    return ret_val;
} /* alngam_ */

