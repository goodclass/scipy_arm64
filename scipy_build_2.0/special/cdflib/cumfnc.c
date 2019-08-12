/* cumfnc.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int cumfnc_(doublereal *f, doublereal *dfn, doublereal *dfd, 
	doublereal *pnonc, doublereal *cum, doublereal *ccum, integer *status)
{
    /* Initialized data */

    static doublereal eps = 1e-4;
    static doublereal abstol = 1e-300;

    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double log(doublereal), exp(doublereal);

    /* Local variables */
    static doublereal b;
    static integer i__;
    static doublereal xx, yy, adn, aup, sum;
    extern /* Subroutine */ int cumf_(doublereal *, doublereal *, doublereal *
	    , doublereal *, doublereal *);
    static integer ierr;
    static doublereal prod, dsum, betdn;
    static integer icent;
    static doublereal betup, xnonc, dummy, xmult;
    extern doublereal alngam_(doublereal *), betaln_(doublereal *, doublereal 
	    *);
    extern /* Subroutine */ int bratio_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    ;
    static doublereal dnterm, centwt, upterm;

/* ********************************************************************** */

/*               F -NON- -C-ENTRAL F DISTRIBUTION */



/*                              Function */


/*     COMPUTES NONCENTRAL F DISTRIBUTION WITH DFN AND DFD */
/*     DEGREES OF FREEDOM AND NONCENTRALITY PARAMETER PNONC */


/*                              Arguments */


/*     X --> UPPER LIMIT OF INTEGRATION OF NONCENTRAL F IN EQUATION */

/*     DFN --> DEGREES OF FREEDOM OF NUMERATOR */

/*     DFD -->  DEGREES OF FREEDOM OF DENOMINATOR */

/*     PNONC --> NONCENTRALITY PARAMETER. */

/*     CUM <-- CUMULATIVE NONCENTRAL F DISTRIBUTION */

/*     CCUM <-- COMPLIMENT OF CUMMULATIVE */


/*                              Method */


/*     USES FORMULA 26.6.20 OF REFERENCE FOR INFINITE SERIES. */
/*     SERIES IS CALCULATED BACKWARD AND FORWARD FROM J = LAMBDA/2 */
/*     (THIS IS THE TERM WITH THE LARGEST POISSON WEIGHT) UNTIL */
/*     THE CONVERGENCE CRITERION IS MET. */

/*     FOR SPEED, THE INCOMPLETE BETA FUNCTIONS ARE EVALUATED */
/*     BY FORMULA 26.5.16. */


/*               REFERENCE */


/*     HANDBOOD OF MATHEMATICAL FUNCTIONS */
/*     EDITED BY MILTON ABRAMOWITZ AND IRENE A. STEGUN */
/*     NATIONAL BUREAU OF STANDARDS APPLIED MATEMATICS SERIES - 55 */
/*     MARCH 1965 */
/*     P 947, EQUATIONS 26.6.17, 26.6.18 */


/*                              Note */


/*     THE SUM CONTINUES UNTIL A SUCCEEDING TERM IS LESS THAN EPS */
/*     TIMES THE SUM (OR THE SUM IS LESS THAN 1.0E-20).  EPS IS */
/*     SET TO 1.0E-4 IN A DATA STATEMENT WHICH CAN BE CHANGED. */

/* ********************************************************************** */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Statement Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Parameters .. */
/*     .. */
/*     .. Data statements .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */
/*     .. Executable Statements .. */

    *status = 0;
    if (! (*f <= 0.)) {
	goto L10;
    }
    *cum = 0.;
    *ccum = 1.;
    return 0;
L10:
    if (! (*pnonc < 1e-10)) {
	goto L20;
    }

/*     Handle case in which the non-centrality parameter is */
/*     (essentially) zero. */
    cumf_(f, dfn, dfd, cum, ccum);
    return 0;
L20:
    xnonc = *pnonc / 2.;
/*     Calculate the central term of the poisson weighting factor. */
    icent = (integer) xnonc;
    if (! ((d__1 = xnonc - icent, abs(d__1)) < 1.)) {
	*status = 1;
	return 0;
    }
    if (icent == 0) {
	icent = 1;
    }
/*     Compute central weight term */
    d__1 = (doublereal) (icent + 1);
    centwt = exp(-xnonc + icent * log(xnonc) - alngam_(&d__1));
/*     Compute central incomplete beta term */
/*     Assure that minimum of arg to beta and 1 - arg is computed */
/*          accurately. */
    prod = *dfn * *f;
    dsum = *dfd + prod;
    yy = *dfd / dsum;
    if (yy > .5) {
	xx = prod / dsum;
	yy = 1. - xx;
    } else {
	xx = 1. - yy;
    }
    d__1 = *dfn * .5 + (doublereal) icent;
    d__2 = *dfd * .5;
    bratio_(&d__1, &d__2, &xx, &yy, &betdn, &dummy, &ierr);
    adn = *dfn / 2. + (doublereal) icent;
    aup = adn;
    b = *dfd / 2.;
    betup = betdn;
    sum = centwt * betdn;
/*     Now sum terms backward from icent until convergence or all done */
    xmult = centwt;
    i__ = icent;
    if (adn < 2.) {
	d__1 = adn + b;
	d__2 = adn + 1.;
	dnterm = exp(alngam_(&d__1) - alngam_(&d__2) - alngam_(&b) + adn * 
		log(xx) + b * log(yy));
    } else {
/*         Same expression, but avoid problems for large adn */
	dnterm = exp(-betaln_(&adn, &b) - log(adn) + adn * log(xx) + b * log(
		yy));
    }
L30:
    d__1 = xmult * betdn;
    if (! (sum >= abstol && d__1 >= eps * sum) || i__ <= 0) {
	goto L40;
    }
    xmult *= i__ / xnonc;
    --i__;
    adn += -1;
    dnterm = (adn + 1) / ((adn + b) * xx) * dnterm;
    betdn += dnterm;
    sum += xmult * betdn;
    goto L30;
L40:
    i__ = icent + 1;
/*     Now sum forwards until convergence */
    xmult = centwt;
    if (aup - 1 + b == 0.) {
	upterm = exp(-alngam_(&aup) - alngam_(&b) + (aup - 1) * log(xx) + b * 
		log(yy));
    } else {
	if (aup < 2.) {
	    d__1 = aup - 1 + b;
	    upterm = exp(alngam_(&d__1) - alngam_(&aup) - alngam_(&b) + (aup 
		    - 1) * log(xx) + b * log(yy));
	} else {
/*             Same expression, but avoid problems for large aup */
	    d__1 = aup - 1;
	    upterm = exp(-betaln_(&d__1, &b) - log(aup - 1) + (aup - 1) * log(
		    xx) + b * log(yy));
	}
    }
    goto L60;
L50:
    d__1 = xmult * betup;
    if (! (sum >= abstol && d__1 >= eps * sum)) {
	goto L70;
    }
L60:
    xmult *= xnonc / i__;
    ++i__;
    aup += 1;
    upterm = (aup + b - 2.) * xx / (aup - 1) * upterm;
    betup -= upterm;
    sum += xmult * betup;
    goto L50;
L70:
    *cum = sum;
    *ccum = .5 - *cum + .5;
    return 0;
} /* cumfnc_ */

