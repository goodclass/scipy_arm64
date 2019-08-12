/* cumchn.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int cumchn_(doublereal *x, doublereal *df, doublereal *pnonc,
	 doublereal *cum, doublereal *ccum)
{
    /* Initialized data */

    static doublereal eps = 1e-15;
    static doublereal abstol = 1e-300;

    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double log(doublereal), exp(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal wt, adj, sum, dfd2, term, chid2, lfact;
    static integer icent;
    static doublereal pcent, xnonc, pterm;
    extern doublereal alngam_(doublereal *);
    static doublereal centaj;
    extern /* Subroutine */ int cumchi_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal lcntaj, sumadj, centwt, lcntwt;

/* *********************************************************************** */

/*     SUBROUTINE CUMCHN(X,DF,PNONC,CUM,CCUM) */
/*             CUMulative of the Non-central CHi-square distribution */

/*                               Function */

/*     Calculates     the       cumulative      non-central    chi-square */
/*     distribution, i.e.,  the probability   that  a   random   variable */
/*     which    follows  the  non-central chi-square  distribution,  with */
/*     non-centrality  parameter    PNONC  and   continuous  degrees   of */
/*     freedom DF, is less than or equal to X. */

/*                              Arguments */

/*     X       --> Upper limit of integration of the non-central */
/*                 chi-square distribution. */
/*                                                 X is DOUBLE PRECISION */

/*     DF      --> Degrees of freedom of the non-central */
/*                 chi-square distribution. */
/*                                                 DF is DOUBLE PRECISION */

/*     PNONC   --> Non-centrality parameter of the non-central */
/*                 chi-square distribution. */
/*                                                 PNONC is DOUBLE PRECIS */

/*     CUM <-- Cumulative non-central chi-square distribution. */
/*                                                 CUM is DOUBLE PRECISIO */

/*     CCUM <-- Compliment of Cumulative non-central chi-square distribut */
/*                                                 CCUM is DOUBLE PRECISI */


/*                                Method */

/*     Uses  formula  26.4.25   of  Abramowitz  and  Stegun, Handbook  of */
/*     Mathematical    Functions,  US   NBS   (1966)    to calculate  the */
/*     non-central chi-square. */

/*                                Variables */

/*     EPS     --- Convergence criterion.  The sum stops when a */
/*                 term is less than EPS*SUM. */
/*                                                 EPS is DOUBLE PRECISIO */

/*     CCUM <-- Compliment of Cumulative non-central */
/*              chi-square distribution. */
/*                                                 CCUM is DOUBLE PRECISI */

/* *********************************************************************** */


/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Statement Functions .. */
/*     .. */
/*     .. Data statements .. */
/*     .. */
/*     .. Statement Function definitions .. */
/*     .. */

    if (! (*x <= 0.)) {
	goto L10;
    }
    *cum = 0.;
    *ccum = 1.;
    return 0;
L10:
    if (! (*pnonc <= 1e-10)) {
	goto L20;
    }


/*     When non-centrality parameter is (essentially) zero, */
/*     use cumulative chi-square distribution */


    cumchi_(x, df, cum, ccum);
    return 0;
L20:
    xnonc = *pnonc / 2.;
/* *********************************************************************** */

/*     The following code calculates the weight, chi-square, and */
/*     adjustment term for the central term in the infinite series. */
/*     The central term is the one in which the poisson weight is */
/*     greatest.  The adjustment term is the amount that must */
/*     be subtracted from the chi-square to move up two degrees */
/*     of freedom. */

/* *********************************************************************** */
    icent = (integer) xnonc;
    if (icent == 0) {
	icent = 1;
    }
    chid2 = *x / 2.;


/*     Calculate central weight term */


    d__1 = (doublereal) (icent + 1);
    lfact = alngam_(&d__1);
    lcntwt = -xnonc + icent * log(xnonc) - lfact;
    centwt = exp(lcntwt);


/*     Calculate central chi-square */


    d__1 = *df + 2. * (doublereal) icent;
    cumchi_(x, &d__1, &pcent, ccum);


/*     Calculate central adjustment term */


    dfd2 = (*df + 2. * (doublereal) icent) / 2.;
    d__1 = dfd2 + 1.;
    lfact = alngam_(&d__1);
    lcntaj = dfd2 * log(chid2) - chid2 - lfact;
    centaj = exp(lcntaj);
    sum = centwt * pcent;
/* *********************************************************************** */

/*     Sum backwards from the central term towards zero. */
/*     Quit whenever either */
/*     (1) the zero term is reached, or */
/*     (2) the term gets small relative to the sum, or */

/* *********************************************************************** */
    sumadj = 0.;
    adj = centaj;
    wt = centwt;
    i__ = icent;

    goto L40;
L30:
    if (! (sum >= abstol && term >= eps * sum) || i__ == 0) {
	goto L50;
    }
L40:
    dfd2 = (*df + 2. * (doublereal) i__) / 2.;


/*     Adjust chi-square for two fewer degrees of freedom. */
/*     The adjusted value ends up in PTERM. */


    adj = adj * dfd2 / chid2;
    sumadj += adj;
    pterm = pcent + sumadj;


/*     Adjust poisson weight for J decreased by one */


    wt *= i__ / xnonc;
    term = wt * pterm;
    sum += term;
    --i__;
    goto L30;
L50:
    sumadj = centaj;
/* *********************************************************************** */

/*     Now sum forward from the central term towards infinity. */
/*     Quit when either */
/*     (1) the term gets small relative to the sum, or */

/* *********************************************************************** */
    adj = centaj;
    wt = centwt;
    i__ = icent;

    goto L70;
L60:
    if (! (sum >= abstol && term >= eps * sum)) {
	goto L80;
    }


/*     Update weights for next higher J */


L70:
    wt *= xnonc / (i__ + 1);


/*     Calculate PTERM and add term to sum */


    pterm = pcent - sumadj;
    term = wt * pterm;
    sum += term;


/*     Update adjustment term for DF for next iteration */


    ++i__;
    dfd2 = (*df + 2. * (doublereal) i__) / 2.;
    adj = adj * chid2 / dfd2;
    sumadj += adj;
    goto L60;
L80:
    *cum = sum;
    *ccum = .5 - *cum + .5;

    return 0;
} /* cumchn_ */

