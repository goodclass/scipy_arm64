/* cumtnc.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int cumtnc_(doublereal *t, doublereal *df, doublereal *pnonc,
	 doublereal *cum, doublereal *ccum)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double log(doublereal), d_int(doublereal *), exp(doublereal);

    /* Local variables */
    static doublereal b, d__, e, s, x, t2, bb, xi, ss, tt, lnx, omx, dum1, 
	    dum2, cent;
    static integer ierr;
    static doublereal xlnd, term, xlne;
    extern /* Subroutine */ int cumt_(doublereal *, doublereal *, doublereal *
	    , doublereal *);
    static doublereal twoi, bcent, dcent, ecent;
    extern doublereal gamln_(doublereal *);
    static doublereal scent, lnomx;
    static logical qrevs;
    static doublereal pnonc2, lambda, halfdf, alghdf, bbcent;
    extern /* Subroutine */ int bratio_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    ;
    static doublereal dpnonc, sscent;
    extern /* Subroutine */ int cumnor_(doublereal *, doublereal *, 
	    doublereal *);

/* ********************************************************************** */

/*     SUBROUTINE CUMTNC(T,DF,PNONC,CUM,CCUM) */

/*                 CUMulative Non-Central T-distribution */


/*                              Function */


/*     Computes the integral from -infinity to T of the non-central */
/*     t-density. */


/*                              Arguments */


/*     T --> Upper limit of integration of the non-central t-density. */
/*                                                  T is DOUBLE PRECISION */

/*     DF --> Degrees of freedom of the non-central t-distribution. */
/*                                                  DF is DOUBLE PRECISIO */

/*     PNONC --> Non-centrality parameter of the non-central t distibutio */
/*                                                  PNONC is DOUBLE PRECI */

/*     CUM <-- Cumulative t-distribution. */
/*                                                  CCUM is DOUBLE PRECIS */

/*     CCUM <-- Compliment of Cumulative t-distribution. */
/*                                                  CCUM is DOUBLE PRECIS */


/*                              Method */

/*     Upper tail    of  the  cumulative  noncentral t   using */
/*     formulae from page 532  of Johnson, Kotz,  Balakrishnan, Coninuous */
/*     Univariate Distributions, Vol 2, 2nd Edition.  Wiley (1995) */

/*     This implementation starts the calculation at i = lambda, */
/*     which is near the largest Di.  It then sums forward and backward. */
/* *********************************************************************** */
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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     Case pnonc essentially zero */
    if (abs(*pnonc) <= 1e-10) {
	cumt_(t, df, cum, ccum);
	return 0;
    }
    qrevs = *t < 0.;
    if (qrevs) {
	tt = -(*t);
	dpnonc = -(*pnonc);
    } else {
	tt = *t;
	dpnonc = *pnonc;
    }
    pnonc2 = dpnonc * dpnonc;
    t2 = tt * tt;
    if (abs(tt) <= 1e-10) {
	d__1 = -(*pnonc);
	cumnor_(&d__1, cum, ccum);
	return 0;
    }
    lambda = pnonc2 * .5;
    x = *df / (*df + t2);
    omx = 1. - x;
    lnx = log(x);
    lnomx = log(omx);
    halfdf = *df * .5;
    alghdf = gamln_(&halfdf);
/*     ******************** Case i = lambda */
    cent = d_int(&lambda);
    if (cent < 1.) {
	cent = 1.;
    }
/*     Compute d=T(2i) in log space and offset by exp(-lambda) */
    d__1 = cent + 1.;
    xlnd = cent * log(lambda) - gamln_(&d__1) - lambda;
    dcent = exp(xlnd);
/*     Compute e=t(2i+1) in log space offset by exp(-lambda) */
    d__1 = cent + 1.5;
    xlne = (cent + .5) * log(lambda) - gamln_(&d__1) - lambda;
    ecent = exp(xlne);
    if (dpnonc < 0.) {
	ecent = -ecent;
    }
/*     Compute bcent=B(2*cent) */
    d__1 = cent + .5;
    bratio_(&halfdf, &d__1, &x, &omx, &bcent, &dum1, &ierr);
/*     compute bbcent=B(2*cent+1) */
    d__1 = cent + 1.;
    bratio_(&halfdf, &d__1, &x, &omx, &bbcent, &dum2, &ierr);
/*     Case bcent and bbcent are essentially zero */
/*     Thus t is effectively infinite */
    if (bcent + bbcent < 1e-10) {
	if (qrevs) {
	    *cum = 0.;
	    *ccum = 1.;
	} else {
	    *cum = 1.;
	    *ccum = 0.;
	}
	return 0;
    }
/*     Case bcent and bbcent are essentially one */
/*     Thus t is effectively zero */
    if (dum1 + dum2 < 1e-10) {
	d__1 = -(*pnonc);
	cumnor_(&d__1, cum, ccum);
	return 0;
    }
/*     First term in ccum is D*B + E*BB */
    *ccum = dcent * bcent + ecent * bbcent;
/*     compute s(cent) = B(2*(cent+1)) - B(2*cent)) */
    d__1 = halfdf + cent + .5;
    d__2 = cent + 1.5;
    scent = gamln_(&d__1) - gamln_(&d__2) - alghdf + halfdf * lnx + (cent + 
	    .5) * lnomx;
    scent = exp(scent);
/*     compute ss(cent) = B(2*cent+3) - B(2*cent+1) */
    d__1 = halfdf + cent + 1.;
    d__2 = cent + 2.;
    sscent = gamln_(&d__1) - gamln_(&d__2) - alghdf + halfdf * lnx + (cent + 
	    1.) * lnomx;
    sscent = exp(sscent);
/*     ******************** Sum Forward */
    xi = cent + 1.;
    twoi = xi * 2.;
    d__ = dcent;
    e = ecent;
    b = bcent;
    bb = bbcent;
    s = scent;
    ss = sscent;
L10:
    b += s;
    bb += ss;
    d__ = lambda / xi * d__;
    e = lambda / (xi + .5) * e;
    term = d__ * b + e * bb;
    *ccum += term;
    s = s * omx * (*df + twoi - 1.) / (twoi + 1.);
    ss = ss * omx * (*df + twoi) / (twoi + 2.);
    xi += 1.;
    twoi = xi * 2.;
    if (abs(term) > *ccum * 1e-7) {
	goto L10;
    }
/*     ******************** Sum Backward */
    xi = cent;
    twoi = xi * 2.;
    d__ = dcent;
    e = ecent;
    b = bcent;
    bb = bbcent;
    s = scent * (twoi + 1.) / ((*df + twoi - 1.) * omx);
    ss = sscent * (twoi + 2.) / ((*df + twoi) * omx);
L20:
    b -= s;
    bb -= ss;
    d__ *= xi / lambda;
    e *= (xi + .5) / lambda;
    term = d__ * b + e * bb;
    *ccum += term;
    xi += -1.;
    if (xi < .5) {
	goto L30;
    }
    twoi = xi * 2.;
    s = s * (twoi + 1.) / ((*df + twoi - 1.) * omx);
    ss = ss * (twoi + 2.) / ((*df + twoi) * omx);
    if (abs(term) > *ccum * 1e-7) {
	goto L20;
    }
L30:
    if (qrevs) {
	*cum = *ccum * .5;
	*ccum = 1. - *cum;
    } else {
	*ccum *= .5;
	*cum = 1. - *ccum;
    }
/*     Due to roundoff error the answer may not lie between zero and one */
/*     Force it to do so */
/* Computing MAX */
    d__1 = min(*cum,1.);
    *cum = max(d__1,0.);
/* Computing MAX */
    d__1 = min(*ccum,1.);
    *ccum = max(d__1,0.);
    return 0;
} /* cumtnc_ */

