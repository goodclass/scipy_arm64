/* dqk15i.f -- translated by f2c (version 20190311).
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

static integer c__4 = 4;
static integer c__1 = 1;
static doublereal c_b6 = 1.5;

/* Subroutine */ int dqk15i_(D_fp f, doublereal *boun, integer *inf, 
	doublereal *a, doublereal *b, doublereal *result, doublereal *abserr, 
	doublereal *resabs, doublereal *resasc)
{
    /* Initialized data */

    static doublereal wg[8] = { 0.,.129484966168869693270611432679082,0.,
	    .27970539148927666790146777142378,0.,
	    .381830050505118944950369775488975,0.,
	    .417959183673469387755102040816327 };
    static doublereal xgk[8] = { .991455371120812639206854697526329,
	    .949107912342758524526189684047851,
	    .864864423359769072789712788640926,
	    .741531185599394439863864773280788,
	    .58608723546769113029414483825873,
	    .405845151377397166906606412076961,
	    .207784955007898467600689403773245,0. };
    static doublereal wgk[8] = { .02293532201052922496373200805897,
	    .063092092629978553290700663189204,
	    .104790010322250183839876322541518,
	    .140653259715525918745189590510238,
	    .16900472663926790282658342659855,
	    .190350578064785409913256402421014,
	    .204432940075298892414161999234649,
	    .209482141084727828012999174891714 };

    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static integer j;
    static doublereal fc, fv1[7], fv2[7], absc, dinf, resg, resk, fsum, absc1,
	     absc2, fval1, fval2, hlgth, centr, reskh, uflow;
    extern doublereal d1mach_(integer *);
    static doublereal tabsc1, tabsc2, epmach;

/* ***begin prologue  dqk15i */
/* ***date written   800101   (yymmdd) */
/* ***revision date  830518   (yymmdd) */
/* ***category no.  h2a3a2,h2a4a2 */
/* ***keywords  15-point transformed gauss-kronrod rules */
/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
/* ***purpose  the original (infinite integration range is mapped */
/*            onto the interval (0,1) and (a,b) is a part of (0,1). */
/*            it is the purpose to compute */
/*            i = integral of transformed integrand over (a,b), */
/*            j = integral of abs(transformed integrand) over (a,b). */
/* ***description */

/*           integration rule */
/*           standard fortran subroutine */
/*           double precision version */

/*           parameters */
/*            on entry */
/*              f      - double precision */
/*                       function subprogram defining the integrand */
/*                       function f(x). the actual name for f needs to be */
/*                       declared e x t e r n a l in the calling program. */

/*              boun   - double precision */
/*                       finite bound of original integration */
/*                       range (set to zero if inf = +2) */

/*              inf    - integer */
/*                       if inf = -1, the original interval is */
/*                                   (-infinity,bound), */
/*                       if inf = +1, the original interval is */
/*                                   (bound,+infinity), */
/*                       if inf = +2, the original interval is */
/*                                   (-infinity,+infinity) and */
/*                       the integral is computed as the sum of two */
/*                       integrals, one over (-infinity,0) and one over */
/*                       (0,+infinity). */

/*              a      - double precision */
/*                       lower limit for integration over subrange */
/*                       of (0,1) */

/*              b      - double precision */
/*                       upper limit for integration over subrange */
/*                       of (0,1) */

/*            on return */
/*              result - double precision */
/*                       approximation to the integral i */
/*                       result is computed by applying the 15-point */
/*                       kronrod rule(resk) obtained by optimal addition */
/*                       of abscissae to the 7-point gauss rule(resg). */

/*              abserr - double precision */
/*                       estimate of the modulus of the absolute error, */
/*                       which should equal or exceed abs(i-result) */

/*              resabs - double precision */
/*                       approximation to the integral j */

/*              resasc - double precision */
/*                       approximation to the integral of */
/*                       abs((transformed integrand)-i/(b-a)) over (a,b) */

/* ***references  (none) */
/* ***routines called  d1mach */
/* ***end prologue  dqk15i */



/*           the abscissae and weights are supplied for the interval */
/*           (-1,1).  because of symmetry only the positive abscissae and */
/*           their corresponding weights are given. */

/*           xgk    - abscissae of the 15-point kronrod rule */
/*                    xgk(2), xgk(4), ... abscissae of the 7-point */
/*                    gauss rule */
/*                    xgk(1), xgk(3), ...  abscissae which are optimally */
/*                    added to the 7-point gauss rule */

/*           wgk    - weights of the 15-point kronrod rule */

/*           wg     - weights of the 7-point gauss rule, corresponding */
/*                    to the abscissae xgk(2), xgk(4), ... */
/*                    wg(1), wg(3), ... are set to zero. */





/*           list of major variables */
/*           ----------------------- */

/*           centr  - mid point of the interval */
/*           hlgth  - half-length of the interval */
/*           absc*  - abscissa */
/*           tabsc* - transformed abscissa */
/*           fval*  - function value */
/*           resg   - result of the 7-point gauss formula */
/*           resk   - result of the 15-point kronrod formula */
/*           reskh  - approximation to the mean value of the transformed */
/*                    integrand over (a,b), i.e. to i/(b-a) */

/*           machine dependent constants */
/*           --------------------------- */

/*           epmach is the largest relative spacing. */
/*           uflow is the smallest positive magnitude. */

/* ***first executable statement  dqk15i */
    epmach = d1mach_(&c__4);
    uflow = d1mach_(&c__1);
    dinf = (doublereal) min(1,*inf);

    centr = (*a + *b) * .5;
    hlgth = (*b - *a) * .5;
    tabsc1 = *boun + dinf * (1. - centr) / centr;
    fval1 = (*f)(&tabsc1);
    if (*inf == 2) {
	d__1 = -tabsc1;
	fval1 += (*f)(&d__1);
    }
    fc = fval1 / centr / centr;

/*           compute the 15-point kronrod approximation to */
/*           the integral, and estimate the error. */

    resg = wg[7] * fc;
    resk = wgk[7] * fc;
    *resabs = abs(resk);
    for (j = 1; j <= 7; ++j) {
	absc = hlgth * xgk[j - 1];
	absc1 = centr - absc;
	absc2 = centr + absc;
	tabsc1 = *boun + dinf * (1. - absc1) / absc1;
	tabsc2 = *boun + dinf * (1. - absc2) / absc2;
	fval1 = (*f)(&tabsc1);
	fval2 = (*f)(&tabsc2);
	if (*inf == 2) {
	    d__1 = -tabsc1;
	    fval1 += (*f)(&d__1);
	}
	if (*inf == 2) {
	    d__1 = -tabsc2;
	    fval2 += (*f)(&d__1);
	}
	fval1 = fval1 / absc1 / absc1;
	fval2 = fval2 / absc2 / absc2;
	fv1[j - 1] = fval1;
	fv2[j - 1] = fval2;
	fsum = fval1 + fval2;
	resg += wg[j - 1] * fsum;
	resk += wgk[j - 1] * fsum;
	*resabs += wgk[j - 1] * (abs(fval1) + abs(fval2));
/* L10: */
    }
    reskh = resk * .5;
    *resasc = wgk[7] * (d__1 = fc - reskh, abs(d__1));
    for (j = 1; j <= 7; ++j) {
	*resasc += wgk[j - 1] * ((d__1 = fv1[j - 1] - reskh, abs(d__1)) + (
		d__2 = fv2[j - 1] - reskh, abs(d__2)));
/* L20: */
    }
    *result = resk * hlgth;
    *resasc *= hlgth;
    *resabs *= hlgth;
    *abserr = (d__1 = (resk - resg) * hlgth, abs(d__1));
    if (*resasc != 0. && *abserr != 0.) {
/* Computing MIN */
	d__3 = *abserr * 200. / *resasc;
	d__1 = 1., d__2 = pow_dd(&d__3, &c_b6);
	*abserr = *resasc * min(d__1,d__2);
    }
    if (*resabs > uflow / (epmach * 50.)) {
/* Computing MAX */
	d__1 = epmach * 50. * *resabs;
	*abserr = max(d__1,*abserr);
    }
    return 0;
} /* dqk15i_ */

