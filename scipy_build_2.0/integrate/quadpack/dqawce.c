/* dqawce.f -- translated by f2c (version 20190311).
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

/* Subroutine */ int dqawce_(D_fp f, doublereal *a, doublereal *b, doublereal 
	*c__, doublereal *epsabs, doublereal *epsrel, integer *limit, 
	doublereal *result, doublereal *abserr, integer *neval, integer *ier, 
	doublereal *alist__, doublereal *blist, doublereal *rlist, doublereal 
	*elist, integer *iord, integer *last)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer k;
    static doublereal a1, a2, b1, b2, aa, bb;
    static integer nev;
    static doublereal area, area1, area2, area12;
    extern /* Subroutine */ int dqc25c_(D_fp, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *);
    static doublereal erro12;
    static integer krule, nrmax;
    static doublereal uflow;
    extern doublereal d1mach_(integer *);
    static integer iroff1, iroff2;
    static doublereal error1, error2, epmach, errbnd, errmax;
    static integer maxerr;
    static doublereal errsum;
    extern /* Subroutine */ int dqpsrt_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *);

/* ***begin prologue  dqawce */
/* ***date written   800101   (yymmdd) */
/* ***revision date  830518   (yymmdd) */
/* ***category no.  h2a2a1,j4 */
/* ***keywords  automatic integrator, special-purpose, */
/*             cauchy principal value, clenshaw-curtis method */
/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
/* ***  purpose  the routine calculates an approximation result to a */
/*              cauchy principal value i = integral of f*w over (a,b) */
/*              (w(x) = 1/(x-c), (c.ne.a, c.ne.b), hopefully satisfying */
/*              following claim for accuracy */
/*              abs(i-result).le.max(epsabs,epsrel*abs(i)) */
/* ***description */

/*        computation of a cauchy principal value */
/*        standard fortran subroutine */
/*        double precision version */

/*        parameters */
/*         on entry */
/*            f      - double precision */
/*                     function subprogram defining the integrand */
/*                     function f(x). the actual name for f needs to be */
/*                     declared e x t e r n a l in the driver program. */

/*            a      - double precision */
/*                     lower limit of integration */

/*            b      - double precision */
/*                     upper limit of integration */

/*            c      - double precision */
/*                     parameter in the weight function, c.ne.a, c.ne.b */
/*                     if c = a or c = b, the routine will end with */
/*                     ier = 6. */

/*            epsabs - double precision */
/*                     absolute accuracy requested */
/*            epsrel - double precision */
/*                     relative accuracy requested */
/*                     if  epsabs.le.0 */
/*                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
/*                     the routine will end with ier = 6. */

/*            limit  - integer */
/*                     gives an upper bound on the number of subintervals */
/*                     in the partition of (a,b), limit.ge.1 */

/*         on return */
/*            result - double precision */
/*                     approximation to the integral */

/*            abserr - double precision */
/*                     estimate of the modulus of the absolute error, */
/*                     which should equal or exceed abs(i-result) */

/*            neval  - integer */
/*                     number of integrand evaluations */

/*            ier    - integer */
/*                     ier = 0 normal and reliable termination of the */
/*                             routine. it is assumed that the requested */
/*                             accuracy has been achieved. */
/*                     ier.gt.0 abnormal termination of the routine */
/*                             the estimates for integral and error are */
/*                             less reliable. it is assumed that the */
/*                             requested accuracy has not been achieved. */
/*            error messages */
/*                     ier = 1 maximum number of subdivisions allowed */
/*                             has been achieved. one can allow more sub- */
/*                             divisions by increasing the value of */
/*                             limit. however, if this yields no */
/*                             improvement it is advised to analyze the */
/*                             the integrand, in order to determine the */
/*                             the integration difficulties. if the */
/*                             position of a local difficulty can be */
/*                             determined (e.g. singularity, */
/*                             discontinuity within the interval) one */
/*                             will probably gain from splitting up the */
/*                             interval at this point and calling */
/*                             appropriate integrators on the subranges. */
/*                         = 2 the occurrence of roundoff error is detec- */
/*                             ted, which prevents the requested */
/*                             tolerance from being achieved. */
/*                         = 3 extremely bad integrand behaviour */
/*                             occurs at some interior points of */
/*                             the integration interval. */
/*                         = 6 the input is invalid, because */
/*                             c = a or c = b or */
/*                             (epsabs.le.0 and */
/*                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28)) */
/*                             or limit.lt.1. */
/*                             result, abserr, neval, rlist(1), elist(1), */
/*                             iord(1) and last are set to zero. alist(1) */
/*                             and blist(1) are set to a and b */
/*                             respectively. */

/*            alist   - double precision */
/*                      vector of dimension at least limit, the first */
/*                       last  elements of which are the left */
/*                      end points of the subintervals in the partition */
/*                      of the given integration range (a,b) */

/*            blist   - double precision */
/*                      vector of dimension at least limit, the first */
/*                       last  elements of which are the right */
/*                      end points of the subintervals in the partition */
/*                      of the given integration range (a,b) */

/*            rlist   - double precision */
/*                      vector of dimension at least limit, the first */
/*                       last  elements of which are the integral */
/*                      approximations on the subintervals */

/*            elist   - double precision */
/*                      vector of dimension limit, the first  last */
/*                      elements of which are the moduli of the absolute */
/*                      error estimates on the subintervals */

/*            iord    - integer */
/*                      vector of dimension at least limit, the first k */
/*                      elements of which are pointers to the error */
/*                      estimates over the subintervals, so that */
/*                      elist(iord(1)), ..., elist(iord(k)) with k = last */
/*                      if last.le.(limit/2+2), and k = limit+1-last */
/*                      otherwise, form a decreasing sequence */

/*            last    - integer */
/*                      number of subintervals actually produced in */
/*                      the subdivision process */

/* ***references  (none) */
/* ***routines called  d1mach,dqc25c,dqpsrt */
/* ***end prologue  dqawce */




/*            list of major variables */
/*            ----------------------- */

/*           alist     - list of left end points of all subintervals */
/*                       considered up to now */
/*           blist     - list of right end points of all subintervals */
/*                       considered up to now */
/*           rlist(i)  - approximation to the integral over */
/*                       (alist(i),blist(i)) */
/*           elist(i)  - error estimate applying to rlist(i) */
/*           maxerr    - pointer to the interval with largest */
/*                       error estimate */
/*           errmax    - elist(maxerr) */
/*           area      - sum of the integrals over the subintervals */
/*           errsum    - sum of the errors over the subintervals */
/*           errbnd    - requested accuracy max(epsabs,epsrel* */
/*                       abs(result)) */
/*           *****1    - variable for the left subinterval */
/*           *****2    - variable for the right subinterval */
/*           last      - index for subdivision */


/*            machine dependent constants */
/*            --------------------------- */

/*           epmach is the largest relative spacing. */
/*           uflow is the smallest positive magnitude. */

/* ***first executable statement  dqawce */
    /* Parameter adjustments */
    --iord;
    --elist;
    --rlist;
    --blist;
    --alist__;

    /* Function Body */
    epmach = d1mach_(&c__4);
    uflow = d1mach_(&c__1);


/*           test on validity of parameters */
/*           ------------------------------ */

    *ier = 6;
    *neval = 0;
    *last = 0;
    alist__[1] = *a;
    blist[1] = *b;
    rlist[1] = 0.;
    elist[1] = 0.;
    iord[1] = 0;
    *result = 0.;
    *abserr = 0.;
/* Computing MAX */
    d__1 = epmach * 50.;
    if (*c__ == *a || *c__ == *b || *epsabs <= 0. && *epsrel < max(d__1,5e-29)
	    ) {
	goto L999;
    }

/*           first approximation to the integral */
/*           ----------------------------------- */

    aa = *a;
    bb = *b;
    if (*a <= *b) {
	goto L10;
    }
    aa = *b;
    bb = *a;
L10:
    *ier = 0;
    krule = 1;
    dqc25c_((D_fp)f, &aa, &bb, c__, result, abserr, &krule, neval);
    *last = 1;
    rlist[1] = *result;
    elist[1] = *abserr;
    iord[1] = 1;
    alist__[1] = *a;
    blist[1] = *b;

/*           test on accuracy */

/* Computing MAX */
    d__1 = *epsabs, d__2 = *epsrel * abs(*result);
    errbnd = max(d__1,d__2);
    if (*limit == 1) {
	*ier = 1;
    }
/* Computing MIN */
    d__1 = abs(*result) * .01;
    if (*abserr < min(d__1,errbnd) || *ier == 1) {
	goto L70;
    }

/*           initialization */
/*           -------------- */

    alist__[1] = aa;
    blist[1] = bb;
    rlist[1] = *result;
    errmax = *abserr;
    maxerr = 1;
    area = *result;
    errsum = *abserr;
    nrmax = 1;
    iroff1 = 0;
    iroff2 = 0;

/*           main do-loop */
/*           ------------ */

    i__1 = *limit;
    for (*last = 2; *last <= i__1; ++(*last)) {

/*           bisect the subinterval with nrmax-th largest */
/*           error estimate. */

	a1 = alist__[maxerr];
	b1 = (alist__[maxerr] + blist[maxerr]) * .5;
	b2 = blist[maxerr];
	if (*c__ <= b1 && *c__ > a1) {
	    b1 = (*c__ + b2) * .5;
	}
	if (*c__ > b1 && *c__ < b2) {
	    b1 = (a1 + *c__) * .5;
	}
	a2 = b1;
	krule = 2;
	dqc25c_((D_fp)f, &a1, &b1, c__, &area1, &error1, &krule, &nev);
	*neval += nev;
	dqc25c_((D_fp)f, &a2, &b2, c__, &area2, &error2, &krule, &nev);
	*neval += nev;

/*           improve previous approximations to integral */
/*           and error and test for accuracy. */

	area12 = area1 + area2;
	erro12 = error1 + error2;
	errsum = errsum + erro12 - errmax;
	area = area + area12 - rlist[maxerr];
	if ((d__1 = rlist[maxerr] - area12, abs(d__1)) < abs(area12) * 1e-5 &&
		 erro12 >= errmax * .99 && krule == 0) {
	    ++iroff1;
	}
	if (*last > 10 && erro12 > errmax && krule == 0) {
	    ++iroff2;
	}
	rlist[maxerr] = area1;
	rlist[*last] = area2;
/* Computing MAX */
	d__1 = *epsabs, d__2 = *epsrel * abs(area);
	errbnd = max(d__1,d__2);
	if (errsum <= errbnd) {
	    goto L15;
	}

/*           test for roundoff error and eventually set error flag. */

	if (iroff1 >= 6 && iroff2 > 20) {
	    *ier = 2;
	}

/*           set error flag in the case that number of interval */
/*           bisections exceeds limit. */

	if (*last == *limit) {
	    *ier = 1;
	}

/*           set error flag in the case of bad integrand behaviour */
/*           at a point of the integration range. */

/* Computing MAX */
	d__1 = abs(a1), d__2 = abs(b2);
	if (max(d__1,d__2) <= (epmach * 100. + 1.) * (abs(a2) + uflow * 1e3)) 
		{
	    *ier = 3;
	}

/*           append the newly-created intervals to the list. */

L15:
	if (error2 > error1) {
	    goto L20;
	}
	alist__[*last] = a2;
	blist[maxerr] = b1;
	blist[*last] = b2;
	elist[maxerr] = error1;
	elist[*last] = error2;
	goto L30;
L20:
	alist__[maxerr] = a2;
	alist__[*last] = a1;
	blist[*last] = b1;
	rlist[maxerr] = area2;
	rlist[*last] = area1;
	elist[maxerr] = error2;
	elist[*last] = error1;

/*           call subroutine dqpsrt to maintain the descending ordering */
/*           in the list of error estimates and select the subinterval */
/*           with nrmax-th largest error estimate (to be bisected next). */

L30:
	dqpsrt_(limit, last, &maxerr, &errmax, &elist[1], &iord[1], &nrmax);
/* ***jump out of do-loop */
	if (*ier != 0 || errsum <= errbnd) {
	    goto L50;
	}
/* L40: */
    }

/*           compute final result. */
/*           --------------------- */

L50:
    *result = 0.;
    i__1 = *last;
    for (k = 1; k <= i__1; ++k) {
	*result += rlist[k];
/* L60: */
    }
    *abserr = errsum;
L70:
    if (aa == *b) {
	*result = -(*result);
    }
L999:
    return 0;
} /* dqawce_ */

