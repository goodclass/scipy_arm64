/* dqawo.f -- translated by f2c (version 20190311).
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

static integer c__1 = 1;
static integer c__26 = 26;

/* Subroutine */ int dqawo_(D_fp f, doublereal *a, doublereal *b, doublereal *
	omega, integer *integr, doublereal *epsabs, doublereal *epsrel, 
	doublereal *result, doublereal *abserr, integer *neval, integer *ier, 
	integer *leniw, integer *maxp1, integer *lenw, integer *last, integer 
	*iwork, doublereal *work)
{
    static integer l1, l2, l3, l4, lvl, limit;
    extern /* Subroutine */ int dqawoe_(D_fp, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, doublereal *);
    static integer momcom;
    extern /* Subroutine */ int xerror_(char *, integer *, integer *, integer 
	    *, ftnlen);

/* ***begin prologue  dqawo */
/* ***date written   800101   (yymmdd) */
/* ***revision date  830518   (yymmdd) */
/* ***category no.  h2a2a1 */
/* ***keywords  automatic integrator, special-purpose, */
/*             integrand with oscillatory cos or sin factor, */
/*             clenshaw-curtis method, (end point) singularities, */
/*             extrapolation, globally adaptive */
/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
/* ***purpose  the routine calculates an approximation result to a given */
/*            definite integral i=integral of f(x)*w(x) over (a,b) */
/*            where w(x) = cos(omega*x) */
/*            or w(x) = sin(omega*x), */
/*            hopefully satisfying following claim for accuracy */
/*            abs(i-result).le.max(epsabs,epsrel*abs(i)). */
/* ***description */

/*        computation of oscillatory integrals */
/*        standard fortran subroutine */
/*        double precision version */

/*        parameters */
/*         on entry */
/*            f      - double precision */
/*                     function subprogram defining the function */
/*                     f(x).  the actual name for f needs to be */
/*                     declared e x t e r n a l in the driver program. */

/*            a      - double precision */
/*                     lower limit of integration */

/*            b      - double precision */
/*                     upper limit of integration */

/*            omega  - double precision */
/*                     parameter in the integrand weight function */

/*            integr - integer */
/*                     indicates which of the weight functions is used */
/*                     integr = 1      w(x) = cos(omega*x) */
/*                     integr = 2      w(x) = sin(omega*x) */
/*                     if integr.ne.1.and.integr.ne.2, the routine will */
/*                     end with ier = 6. */

/*            epsabs - double precision */
/*                     absolute accuracy requested */
/*            epsrel - double precision */
/*                     relative accuracy requested */
/*                     if epsabs.le.0 and */
/*                     epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
/*                     the routine will end with ier = 6. */

/*         on return */
/*            result - double precision */
/*                     approximation to the integral */

/*            abserr - double precision */
/*                     estimate of the modulus of the absolute error, */
/*                     which should equal or exceed abs(i-result) */

/*            neval  - integer */
/*                     number of  integrand evaluations */

/*            ier    - integer */
/*                     ier = 0 normal and reliable termination of the */
/*                             routine. it is assumed that the requested */
/*                             accuracy has been achieved. */
/*                   - ier.gt.0 abnormal termination of the routine. */
/*                             the estimates for integral and error are */
/*                             less reliable. it is assumed that the */
/*                             requested accuracy has not been achieved. */
/*            error messages */
/*                     ier = 1 maximum number of subdivisions allowed */
/*                             (= leniw/2) has been achieved. one can */
/*                             allow more subdivisions by increasing the */
/*                             value of leniw (and taking the according */
/*                             dimension adjustments into account). */
/*                             however, if this yields no improvement it */
/*                             is advised to analyze the integrand in */
/*                             order to determine the integration */
/*                             difficulties. if the position of a local */
/*                             difficulty can be determined (e.g. */
/*                             singularity, discontinuity within the */
/*                             interval) one will probably gain from */
/*                             splitting up the interval at this point */
/*                             and calling the integrator on the */
/*                             subranges. if possible, an appropriate */
/*                             special-purpose integrator should be used */
/*                             which is designed for handling the type of */
/*                             difficulty involved. */
/*                         = 2 the occurrence of roundoff error is */
/*                             detected, which prevents the requested */
/*                             tolerance from being achieved. */
/*                             the error may be under-estimated. */
/*                         = 3 extremely bad integrand behaviour occurs */
/*                             at some interior points of the */
/*                             integration interval. */
/*                         = 4 the algorithm does not converge. */
/*                             roundoff error is detected in the */
/*                             extrapolation table. it is presumed that */
/*                             the requested tolerance cannot be achieved */
/*                             due to roundoff in the extrapolation */
/*                             table, and that the returned result is */
/*                             the best which can be obtained. */
/*                         = 5 the integral is probably divergent, or */
/*                             slowly convergent. it must be noted that */
/*                             divergence can occur with any other value */
/*                             of ier. */
/*                         = 6 the input is invalid, because */
/*                             (epsabs.le.0 and */
/*                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28)) */
/*                             or (integr.ne.1 and integr.ne.2), */
/*                             or leniw.lt.2 or maxp1.lt.1 or */
/*                             lenw.lt.leniw*2+maxp1*25. */
/*                             result, abserr, neval, last are set to */
/*                             zero. except when leniw, maxp1 or lenw are */
/*                             invalid, work(limit*2+1), work(limit*3+1), */
/*                             iwork(1), iwork(limit+1) are set to zero, */
/*                             work(1) is set to a and work(limit+1) to */
/*                             b. */

/*         dimensioning parameters */
/*            leniw  - integer */
/*                     dimensioning parameter for iwork. */
/*                     leniw/2 equals the maximum number of subintervals */
/*                     allowed in the partition of the given integration */
/*                     interval (a,b), leniw.ge.2. */
/*                     if leniw.lt.2, the routine will end with ier = 6. */

/*            maxp1  - integer */
/*                     gives an upper bound on the number of chebyshev */
/*                     moments which can be stored, i.e. for the */
/*                     intervals of lengths abs(b-a)*2**(-l), */
/*                     l=0,1, ..., maxp1-2, maxp1.ge.1 */
/*                     if maxp1.lt.1, the routine will end with ier = 6. */

/*            lenw   - integer */
/*                     dimensioning parameter for work */
/*                     lenw must be at least leniw*2+maxp1*25. */
/*                     if lenw.lt.(leniw*2+maxp1*25), the routine will */
/*                     end with ier = 6. */

/*            last   - integer */
/*                     on return, last equals the number of subintervals */
/*                     produced in the subdivision process, which */
/*                     determines the number of significant elements */
/*                     actually in the work arrays. */

/*         work arrays */
/*            iwork  - integer */
/*                     vector of dimension at least leniw */
/*                     on return, the first k elements of which contain */
/*                     pointers to the error estimates over the */
/*                     subintervals, such that work(limit*3+iwork(1)), .. */
/*                     work(limit*3+iwork(k)) form a decreasing */
/*                     sequence, with limit = lenw/2 , and k = last */
/*                     if last.le.(limit/2+2), and k = limit+1-last */
/*                     otherwise. */
/*                     furthermore, iwork(limit+1), ..., iwork(limit+ */
/*                     last) indicate the subdivision levels of the */
/*                     subintervals, such that iwork(limit+i) = l means */
/*                     that the subinterval numbered i is of length */
/*                     abs(b-a)*2**(1-l). */

/*            work   - double precision */
/*                     vector of dimension at least lenw */
/*                     on return */
/*                     work(1), ..., work(last) contain the left */
/*                      end points of the subintervals in the */
/*                      partition of (a,b), */
/*                     work(limit+1), ..., work(limit+last) contain */
/*                      the right end points, */
/*                     work(limit*2+1), ..., work(limit*2+last) contain */
/*                      the integral approximations over the */
/*                      subintervals, */
/*                     work(limit*3+1), ..., work(limit*3+last) */
/*                      contain the error estimates. */
/*                     work(limit*4+1), ..., work(limit*4+maxp1*25) */
/*                      provide space for storing the chebyshev moments. */
/*                     note that limit = lenw/2. */

/* ***references  (none) */
/* ***routines called  dqawoe,xerror */
/* ***end prologue  dqawo */




/*         check validity of leniw, maxp1 and lenw. */

/* ***first executable statement  dqawo */
    /* Parameter adjustments */
    --iwork;
    --work;

    /* Function Body */
    *ier = 6;
    *neval = 0;
    *last = 0;
    *result = 0.;
    *abserr = 0.;
    if (*leniw < 2 || *maxp1 < 1 || *lenw < (*leniw << 1) + *maxp1 * 25) {
	goto L10;
    }

/*         prepare call for dqawoe */

    limit = *leniw / 2;
    l1 = limit + 1;
    l2 = limit + l1;
    l3 = limit + l2;
    l4 = limit + l3;
    dqawoe_((D_fp)f, a, b, omega, integr, epsabs, epsrel, &limit, &c__1, 
	    maxp1, result, abserr, neval, ier, last, &work[1], &work[l1], &
	    work[l2], &work[l3], &iwork[1], &iwork[l1], &momcom, &work[l4]);

/*         call error handler if necessary */

    lvl = 0;
L10:
    if (*ier == 6) {
	lvl = 0;
    }
    if (*ier != 0) {
	xerror_("abnormal return from dqawo", &c__26, ier, &lvl, (ftnlen)26);
    }
    return 0;
} /* dqawo_ */

